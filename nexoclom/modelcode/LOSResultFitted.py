import os

import numpy as np
import pandas as pd
import pickle
import astropy.units as u
import sqlalchemy as sqla
from nexoclom import Output, engine
from nexoclom.modelcode.LOSResult import LOSResult
from nexoclom.modelcode.compute_iteration import IterationResultFitted


xcols = ['x', 'y', 'z']
borecols = ['xbore', 'ybore', 'zbore']


class LOSResultFitted(LOSResult):
    def __init__(self, scdata, label_for_fitted, params=None, dphi=1*u.deg, **kwargs):
        inputs = scdata.model_result[label_for_fitted].inputs
        inputs.options.fitted = True
        super().__init__(scdata, inputs, params=params, dphi=dphi, **kwargs)

        self.unfitted_label = label_for_fitted
        self.unfit_outid = None
        self.unfit_outputfiles = None

    # Helper functions
    @staticmethod
    def _should_add_weight(index, saved):
        return index in saved

    @staticmethod
    def _add_weight(x, ratio):
        return np.append(x, ratio)

    @staticmethod
    def _add_index(x, i):
        return np.append(x, i)
    
    def fitted_iteration_search(self, ufit_id):
        metadata_obj = sqla.MetaData()
        table = sqla.Table("uvvsmodels", metadata_obj, autoload_with=engine)

        query = sqla.select(table).where(
            table.columns.unfit_idnum == ufit_id,
            table.columns.quantity == self.quantity,
            table.columns.query == self.query,
            table.columns.dphi == self.dphi,
            table.columns.mechanism == self.mechanism,
            table.columns.wavelength == [w.value for w in self.wavelength],
            table.columns.fitted)
        
        with engine.connect() as con:
            result = pd.DataFrame(con.execute(query))

        # Should only have one match per outputfile
        # if len(result) != 1:
        #     from inspect import currentframe, getframeinfo
        #     frameinfo = getframeinfo(currentframe())
        #     print(frameinfo.filename, frameinfo.lineno)
        #     from IPython import embed; embed()
        #     import sys; sys.exit()
        if len(result) == 1:
            return result.loc[0, 'idnum'], ufit_id, result.loc[0, 'filename']
        elif len(result) == 0:
            return None
        else:
            assert False, 'Error'

    def determine_source_from_data(self, scdata, overwrite=False):
        """Determine the source using a previous LOSResult
        scdata = spacecraft data with at least one model result saved
        """
        unfit_model_result = scdata.model_result[self.unfitted_label]
        data = scdata.data

        if overwrite:
            self.inputs.delete_files()
        else:
            pass
        
        fitted_iteration_results = []
        print(f'LOSResultFitted: {len(unfit_model_result.outid)} unfitted files.')
        for ufit_id, ufit_outfile in zip(unfit_model_result.outid,
                                         unfit_model_result.outputfiles):
            # Check to see if there is already a result for this
            search_result = self.fitted_iteration_search(ufit_id)
            if search_result is None:
                # Need to compute for this unfit output file
                output = Output.restore(ufit_outfile)
                unfit_modelfile = unfit_model_result.modelfiles[ufit_outfile]
                with open(unfit_modelfile, 'rb') as file:
                    iteration_unfit = pickle.load(file)
                    
                # Remove packets that don't intersect the line of sight
                all_used_packets = set(np.concatenate(
                    iteration_unfit.used_packets.apply(list).values))
                # all_used_packets0 = set(np.concatenate(
                    # iteration_unfit.used_packets0.apply(list).values))

                output.X.frac *= output.X.index.isin(all_used_packets)
                output.X = output.X[output.X.frac > 0]

                packets = output.X.copy()
                packets0 = output.X0.copy()
                
                # radiance = fitted radiance
                # weighting = weighting factor for each X0
                # included = number of times trajectory included in lines of sight
                radiance = pd.Series(np.zeros(data.shape[0]), index=data.index)
                weighting = pd.Series(np.zeros(packets0.shape[0]),
                                      index=packets0.index)
                included = pd.Series(np.zeros(packets0.shape[0]),
                                     index=packets0.index)
                
                ratio = data.radiance / unfit_model_result.radiance
                ratio.fillna(0, inplace=True)

                # def find_where_used(ind):
                #     which = iteration_unfit.used_packets0.apply(lambda x:ind in x)
                #     return list(iteration_unfit.used_packets0.index[which])

                # packets0_used_by_spectrum = pd.Series(packets0.index).apply(
                #     find_where_used)
                # weight = packets0_used_by_spectrum.apply(lambda x:ratio[x].values)
                # weighting_ = weight.apply(lambda x: np.mean(x) if len(x) > 0 else 0)

                for spnum, spectrum in data.iterrows():
                    used = list(iteration_unfit.used_packets.loc[spnum])
                    cts = packets.loc[used, 'Index'].value_counts()
                    weighting.loc[cts.index] += cts.values * ratio[spnum]
                    included.loc[cts.index] += cts.values

                used = included > 0
                weighting[used] = weighting[used] / included[used]
                weighting /= weighting[used].mean()
                assert np.all(np.isfinite(weighting))

                multiplier = weighting.loc[output.X['Index']].values
                output.X.loc[:, 'frac'] = output.X.loc[:, 'frac'] * multiplier
                output.X0.loc[:, 'frac'] = output.X0.loc[:, 'frac'] * weighting
                output.totalsource = output.X0['frac'].sum() * output.nsteps
                packets = output.X.copy()
                packets['radvel_sun'] = (packets['vy'] +
                                         output.vrplanet.to(self.unit / u.s).value)
                
                self.packet_weighting(packets, output.aplanet)
                
                for spnum, spectrum in data.iterrows():
                    used = list(iteration_unfit.used_packets.loc[spnum])
                    
                    if len(used) > 0:
                        subset = packets.loc[used]
                        x_sc = spectrum[xcols].values.astype(float)
                        subset_rel_sc = subset[xcols].values - x_sc[np.newaxis, :]
                        subset_dist_sc = np.linalg.norm(subset_rel_sc, axis=1)
                        
                        Apix = np.pi * (subset_dist_sc * np.sin(self.dphi))**2 * (
                            self.unit.to(u.cm))**2
                        wtemp = subset['weight'] / Apix
                        radiance.loc[spnum] = wtemp.sum()
                    else:
                        pass

                # Save the fitted output
                output.inputs = self.inputs
                output.save()
                
                iteration = {'radiance': radiance.values,
                             'npackets': output.X0.frac.sum(),
                             'totalsource': output.totalsource,
                             'outputfile': output.filename,
                             'out_idnum': output.idnum,
                             'unfit_outputfile': ufit_outfile,
                             'unfit_outid': ufit_id,
                             'unfit_modelfile': unfit_modelfile}
                iteration_result = IterationResultFitted(iteration, self)
                iteration_result.save_iteration()
                fitted_iteration_results.append(iteration_result)

                del output, packets, packets0
            else:
                print(f'Using saved file {search_result[1]}')
                iteration_result = self.restore_iteration(search_result)
                assert len(iteration_result.radiance) == len(data)
                iteration_result.model_idnum = search_result[0]
                iteration_result.modelfile = search_result[2]
                fitted_iteration_results.append(iteration_result)

        self.modelfiles = {}
        self.outputfiles = []
        for iteration_result in fitted_iteration_results:
            self.radiance += iteration_result.radiance
            self.totalsource += iteration_result.totalsource
            self.modelfiles[iteration_result.outputfile] = (
                iteration_result.modelfile)
            self.outputfiles.append(iteration_result.outputfile)

        model_rate = self.totalsource/self.inputs.options.endtime.value
        self.atoms_per_packet = 1e23 / model_rate
        self.radiance *= self.atoms_per_packet/1e3*u.kR
        self.determine_source_rate(scdata)
        self.atoms_per_packet *= self.sourcerate.unit
        self.unfit_outputfiles = list(self.modelfiles.keys())

        print(self.totalsource, self.atoms_per_packet)
