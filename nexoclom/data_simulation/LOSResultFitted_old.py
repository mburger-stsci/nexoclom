import copy

import numpy as np
import pandas as pd
import pickle
import astropy.units as u
import sqlalchemy as sqla
from nexoclom import engine
from nexoclom.particle_tracking.Output import Output
from nexoclom.data_simulation.LOSResult import LOSResult
from nexoclom.data_simulation.compute_iteration import IterationResultFitted


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
        if len(result) == 1:
            return result.loc[0, 'idnum'], ufit_id, result.loc[0, 'filename']
        elif len(result) == 0:
            return None
        else:
            assert False, 'Error'

    def determine_source_from_data(self, scdata, overwrite=False,
                                   use_weight=False):
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
                if 'Index' not in output.X.columns:
                    output.X['Index'] = output.X.index
                else:
                    pass

                # if output.inputs.options.step_size != 0:
                #     # Save one random packet from each trajectory
                #     Xorig = output.X.copy()
                #     output.X['ind_'] = output.X.index
                #     times = output.X.time.unique()
                #     output.X.set_index(['Index', 'time'], inplace=True)
                #     index = set(output.X.index)
                #     steps = set(zip(np.arange(output.npackets),
                #                     output.randgen.choice(times, output.npackets)))
                #     steps = steps.intersection(index)
                #
                #     ind = pd.MultiIndex.from_tuples(steps,
                #                                     names=['Index', 'time'])
                #     output.X = output.X.loc[ind]
                #     output.X.reset_index(inplace=True)
                #     output.X.index = output.X.ind_
                #     output.X.drop(columns=['ind_'], inplace=True)
                # else:
                #     pass
                
                packets = output.X.copy()
                packets0 = output.X0.copy()

                unfit_modelfile = unfit_model_result.modelfiles[ufit_outfile]
                with open(unfit_modelfile, 'rb') as file:
                    iteration_unfit = pickle.load(file)
                    
                # radiance = fitted radiance
                radiance = pd.Series(np.zeros(data.shape[0]), index=data.index)
                
                # These are used if doing a weighted mean
                ratio_x_sigma = pd.Series(np.zeros(packets0.shape[0]),
                                          index=packets0.index)
                sigma = pd.Series(np.zeros(packets0.shape[0]),
                                  index=packets0.index)

                # These are used if not doing a weighted mean
                # weighting = weighting factor for each X0
                # included = number of times trajectory included in lines of sight
                weighting = pd.Series(np.zeros(packets0.shape[0]),
                                      index=packets0.index)
                included = pd.Series(np.zeros(packets0.shape[0]),
                                     index=packets0.index)
                used_by = {}
                ratio = data.radiance / unfit_model_result.radiance
                ratio.fillna(0, inplace=True)

                # Using weighted mean to determine adjustment
                mask = data[f'mask_{self.unfitted_label}']
                # import matplotlib.pyplot as plt
                # plt.scatter(packets.x, packets.y, marker='.', s=0.1)
                for spnum, spectrum in data[mask].iterrows():
                    to_use = [x for x in iteration_unfit.used_packets.loc[spnum]
                              if x in packets.index]
                    for x in to_use:
                        if x in used_by:
                            used_by[x].append(spnum)
                        else:
                            used_by[x] = [spnum]
                    # plt.scatter(packets.loc[to_use, 'x'], packets.loc[to_use, 'y'],
                    #             marker='.')
                    # plt.scatter(spectrum.x, spectrum.y, color='black')
                    # plt.plot([spectrum.x, spectrum.x + 10*spectrum.xbore],
                    #          [spectrum.y, spectrum.y + 10*spectrum.ybore],
                    #          color='black')
                    # print(len(to_use))
                    # plt.show(block=False)

                    cts = packets.loc[to_use, 'Index'].value_counts()

                    if use_weight:
                        ratio_x_sigma.loc[cts.index] += (cts.values * ratio[spnum] /
                                                         spectrum.sigma**2)
                        sigma.loc[cts.index] += cts.values / spectrum.sigma**2
                    else:
                        weighting.loc[cts.index] += cts.values * ratio[spnum]
                        included.loc[cts.index] += cts.values

                if use_weight:
                    used = sigma > 0
                    ratio_x_sigma[used] = ratio_x_sigma[used] / sigma[used]
                    weighting = ratio_x_sigma / ratio_x_sigma[used].mean()
                else:
                    used = included > 0
                    weighting[used] = weighting[used] / included[used]
                    weighting /= weighting[used].mean()
                assert np.all(np.isfinite(weighting))
                
                multiplier = weighting.loc[output.X['Index']].values
                output.X.loc[:, 'frac'] = output.X.loc[:, 'frac'] * multiplier
                output.X0.loc[:, 'frac'] = output.X0.loc[:, 'frac'] * weighting
                # if 'nsteps' in output.__dict__:
                #     output.totalsource = output.X0['frac'].sum() * output.nsteps
                # else:
                #     output.totalsource = output.X0['frac'].sum()
                output.totalsource = output.X0['frac'].sum()
                packets = output.X.copy()
                packets['radvel_sun'] = (packets['vy'] +
                                         output.vrplanet.to(self.unit / u.s).value)
                self.packet_weighting(packets, output.aplanet)

                for spnum, spectrum in data.iterrows():
                    to_use = [x for x in iteration_unfit.used_packets.loc[spnum]
                              if x in packets.index]

                    if len(to_use) > 0:
                        subset = packets.loc[to_use]
                        x_sc = spectrum[xcols].values.astype(float)
                        subset_rel_sc = subset[xcols].values - x_sc[np.newaxis, :]
                        subset_dist_sc = np.linalg.norm(subset_rel_sc, axis=1)

                        Apix = np.pi * (subset_dist_sc * np.sin(self.dphi))**2 * (
                            self.unit.to(u.cm))**2
                        wtemp = subset['weight'] / Apix
                        radiance.loc[spnum] = wtemp.sum()
                    else:
                        pass
                
                # if use_weight:
                #     output_weight = copy.deepcopy(output)
                #     radiance_weight = copy.deepcopy(radiance)
                #     del output
                # else:
                #     output_noweight = copy.deepcopy(output)
                #     radiance_noweight = copy.deepcopy(radiance)
                #     del output
                
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
                             'unfit_modelfile': unfit_modelfile,
                             'included': included > 0}
                iteration_result = IterationResultFitted(iteration, self)
                iteration_result.save_iteration()
                fitted_iteration_results.append(iteration_result)

                del output
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
