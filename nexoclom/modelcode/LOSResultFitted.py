import numpy as np
import pandas as pd
import pickle
import astropy.units as u
import sqlalchemy as sqla
from sklearn.neighbors import KDTree
from nexoclom import Output
from nexoclom.modelcode.LOSResult import LOSResult, IterationResult


xcols = ['x', 'y', 'z']
borecols = ['xbore', 'ybore', 'zbore']


class IterationResultFitted(IterationResult):
    def __init__(self, iteration):
        super().__init__(iteration)
        
        self.unfit_outputfile = iteration['unfit_outputfile']
        self.unfit_outid = iteration['unfit_outid']
        self.fitted = True

class LOSResultFitted(LOSResult):
    def __init__(self, scdata, inputs, params=None, dphi=1*u.deg, **kwargs):
        super().__init__(scdata, inputs, params=params, dphi=dphi, **kwargs)

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
        engine = self.inputs.config.create_engine()
        metadata_obj = sqla.MetaData()
        table = sqla.Table("uvvsmodels", metadata_obj, autoload_with=engine)

        query = sqla.select(table).where(
            table.columns.unfit_idnum == ufit_id,
            table.columns.quantity == self.quantity,
            table.columns.query == self.query,
            table.columns.dphi == self.dphi,
            table.columns.mechanism == self.mechanism,
            table.columns.wavelength == [w.value for w in self.wavelength],
            table.columns.fitted is True)
        
        with engine.connect() as con:
            result = pd.DataFrame(con.execute(query))

        # Should only have one match per outputfile
        if len(result) == 1:
            return result.loc[0, 'idnum'], ufit_id, result.loc[0, 'filename']
        elif len(result) == 0:
            return None
        else:
            assert False, 'Error'

    def determine_source_from_data(self, scdata, label):
        """Determine the source using a previous LOSResult
        scdata = spacecraft data with at least one model result saved
        label = Label for the model result that should be used as the starting point
        """
        unfit_model_result = scdata.model_result[label]
        data = scdata.data
        
        fitted_iteration_results = []
        for ufit_id, ufit_outfile in zip(unfit_model_result.outid,
                                         unfit_model_result.outputfiles):
            # Check to see if there is already a result for this
            search_result = self.fitted_iteration_search(ufit_id)
            
            if search_result is None:
                # Need to compute for this unfit output file
                output = Output.restore(ufit_outfile)
                with open(unfit_model_result.modelfiles[ufit_outfile], 'rb') as file:
                    iteration_unfit = pickle.load(file)
                assert output.compress is False, ('nexoclom.LOSResult: '
                                                  'Fitted results must start from uncompressed outputs')
                
                packets = output.X.copy()
                packets0 = output.X0.copy()
                
                radiance = pd.Series(np.zeros(data.shape[0]), index=data.index)
                weighting = pd.Series(np.zeros(packets0.shape[0]),
                                      index=packets0.index)
                included = pd.Series(np.zeros(packets0.shape[0]),
                                     index=packets0.index)
                
                ratio = data.radiance / unfit_model_result.radiance
                ratio.fillna(0, inplace=True)
                
                for spnum, spectrum in data.iterrows():
                    used = list(iteration_unfit.used_packets.loc[spnum])
                    # used0 = list(iteration_unfit.used_packets0.loc[spnum])
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
                             'unfit_outid': ufit_id}
                iteration_result = IterationResultFitted(iteration)
                
                modelfile = self.save(iteration_result, ufit_id=ufit_id)
                iteration_result.modelfile = modelfile
                fitted_iteration_results.append(iteration_result)
                
                del output, packets, packets0
            else:
                print(f'Using saved file {search_result[1]}')
                iteration_result = self.restore(search_result)
                assert len(iteration_result.radiance) == len(data)
                iteration_result.model_idnum = search_result[0]
                iteration_result.modelfile = search_result[2]
                fitted_iteration_results.append(iteration_result)
            self.modelfiles = {}
            
        for iteration_result in fitted_iteration_results:
            self.radiance += iteration_result.radiance
            self.totalsource += iteration_result.totalsource
            self.modelfiles[iteration_result.outputfile] = iteration_result.modelfile
        
        model_rate = self.totalsource/self.inputs.options.endtime.value
        self.atoms_per_packet = 1e23 / model_rate
        self.radiance *= self.atoms_per_packet/1e3*u.kR
        self.determine_source_rate(scdata)
        self.atoms_per_packet *= self.sourcerate.unit
        self.outputfiles = list(self.modelfiles.keys())

        print(self.totalsource, self.atoms_per_packet)
