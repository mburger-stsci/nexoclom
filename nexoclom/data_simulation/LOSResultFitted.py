import numpy as np
import pandas as pd
import pickle
import astropy.units as u
import sqlalchemy as sqla
from scipy.spatial import distance_matrix
from sklearn.neighbors import BallTree
from nexoclom import engine
from nexoclom.particle_tracking.Output import Output
from nexoclom.data_simulation.LOSResult import LOSResult
from nexoclom.data_simulation.compute_iteration import IterationResultFitted


xcols = ['x', 'y', 'z']
borecols = ['xbore', 'ybore', 'zbore']


class LOSResultFitted(LOSResult):
    def __init__(self, scdata, label_for_fitted, params=None, dphi=1*u.deg,
                 **kwargs):
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
                                   use_weight=None):
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
        ct = 0
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

                if output.inputs.options.step_size != 0:
                    # Save one random packet from each trajectory
                    Xorig = output.X.copy()
                    output.X['ind_'] = output.X.index
                    times = output.X.time.unique()
                    output.X.set_index(['Index', 'time'], inplace=True)
                    index = set(output.X.index)
                    steps = set(zip(np.arange(output.npackets),
                                    output.randgen.choice(times, output.npackets)))
                    steps = steps.intersection(index)

                    ind = pd.MultiIndex.from_tuples(steps,
                                                    names=['Index', 'time'])
                    output.X = output.X.loc[ind]
                    output.X.reset_index(inplace=True)
                    output.X.index = output.X.ind_
                    output.X.drop(columns=['ind_'], inplace=True)
                else:
                    pass
                
                packets = output.X.copy()
                packets0 = output.X0.copy()

                unfit_modelfile = unfit_model_result.modelfiles[ufit_outfile]
                with open(unfit_modelfile, 'rb') as file:
                    iteration_unfit = pickle.load(file)
                    
                # radiance = fitted radiance
                radiance = pd.Series(np.zeros(data.shape[0]), index=data.index)
                
                # weighting = (distance from s/c)**(-2)
                # Final multiplier:
                #   f_i = sum_j(w_ij * m_ij)/sum_ij(w_ij) for packet0_i and
                #       specutrm_j
                
                ratio_x_sigma = pd.Series(np.zeros(packets0.shape[0]),
                                          index=packets0.index)
                sigma = pd.Series(np.zeros(packets0.shape[0]),
                                  index=packets0.index)

                ratio = data.radiance / unfit_model_result.radiance
                ratio.fillna(0, inplace=True)

                # Using weighted mean to determine adjustment
                mask = data[f'mask_{self.unfitted_label}']
                
                used_by, used_by0 = {}, {}
                # print(f'Weighting = {use_weight}')
                for spnum, spectrum in data[mask].iterrows():
                    to_use = [x for x in iteration_unfit.used_packets.loc[spnum]
                              if x in packets.index]
                    for x in to_use:
                        if x in used_by:
                            used_by[x].append(spnum)
                        else:
                            used_by[x] = [spnum]

                        x0 = packets.loc[x, 'Index']
                        if x0 in used_by0:
                            used_by0[packets.loc[x, 'Index']].append(spnum)
                        else:
                            used_by0[x0] = [spnum]

                    if use_weight == 'dist2':
                        sc_dist = np.sqrt((packets.loc[to_use, 'x'] - spectrum.x)**2 +
                                          (packets.loc[to_use, 'y'] - spectrum.y)**2 +
                                          (packets.loc[to_use, 'z'] - spectrum.z)**2)
                        weight = 1/sc_dist**2
                    elif use_weight == 'dist':
                        sc_dist = np.sqrt((packets.loc[to_use, 'x'] - spectrum.x)**2 +
                                          (packets.loc[to_use, 'y'] - spectrum.y)**2 +
                                          (packets.loc[to_use, 'z'] - spectrum.z)**2)
                        weight = 1/sc_dist
                    elif use_weight == 'sigma':
                        weight = pd.Series(np.ones(len(to_use))/spectrum.sigma*2,
                                           index=packets.loc[to_use].index)
                    else:
                        weight = pd.Series(np.ones(len(to_use)),
                                           index=packets.loc[to_use].index)
                    
                    for tu in to_use:
                        ind0 = packets.loc[tu, 'Index']
                        ratio_x_sigma.loc[ind0] += ratio[spnum] * weight[tu]
                        sigma.loc[ind0] += weight[tu]
                    
                used = sigma > 0
                ratio_x_sigma[used] = ratio_x_sigma[used] / sigma[used]
                weighting = ratio_x_sigma / ratio_x_sigma[used].mean()
                
                multiplier = weighting.loc[output.X['Index']].values
                output.X.loc[:, 'frac'] = output.X.loc[:, 'frac'] * multiplier
                output.X0.loc[:, 'frac'] = output.X0.loc[:, 'frac'] * weighting

                output.totalsource = output.X0['frac'].sum() * output.nsteps
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
                             'included': True}
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
            print(f'Completed {ct} of {len(unfit_model_result.outid)}')

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
        self.determine_source_rate(scdata, use_weight=False)
        self.atoms_per_packet *= self.sourcerate.unit * u.s
        self.unfit_outputfiles = list(self.modelfiles.keys())

        print(self.totalsource, self.atoms_per_packet)
