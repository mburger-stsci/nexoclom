import os.path
import numpy as np
import pandas as pd
import pickle
import random
import copy
import astropy.units as u
from sklearn.neighbors import KDTree, BallTree
from scipy.interpolate import RectBivariateSpline
from mathMB import fit_model

from .ModelResults import ModelResult
from .database_connect import database_connect
from .Input import Input
from .Output import Output


xcols = ['x', 'y', 'z']
borecols = ['xbore', 'ybore', 'zbore']
NLONBINS, NLATBINS, NVELBINS = 72, 36, 100


class DataPackets:
    def __init__(self, *args, **kwargs):
        self.data = pd.DataFrame(*args, columns=['outputfile', 'specind',
                                                 'oind', 'weight', 'frac0',
                                                 'Index0'], *kwargs)
        self.totalsource = 0.
        self.atoms_per_packet = 0.
        self.outputfile = None
        self.out_idnum = None
        

class InputError(Exception):
    """Raised when a required parameter is not included."""
    def __init__(self, expression, message):
        self.expression = expression
        self.message = message


class LOSResult(ModelResult):
    """Class to contain the LOS result from multiple outputfiles."""
    def __init__(self, scdata, quantity='radiance', dphi=1*u.deg):
        """Determine column or emission along lines of sight.
        This assumes the model has already been run.
        
        Parameters
        ==========
        scdata
            Spacecraft data object (currently designed for MESSENGERdata object
            but can be faked for other types of data)
            
        quantity
            Quantity to calculate: 'column', 'radiance', 'density'
            
        dphi
            Angular size of the view cone. Default = r deg.
        """
        format_ = {'quantity':quantity}
        super().__init__(format_, species=scdata.species)
        
        # Basic information
        self.scdata = scdata
        self.scdata.set_frame('Model')
        self.type = 'LineOfSight'
        self.unit = None
        self.dphi = dphi.to(u.rad).value
        self.oedge = None
        
        self.fitted = None
        self.outputfiles = []
        nspec = len(self.scdata)
        self.radiance = pd.Series(np.zeros(nspec), index=scdata.data.index)
        self.npackets = pd.Series(np.zeros(nspec), index=scdata.data.index,
                                  dtype=int)
        self.totalsource = 0.
        self.sourcemap = None
        self.modelfiles = {}
        self.data_packets = DataPackets()

    def delete_models(self):
        """Deletes any LOSResult models associated with this data and input
        This may never actually do anything. Overwrite=True will also
        erase the outputfiles (which erases any models that depend on them).
        Unless I put separate outputfile and modelfile delete switches,
        This shouldn't do anything"""
        
        search_results = self.search()
        if len(search_results) != 0:
            print('Warning: LOSResult.delete_models found something to delete')
            for _, search_result in search_results.items():
                if search_result is not None:
                    idnum, modelfile = search_result
                    with database_connect() as con:
                        cur = con.cursor()
                        cur.execute(f'''DELETE from uvvsmodels
                                       WHERE idnum = %s''', (idnum, ))
                    if os.path.exists(modelfile):
                        os.remove(modelfile)
        else:
            pass

    def save(self, data_packets):
        # Insert the model into the database
        # Save is on an outputfile basis
        if self.quantity == 'radiance':
            mech = ', '.join(sorted([m for m in self.mechanism]))
            wave_ = sorted([w.value for w in self.wavelength])
            wave = ', '.join([str(w) for w in wave_])
        else:
            mech = None
            wave = None

        # Save query with all white space removed and lowercase
        tempname = f'temp_{str(random.randint(0, 1000000))}'

        with database_connect() as con:
            cur = con.cursor()
            cur.execute(f'''INSERT into uvvsmodels (out_idnum, quantity,
                            query, dphi, mechanism, wavelength,
                            fitted, filename)
                            values (%s, %s, %s, %s, %s, %s, %s, %s)''',
                        (data_packets.out_idnum, self.quantity,
                         self.scdata.query, self.dphi, mech, wave, self.fitted,
                         tempname))

            # Determine the savefile name
            idnum_ = pd.read_sql(f'''SELECT idnum
                                     FROM uvvsmodels
                                     WHERE filename='{tempname}';''', con)
            assert len(idnum_) == 1
            idnum = int(idnum_.idnum[0])

            savefile = os.path.join(os.path.dirname(data_packets.outputfile),
                                    f'model.{idnum}.pkl')

            cur.execute(f'''UPDATE uvvsmodels
                            SET filename=%s
                            WHERE idnum=%s''', (savefile, idnum))
            
        with open(savefile, 'wb') as f:
            pickle.dump(data_packets, f)
            
        return savefile

    def search(self):
        """
        :return: dictionary containing search results:
                 {outputfilename: (modelfile_id, modelfile_name)}
        """
        search_results = {}
        for outputfile in self.outputfiles:
            with database_connect() as con:
                # Determine the id of the outputfile
                idnum_ = pd.read_sql(
                    f'''SELECT idnum
                        FROM outputfile
                        WHERE filename='{outputfile}' ''', con)
                oid = idnum_.idnum[0]
            
                if self.quantity == 'radiance':
                    mech = ("mechanism = '" +
                            ", ".join(sorted([m for m in self.mechanism])) +
                            "'")
                    wave_ = sorted([w.value for w in self.wavelength])
                    wave = ("wavelength = '" +
                            ", ".join([str(w) for w in wave_]) +
                            "'")
                else:
                    mech = 'mechanism is NULL'
                    wave = 'wavelength is NULL'
            
                result = pd.read_sql(
                    f'''SELECT idnum, filename FROM uvvsmodels
                        WHERE out_idnum={oid} and
                              quantity = '{self.quantity}' and
                              query = '{self.scdata.query}' and
                              dphi = {self.dphi} and
                              {mech} and
                              {wave} and
                              fitted = {self.fitted}''', con)
            
                # Should only have one match per outputfile
                assert len(result) <= 1
                
                if len(result) == 0:
                    search_results[outputfile] = None
                else:
                    search_results[outputfile] = (result.iloc[0, 0],
                                                  result.iloc[0, 1])
                
        return search_results
    
    @staticmethod
    def restore(modelfile):
        # Restore is on an outputfile basis
        with open(modelfile, 'rb') as f:
            data_packets = pickle.load(f)

        return data_packets
    
    def _data_setup(self):
        # distance of s/c from planet
        data = self.scdata.data
        dist_from_plan = np.sqrt(data.x**2 + data.y**2 + data.z**2)
    
        # Angle between look direction and planet.
        ang = np.arccos((-data.x*data.xbore - data.y*data.ybore -
                         data.z*data.zbore)/dist_from_plan)
        
        # Check to see if look direction intersects the planet anywhere
        asize_plan = np.arcsin(1./dist_from_plan)

        # Don't worry about lines of sight that don't hit the planet
        dist_from_plan.loc[ang > asize_plan] = 1e30
        
        return dist_from_plan
    
    def _spectrum_process(self, spectrum, packets, tree, dist, i, ofile):
        x_sc = spectrum[xcols].values.astype(float)
        bore = spectrum[borecols].values.astype(float)
    
        dd = 30  # Furthest distance we need to look
        x_far = x_sc + bore * dd
        while np.linalg.norm(x_far) > self.oedge:
            dd -= 0.1
            x_far = x_sc + bore * dd
    
        t = [0.05]
        while t[-1] < dd:
            t.append(t[-1] + t[-1] * np.sin(self.dphi))
        t = np.array(t)
        Xbore = x_sc[np.newaxis, :] + bore[np.newaxis, :] * t[:, np.newaxis]
    
        wid = t * np.sin(self.dphi)
        ind = np.concatenate(tree.query_radius(Xbore, wid))
        ilocs = np.unique(ind).astype(int)
        indicies = packets.iloc[ilocs].index
        subset = packets.loc[indicies]
    
        xpr = subset[xcols] - x_sc[np.newaxis, :]
        rpr = np.sqrt(xpr['x'] * xpr['x'] +
                      xpr['y'] * xpr['y'] +
                      xpr['z'] * xpr['z'])
    
        losrad = np.sum(xpr * bore[np.newaxis, :], axis=1)
        inview = rpr < dist
        
        if np.any(inview):
            used_packets = subset[inview].index.to_list()
            used_packets0 = subset.loc[inview, 'Index']

            Apix = np.pi * (rpr[inview] * np.sin(self.dphi))**2 * (
                self.unit.to(u.cm))**2
            wtemp = subset.loc[inview, 'weight'] / Apix
            if self.quantity == 'radiance':
                # Determine if any packets are in shadow
                # Projection of packet onto LOS
                # Point along LOS the packet represents
                losrad_ = losrad[inview].values
                hit = (x_sc[np.newaxis, :] +
                       bore[np.newaxis, :] * losrad_[:, np.newaxis])
                rhohit = np.linalg.norm(hit[:, [0, 2]], axis=1)
                out_of_shadow = (rhohit > 1) | (hit[:, 1] < 0)
                wtemp *= out_of_shadow

                rad = wtemp.sum()
                # Save the weight information
                rat = spectrum.radiance/rad if rad > 0 else 0.
                processed = pd.DataFrame(
                    {'outputfile': [ofile for _ in used_packets],
                     'specind': [i for _ in used_packets],
                     'oind': used_packets,
                     'weight': wtemp.values,
                     'frac0': packets.loc[used_packets, 'frac0'].values,
                     'Index0': used_packets0.values,
                     'ratio': [rat for _ in used_packets]})
            else:
                assert False, 'Other quantities not set up.'
        else:
            processed = None
        
        return processed

    @staticmethod
    def _tree(values, treetype='KDTree'):
        if treetype == 'KDTree':
            return KDTree(values)
        elif treetype == 'BallTree':
            return BallTree(values)

    def determine_source_from_data(self, scdata, modkey, masking=None):
        self.modelfiles = scdata.model_info[modkey]['modelfiles']
        self.outputfiles = self.modelfiles.keys()
        surface_source = pd.DataFrame(columns=['outputfile', 'Index', 'longitude',
                                               'latitude', 'velocity'])
        data = scdata.data
        # Determine which points should be used for the fit
        _, _, mask = fit_model(data.radiance, None, data.sigma,
                               masking=masking, mask_only=True,
                               altitude=data.alttan)

        self.data_packets = DataPackets()
        nsteps = 0
        assert len(self.modelfiles) > 0
        for outputfile, modefile in self.modelfiles.items():
            # Restore the unfit output file
            output = Output.restore(outputfile)
            self.totalsource += output.totalsource
            endtime = output.inputs.options.endtime
            nsteps = output.nsteps
            data_packets_it = self.restore(modefile)
            
            self.data_packets.data = self.data_packets.data.append(
                data_packets_it.data)
            self.data_packets.totalsource += data_packets_it.totalsource
            
            vel_ = np.sqrt(output.X0.vx**2 + output.X0.vy**2 +
                           output.X0.vz**2) * output.inputs.geometry.planet.radius
            surface_source = surface_source.append(pd.DataFrame(
                {'outputfile': [outputfile for _ in range(len(output.X0))],
                 'Index': output.X0.index.values,
                 'longitude': output.X0.longitude.values,
                 'latitude': output.X0.latitude.values,
                 'velocity': vel_}), ignore_index=True)
            
            del output
            del data_packets_it
            
        # Determine the proper weightings
        spec_group = self.data_packets.data.groupby('specind')
        rad_ = spec_group['weight'].sum()
        rad = pd.Series(np.zeros(data.shape[0]), index=data.index)
        rad.loc[rad_.index] = rad_

        spec_ind = self.data_packets.data.set_index('specind')
        spec_ind = spec_ind.loc[mask]
        scale_factor_ = spec_ind.groupby(['outputfile', 'Index0'])['ratio'].mean()
        scale_factor_ /= scale_factor_[scale_factor_ > 0].mean()
        
        self.data_packets.data['scale_factor'] = np.zeros(len(self.data_packets.data))
        ind0_ind = self.data_packets.data.set_index(['outputfile', 'Index0'])
        ind0_ind.loc[scale_factor_.index, 'scale_factor'] = scale_factor_
        self.data_packets.data = ind0_ind.reset_index()
        self.data_packets.data['fit_weight'] = (self.data_packets.data.weight *
                                                self.data_packets.data.scale_factor)
        
        initial = set(zip(self.data_packets.data.outputfile.to_list(),
                          self.data_packets.data.Index0.to_list(),
                          self.data_packets.data.frac0.to_list(),
                          self.data_packets.data.scale_factor.to_list()))
        frac0 = sum(x[2]*x[3] for x in initial)
        self.totalsource = frac0 * nsteps
        mod_rate = self.totalsource / endtime.value
        atoms_per_packet = 1e23 / mod_rate
        
        spec_group = self.data_packets.data.groupby('specind')
        fitrad_ = spec_group['fit_weight'].sum() * atoms_per_packet
        
        self.radiance.loc[fitrad_.index] = fitrad_
        self.radiance *= u.R

        surface_ind = surface_source.set_index(['outputfile', 'Index'])
        surface_ind.loc[scale_factor_.index, 'weight'] = scale_factor_
        surface_ind.dropna(inplace=True)
        surface_source = surface_ind.reset_index()

        # Create a new sourcemap
        source, xx, yy = np.histogram2d(surface_source['longitude'],
                                        surface_source['latitude'],
                                        weights=surface_source['weight'],
                                        range=[[0, 2*np.pi], [-np.pi/2, np.pi/2]],
                                        bins=(NLONBINS, NLATBINS))
        source = source/np.cos(yy+(yy[1]-yy[0])/2.)[np.newaxis,:-1]
        source[:,[0,-1]] = 0
        
        v_source, v = np.histogram(surface_source['velocity'], bins=NVELBINS,
                                   range=[0, surface_source['velocity'].max()],
                                   weights=surface_source['weight'])
        v_source /= np.max(v_source)

        # packets available
        packets, _, _ = np.histogram2d(surface_source['longitude'],
                                       surface_source['latitude'],
                                       range=[[0, 2*np.pi], [-np.pi/2, np.pi/2]],
                                       bins=(NLONBINS, NLATBINS))
        packets = packets/np.cos(yy+(yy[1]-yy[0])/2.)[np.newaxis,:-1]
        packets[:,[0,-1]] = 0
        v_packets, _ = np.histogram(surface_source['velocity'], bins=NVELBINS,
                                    range=[0, surface_source['velocity'].max()])
        v_packets = v_packets / np.max(v_packets)
        
        sourcemap = {'longitude': xx*u.rad,
                     'latitude': yy*u.rad,
                     'abundance': source,
                     'p_available': packets,
                     'velocity': v*u.km/u.s,
                     'vdist': v_source,
                     'v_available': v_packets,
                     'coordinate_system': 'solar-fixed'}
        self.sourcemap = sourcemap

    def simulate_data_from_inputs(self, inputs_, npackets, overwrite=False,
                                  packs_per_it=None):
        """Given a set of inputs, determine what the spacecraft should see.
        
        Parameters
        ==========
        inputs_
            A nexoclom Input object or the name of an inputs file
        """
        if isinstance(inputs_, str):
            self.inputs = Input(inputs_)
        elif isinstance(inputs_, Input):
            self.inputs = copy.deepcopy(inputs_)
        else:
            raise InputError('nexoclom.LOSResult', 'Problem with the inputs.')
    
        # TAA needs to match the data
        self.inputs.geometry.taa = self.scdata.taa
        self.unit = u.def_unit('R_' + self.inputs.geometry.planet.object,
                               self.inputs.geometry.planet.radius)
    
        # If using a planet-fixed source map, need to set subsolarlon
        if ((self.inputs.spatialdist.type == 'surface map') and
            (self.inputs.spatialdist.coordinate_system == 'planet-fixed')):
            self.inputs.spatialdist.subsolarlon = self.scdata.subslong.median() * u.rad
        else:
            pass
    
        # Run the model
        self.fitted = False
        self.inputs.run(npackets, packs_per_it=packs_per_it, overwrite=overwrite)
        self.search_for_outputs()
    
        data = self.scdata.data
        search_results = self.search()
    
        dist_from_plan = (self._data_setup()
                          if None in search_results.values()
                          else None)
        
        for outputfile, search_result in search_results.items():
            if search_result is None:
                # simulate the data
                output = Output.restore(outputfile)
                data_packets_it = DataPackets()
                data_packets_it.outputfile = outputfile
                data_packets_it.out_idnum = output.idnum
                data_packets_it.totalsource = output.totalsource

                # need model rate for this output
                mod_rate = output.totalsource / self.inputs.options.endtime.value
        
                packets = copy.deepcopy(output.X)
                packets['frac0'] = output.X0.loc[packets['Index'], 'frac'].values
                packets['radvel_sun'] = (packets['vy'] +
                                         output.vrplanet.to(self.unit / u.s).value)
                self.oedge = output.inputs.options.outeredge * 2
        
                # Will base shadow on line of sight, not the packets
                out_of_shadow = np.ones(len(packets))
                self.packet_weighting(packets, out_of_shadow, output.aplanet)
        
                # This sets limits on regions where packets might be
                tree = self._tree(packets[xcols].values)
        
                print(f'{data.shape[0]} spectra taken.')
                for i, spectrum in data.iterrows():
                    processed = self._spectrum_process(spectrum, packets, tree,
                                                       dist_from_plan[i],
                                                       i, outputfile)
                    if processed is not None:
                        data_packets_it.data = data_packets_it.data.append(
                            processed, ignore_index=True)
                    else:
                        pass
                    
                    ind = data.index.get_loc(i)
                    if (ind % (len(data)//10)) == 0:
                        print(f'Completed {ind+1} spectra')
                    else:
                        pass

                self.data_packets.data = self.data_packets.data.append(
                    data_packets_it.data)
                self.data_packets.totalsource += data_packets_it.totalsource
                
                modelfile = self.save(data_packets_it)
                self.modelfiles[outputfile] = modelfile
                
                del output
                del data
                del data_packets_it
            else:
                data_packets_it = self.restore(search_result[1])
                self.data_packets.data = self.data_packets.data.append(
                    data_packets_it.data)
                self.data_packets.totalsource += data_packets_it.totalsource
                self.modelfiles[outputfile] = search_result[1]
                
                del data_packets_it

        mod_rate = self.data_packets.totalsource / self.inputs.options.endtime.value
        atoms_per_packet = 1e23 / mod_rate
        
        spec_group = self.data_packets.data.groupby('specind')
        rad_ = spec_group['weight'].sum() * atoms_per_packet
        self.radiance.loc[rad_.index] = rad_
        self.radiance *= u.R


def recompute_radiance(scdata, modnum, new_source, mask=None):
    data = scdata.data
    radiance = pd.Series(index=data.index)

    modkey = f'model{modnum:02d}'
    maskkey = f'mask{modnum:02d}'
    model_info = scdata.model_info[modkey]
    mask = scdata.data[maskkey]

    for outputfile in data.model_info['outputfiles']:
        output = Output.restore(outputfile)
        
        usedpackets = model_info['usedpackets'][outputfile]
        usedpackets0 = model_info['usedpackets0'][outputfile]
        
        # Reset fractional values
        output.X.frac = output.X.frac/output.X0.loc[output.X.Index, 'frac'].values
        output.X0.loc[:,'frac'] = 1.
        
        # Determine new weighting from given map
        new_weight_fn = RectBivariateSpline(new_source['longitude'][:-1].value,
                                            new_source['latitude'][:-1].value,
                                            new_source['abundance'])
        new_weight = pd.Series(new_weight_fn(output.X0.longitude, output.X0.latitude,
                               grid=False), index=output.X0.index)
        
        multiplier = pd.Series(new_weight.loc[output.X['Index']].values,
                               index=output.X.index)
        # Need to determine proper scale factor

        # for i, row in data.iterrows():
            # radiance.loc[i] += np.sum(row.weight * multiplier.loc[].values)
    
    strength, goodness_of_fit, mask = fit_model(data.radiance.values,
                                                radiance,
                                                data.sigma.values,
                                                fit_method='chisq', mask=mask)
    radiance *= strength
    
    return radiance, strength, goodness_of_fit, mask

# import matplotlib.pyplot as plt
# hist, xx, yy = np.histogram2d(output.X0.longitude, output.X0.latitude,
#                               weights=new_weight,
#                               range=[[0, 2*np.pi], [-np.pi/2, np.pi/2]],
#                               bins=(NLONBINS, NLATBINS))
# hist = hist / np.cos(yy + (yy[1] - yy[0]) / 2.)[np.newaxis, :-1]
# hist /= hist.mean()
#
# fig, ax = plt.subplots(1, 2)
# ax[0].imshow(new_source['abundance'], vmin=0, vmax=3)
# ax[1].imshow(hist/hist.mean(), vmin=0, vmax=3)
# plt.show()
#
# plt.plot(np.sum(new_source['abundance'], axis=0))
# plt.plot(np.sum(hist, axis=0))
# plt.show()
# plt.plot(np.sum(new_source['abundance'], axis=1))
# plt.plot(np.sum(hist, axis=1))
# plt.show()
#
# from IPython import embed; embed()
# import sys; sys.exit()
