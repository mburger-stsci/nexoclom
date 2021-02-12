import os.path
import numpy as np
import pandas as pd
import pickle
import random
import astropy.units as u
from datetime import datetime
from sklearn.neighbors import KDTree

from mathMB import fit_model

from .ModelResults import ModelResult
from .database_connect import database_connect
from .Input import Input
from .Output import Output


def determine_model_type(data, species):
    # Determine whether model is one complete orbit or a more complicated query
    from MESSENGERuvvs import MESSENGERdata
    if len(data.data.orbit.unique() == 1):
        orb = data.data.orbit.unique()[0]
        mdata = MESSENGERdata(species, f'orbit = {orb}')
        if len(mdata) == len(data):
            tempname = f'temp_{orb}_{str(random.randint(0, 1000000))}'
            field = 'orbit'
            quant = int(orb)
        else:
            tempname = f'temp_query_{str(random.randint(0, 1000000))}'
            field = 'query'
            quant = data.query
    else:
        tempname = f'temp_query_{str(random.randint(0, 1000000))}'
        field = 'query'
        quant = data.query
        
    return tempname, field, quant


class LOSResult(ModelResult):
    '''Class to contain the LOS result from multiple outputfiles.'''
    def __init__(self, start_from, mdata, quantity, dphi=3*u.deg,
                 filenames=None, overwrite=False, fit_to_data=False, **kwargs):
        """Determine column or emission along lines of sight.
        This assumes the model has already been run.
        
        Parameters
        ==========
        start_from
            Either an Input or Output object
        
        data
            A Pandas DataFrame object with information on the lines of sight.
            
        quantity
            Quantity to calculate: 'column', 'radiance', 'density'
            
        dphi
            Angular size of the view cone. Default = 3 deg.
            
        filenames
            A filename or list of filenames to use. Default = None is to
            find all files created for the inputs. Not used when start_from
            is an Output object.
            
        overwrite
            If True, deletes any images that have already been computed. Not
            used when start_from is an Output object.
            Default = False
        """
        format_ = {'quantity':quantity}
        # super().__init__ sets totalsource, output filenames
        if isinstance(start_from, Input):
            # Given inputs: Will need to search for and/or compute outputs
            inputs = start_from
            output = None
            super().__init__(inputs, format_, filenames=filenames)
        elif isinstance(start_from, Output):
            # Starting from an precomputed output
            output = start_from
            inputs = output.inputs
            super().__init__(inputs, format_, output=output)
        else:
            raise TypeError
        
        # Make sure data is in the right reference frame
        mdata.set_frame('Model')

        tstart = datetime.now()
        
        # Basic information
        self.type = 'LineOfSight'
        self.species = inputs.options.species
        self.origin = inputs.geometry.planet
        self.unit = u.def_unit('R_' + self.origin.object,
                               self.origin.radius)
        self.dphi = dphi.to(u.rad).value
        
        # Just compute lines of sight or also try to fit
        self.fitted = fit_to_data

        nspec = len(mdata)
        self.radiance = np.zeros(nspec)
        self.packets = pd.DataFrame()
        self.ninview = np.zeros(nspec, dtype=int)
        new_total_source = 0.
        saved_files = {}
        
        # self.filenames are the output files appropriate for this dataset
        for j,outfile in enumerate(self.filenames):
            if outfile is None:
                # Output was provided - outfile not known
                if self.fitted:
                    # Fit the outputs to the data
                    masking = kwargs['masking'] if 'masking' in kwargs else None
                    model_result = self.fit_model_to_data(mdata.data, output,
                                                          masking=masking)
                else:
                    # Just compute lines of sight
                    model_result = self.create_model(mdata.data, output)
                print(f'Completed model {j+1} of {len(self.filenames)}')
                del output
            else:
                # Search to see if it is already done
                model_result = self.restore(mdata, outfile)

                if (model_result['radiance'] is None) or overwrite:
                    if (model_result['radiance'] is not None) and overwrite:
                        self.delete_model(mdata, model_result['idnum'])
                    else:
                        pass

                    output = Output.restore(outfile)
                    if self.fitted:
                        masking = kwargs['masking'] if 'masking' in kwargs else None
                        model_result = self.fit_model_to_data(mdata.data,
                                                       output, masking=masking)
                        saved_file = self.save(mdata, outfile, model_result)
                    else:
                        model_result = self.create_model(mdata.data, output)
                        saved_file = self.save(mdata, outfile, model_result)
                    saved_files[outfile] = saved_file
                    print(f'Completed model {j+1} of {len(self.filenames)}')
                    del output
                else:
                    saved_files[outfile] = model_result['savefile']
                    print(f'Model {j+1} of {len(self.filenames)} '
                          'previously completed.')

            self.radiance += model_result['radiance'].values
            new_total_source += model_result['model_total_source']

        self.savefiles = saved_files
        self.radiance *= u.R
        tend = datetime.now()
        print(f'Total time = {tend-tstart}')

    def delete_model(self, data, idnum):
        _, field, _ = determine_model_type(data, self.species)
    
        with database_connect() as con:
            cur = con.cursor()
            cur.execute(f'''SELECT idnum, filename FROM uvvsmodels_{field}
                           WHERE out_idnum = %s''', (int(idnum), ))
            assert cur.rowcount in (0, 1)
            for mid, mfile in cur.fetchall():
                cur.execute(f'''DELETE from uvvsmodels_{field}
                               WHERE idnum = %s''', (mid, ))
                if os.path.exists(mfile):
                    os.remove(mfile)

    def save(self, data, fname, model_result):
        with database_connect() as con:
            cur = con.cursor()

            # Determine the id of the outputfile
            idnum_ = pd.read_sql(f'''SELECT idnum
                                    FROM outputfile
                                    WHERE filename='{fname}' ''', con)
            idnum = int(idnum_.idnum[0])

            # Insert the model into the database
            if self.quantity == 'radiance':
                mech = ', '.join(sorted([m for m in self.mechanism]))
                wave_ = sorted([w.value for w in self.wavelength])
                wave = ', '.join([str(w) for w in wave_])
            else:
                mech = None
                wave = None

            # Determine if model represents a single complete orbit
            tempname, field, quant = determine_model_type(data, self.species)

            cur.execute(f'''INSERT into uvvsmodels_{field} (out_idnum, quantity,
                            {field}, dphi, mechanism, wavelength,
                            fitted, filename)
                            values (%s, %s, %s, %s, %s, %s, %s, %s)''',
                        (idnum, self.quantity, quant, self.dphi,
                         mech, wave, self.fitted, tempname))

            # Determine the savefile name
            idnum_ = pd.read_sql(f'''SELECT idnum
                                     FROM uvvsmodels_{field}
                                     WHERE filename='{tempname}';''', con)
            assert len(idnum_) == 1
            idnum = int(idnum_.idnum[0])

            savefile = os.path.join(os.path.dirname(fname),
                                    f'model.{idnum}.pkl')
            with open(savefile, 'wb') as f:
                pickle.dump(model_result, f)
                
            cur.execute(f'''UPDATE uvvsmodels_{field}
                            SET filename=%s
                            WHERE idnum=%s''', (savefile, idnum))
    
        return savefile

    def restore(self, data, fname):
        tempname, field, quant = determine_model_type(data, self.species)
        
        with database_connect() as con:
            # Determine the id of the outputfile
            idnum_ = pd.read_sql(f'''SELECT idnum
                                    FROM outputfile
                                    WHERE filename='{fname}' ''', con)
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

            if field == 'query':
                quant = f"'{quant}'"
            else:
                pass
            result = pd.read_sql(
                f'''SELECT idnum, filename FROM uvvsmodels_{field}
                    WHERE out_idnum={oid} and
                          quantity = '{self.quantity}' and
                          {field} = {quant} and
                          dphi = {self.dphi} and
                          {mech} and
                          {wave} and
                          fitted = {self.fitted}''', con)

            assert len(result) <= 1
            if len(result) == 1:
                savefile = result.filename[0]

                with open(savefile, 'rb') as f:
                    model_result = pickle.load(f)
                model_result['idnum'] = result.idnum[0]
                model_result['savefile'] = savefile
                
                # This is a check -- I don't think it will ever happen
                if len(model_result['radiance']) != len(data):
                    model_result = {'radiance': None,
                                    'model_total_source': None,
                                    'weighting': None,
                                    'packets': None,
                                    'idnum': None,
                                    'savefile': None}
                else:
                    pass
            else:
                model_result = {'radiance':None,
                                'model_total_source':None,
                                'weighting':None,
                                'packets':None,
                                'idnum':None,
                                'savefile': None}

        return model_result
    
    def create_model(self, data, output):
        # distance of s/c from planet
        dist_from_plan = np.sqrt(data.x**2 + data.y**2 + data.z**2)

        # Angle between look direction and planet.
        ang = np.arccos((-data.x*data.xbore - data.y*data.ybore -
                         data.z*data.zbore)/dist_from_plan)

        # Check to see if look direction intersects the planet anywhere
        asize_plan = np.arcsin(1./dist_from_plan)

        # Don't worry about lines of sight that don't hit the planet
        dist_from_plan.loc[ang > asize_plan] = 1e30

        # Load the outputfile
        packets = output.X
        packets['radvel_sun'] = (packets['vy'] +
                                 output.vrplanet.to(self.unit/u.s).value)

        # Will base shadow on line of sight, not the packets
        out_of_shadow = np.ones(len(packets))
        self.packet_weighting(packets, out_of_shadow, output.aplanet)

        xcols = ['x', 'y', 'z']
        borecols = ['xbore', 'ybore', 'zbore']

        # This sets limits on regions where packets might be
        tree = KDTree(packets[xcols].values)
        # tree = BallTree(packets[xcols].values)

        # This removes the packets that aren't close to the los
        oedge = output.inputs.options.outeredge*2
        
        # This sets limits on regions where packets might be
        rad = pd.Series(np.zeros(data.shape[0]), index=data.index)

        print(f'{data.shape[0]} spectra taken.')
        for i, spectrum in data.iterrows():
            x_sc = spectrum[xcols].values.astype(float)
            bore = spectrum[borecols].values.astype(float)

            dd = 30  # Furthest distance we need to look
            x_far = x_sc+bore*dd
            while np.linalg.norm(x_far) > oedge:
                dd -= 0.1
                x_far = x_sc+bore*dd
                
            t = [0.05]
            while t[-1] < dd:
                t.append(t[-1] + t[-1]*np.sin(self.dphi))
            t = np.array(t)
            Xbore = x_sc[np.newaxis, :]+bore[np.newaxis, :]*t[:, np.newaxis]

            wid = t*np.sin(self.dphi)
            ind = np.concatenate(tree.query_radius(Xbore, wid))
            ilocs = np.unique(ind).astype(int)
            indicies = packets.iloc[ilocs].index
            subset = packets.loc[indicies]

            xpr = subset[xcols]-x_sc[np.newaxis, :]
            rpr = np.sqrt(xpr['x']*xpr['x'] +
                          xpr['y']*xpr['y'] +
                          xpr['z']*xpr['z'])

            losrad = np.sum(xpr*bore[np.newaxis, :], axis=1)
            inview = rpr < dist_from_plan[i]
            
            if np.any(inview):
                Apix = np.pi*(rpr[inview]*np.sin(self.dphi))**2*(
                    self.unit.to(u.cm))**2
                wtemp = subset.loc[inview, 'weight']/Apix*self.atoms_per_packet
                if self.quantity == 'radiance':
                    # Determine if any packets are in shadow
                    # Projection of packet onto LOS
                    # Point along LOS the packet represents
                    losrad_ = losrad[inview].values
                    hit = (x_sc[np.newaxis,:] +
                           bore[np.newaxis,:] * losrad_[:,np.newaxis])
                    rhohit = np.linalg.norm(hit[:,[0,2]], axis=1)
                    out_of_shadow = (rhohit > 1) | (hit[:,1] < 0)
                    wtemp *= out_of_shadow

                    rad.loc[i] = wtemp.sum()
                else:
                    assert False, 'Other quantities not set up.'
            else:
                pass

            if len(data) > 10:
                ind = data.index.get_loc(i)
                if (ind%(len(data)//10)) == 0:
                    print(f'Completed {ind+1} spectra')

        result = {'radiance': rad,
                  'model_total_source': output.totalsource,
                  'weighting': None,
                  'packets': None}

        return result
    
    def fit_model_to_data(self, data, output, masking=None):
        '''Determine the source distribution that best fits the data'''

        # Helper functions
        def should_add_weight(index, saved):
            return index in saved

        def add_weight(x, ratio):
            return np.append(x, ratio)
        
        def add_index(x, i):
            return np.append(x, i)

        # distance of s/c from planet
        dist_from_plan = np.sqrt(data.x**2+data.y**2+data.z**2)

        # Angle between look direction and planet.
        ang = np.arccos((-data.x*data.xbore-data.y*data.ybore-
                         data.z*data.zbore)/dist_from_plan)

        # Check to see if look direction intersects the planet anywhere
        asize_plan = np.arcsin(1./dist_from_plan)

        # Don't worry about lines of sight that don't hit the planet
        dist_from_plan.loc[ang > asize_plan] = 1e30

        # Load the outputfile
        packets = output.X.copy()
        packets['radvel_sun'] = (packets['vy']+
                                 output.vrplanet.to(self.unit/u.s).value)
    
        # Will base shadow on line of sight, not the packets
        out_of_shadow = np.ones(packets.shape[0])
        self.packet_weighting(packets, out_of_shadow, output.aplanet)
    
        xcols = ['x', 'y', 'z']
        borecols = ['xbore', 'ybore', 'zbore']
    
        # This sets limits on regions where packets might be
        tree = KDTree(packets[xcols].values)
        # tree = BallTree(packets[xcols].values)
    
        # Maxiumum line of sight length
        oedge = output.inputs.options.outeredge*2
    
        # rad = modeled radiance
        # saved_packets = list of indicies of the packets used for each spectrum
        # weighting = list of the weights that should be applied
        #   - Final weighting for each packet is mean of weights
        rad = pd.Series(np.zeros(data.shape[0]), index=data.index)
        saved_packets = pd.Series((np.ndarray((0,), dtype=int)
                                   for _ in range(data.shape[0])),
                                  index=data.index)
        ind0 = packets.Index.unique()
        weighting = pd.Series((np.ndarray((0,)) for _ in range(ind0.shape[0])),
                              index=ind0)
        included = pd.Series((np.ndarray((0,), dtype=np.int)
                              for _ in range(ind0.shape[0])), index=ind0)

        # Determine which points should be used for the fit
        _, _, mask = fit_model(data.radiance, None, data.sigma,
                               masking=masking, mask_only=True,
                               altitude=data.alttan)
    
        print(f'{data.shape[0]} spectra taken.')
        for i, spectrum in data.iterrows():
            x_sc = spectrum[xcols].values.astype(float)
            bore = spectrum[borecols].values.astype(float)
        
            dd = 30  # Furthest distance we need to look
            x_far = x_sc+bore*dd
            while np.linalg.norm(x_far) > oedge:
                dd -= 0.1
                x_far = x_sc+bore*dd
        
            t = [0.05]
            while t[-1] < dd:
                t.append(t[-1]+t[-1]*np.sin(self.dphi))
            t = np.array(t)
            Xbore = x_sc[np.newaxis, :]+bore[np.newaxis, :]*t[:, np.newaxis]
        
            wid = t*np.sin(self.dphi)
            ind = np.concatenate(tree.query_radius(Xbore, wid))
            ilocs = np.unique(ind).astype(int)
            indicies = packets.iloc[ilocs].index
            subset = packets.loc[indicies]
        
            xpr = subset[xcols]-x_sc[np.newaxis, :]
            rpr = np.sqrt(xpr['x']*xpr['x']+
                          xpr['y']*xpr['y']+
                          xpr['z']*xpr['z'])
    
            losrad = np.sum(xpr*bore[np.newaxis, :], axis=1)
            inview = rpr < dist_from_plan[i]
        
            if np.any(inview):
                Apix = np.pi*(rpr[inview]*np.sin(self.dphi))**2*(
                    self.unit.to(u.cm))**2
                wtemp = subset.loc[inview, 'weight']/Apix*self.atoms_per_packet
                if self.quantity == 'radiance':
                    # Determine if any packets are in shadow
                    # Projection of packet onto LOS
                    # Point along LOS the packet represents
                    losrad_ = losrad[inview].values
                    hit = (x_sc[np.newaxis, :]+
                           bore[np.newaxis, :]*losrad_[:, np.newaxis])
                    rhohit = np.linalg.norm(hit[:, [0, 2]], axis=1)
                    out_of_shadow = (rhohit > 1) | (hit[:, 1] < 0)
                    wtemp *= out_of_shadow
                
                    rad.loc[i] = wtemp.sum()
                    if (wtemp.sum() > 0) and mask[i]:
                        ratio = spectrum.radiance/wtemp.sum()
                
                        # Save which packets are used for each spectrum
                        saved_packets[i] = subset.loc[inview, 'Index'].unique()
                        should = weighting.index.to_series().apply(
                                should_add_weight, args=(saved_packets[i],))
                        weighting.loc[should] = weighting.loc[should].apply(
                                add_weight, args=(ratio,))
                        included.loc[should] = included.loc[should].apply(
                                add_index, args=(i,))
                    else:
                        pass
                else:
                    assert False, 'Other quantities not set up.'
            else:
                pass
        
            if len(data) > 10:
                ind = data.index.get_loc(i)
                if (ind%(len(data)//10)) == 0:
                    print(f'Completed {ind+1} spectra')

        # Determine the proper weightings
        assert np.all(weighting.apply(len) == included.apply(len))
        new_weight = weighting.apply(lambda x: x.mean() if x.shape[0] > 0 else 0.)
        new_weight /= new_weight[new_weight > 0].mean()
        # from IPython import embed
        # embed()
        # import sys; sys.exit()

        assert np.all(np.isfinite(new_weight))
        if np.any(new_weight > 0):
            multiplier = new_weight.loc[packets['Index']].values
            output.X.loc[:, 'frac'] = packets.loc[:, 'frac']*multiplier
            output.X0.loc[:, 'frac'] = output.X0.loc[:, 'frac']*new_weight
            
            output.X = output.X[output.X.frac > 0]
            output.X0 = output.X0[output.X0.frac > 0]

            # Update the LOSResult and output objects with new values
            output.totalsource = output.X0['frac'].sum()
            self.totalsource = output.totalsource
            self.mod_rate = self.totalsource/self.inputs.options.endtime.value
            self.atoms_per_packet = 1e23/self.mod_rate
            print('In fiited model:')
            print(f'Total source = {self.totalsource} packets')
            print(f'1 packet represents {self.atoms_per_packet} atoms')
            print(f'Model rate = {self.mod_rate} packets/sec')

            # Run the model with updated source
            result_with_fitted = self.create_model(data, output)
            
            weighting = pd.DataFrame({'weight':weighting.values,
                                      'included':included.values})
            result = {'radiance':result_with_fitted['radiance'],
                      'model_total_source': self.totalsource,
                      'weighting':weighting,
                      'packets':saved_packets}
        else:
            result = {'radiance': pd.Series(np.zeros(data.shape[0]),
                                            index=data.index),
                      'model_total_source': 0,
                      'weighting': None,
                      'packets': None}

        return result
