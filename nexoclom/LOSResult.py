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


class LOSResult(ModelResult):
    '''Class to contain the LOS result from multiple outputfiles.'''
    def __init__(self, start_from, mdata, quantity, dphi=3*u.deg,
                 filenames=None, overwrite=False, savepackets=False,
                 fit_to_data=False, **kwargs):
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
        if isinstance(start_from, Input):
            inputs = start_from
            output = None
            super().__init__(inputs, format_, filenames=filenames)
        elif isinstance(start_from, Output):
            output = start_from
            inputs = output.inputs
            super().__init__(inputs, format_, output=output)
        else:
            raise TypeError
        
        # Make sure data is in the right reference frame
        mdata.set_frame('Model')
        data = mdata.data

        tstart = datetime.now()
        self.type = 'LineOfSight'
        self.species = inputs.options.species
        self.origin = inputs.geometry.planet
        self.unit = u.def_unit('R_' + self.origin.object,
                               self.origin.radius)
        self.dphi = dphi.to(u.rad).value
        self.fitted = fit_to_data

        nspec = len(data)
        self.radiance = np.zeros(nspec)
        self.packets = pd.DataFrame()
        self.ninview = np.zeros(nspec, dtype=int)
        for j,outfile in enumerate(self.filenames):
            if outfile is None:  # Output was provided
                if self.fitted:
                    masking = kwargs['masking'] if 'masking' in kwargs else None
                    radiance_, npackets_, output = self.fit_model_to_data(
                            data, output, masking=masking)
                else:
                    radiance_, npackets_ = self.create_model(data,
                        output, savepackets=savepackets)
                print(f'Completed model {j+1} of {len(self.filenames)}')
                del output
            else:
                # Search to see if it is already done
                radiance_, npackets_, idnum, output = self.restore(data, outfile)
                if (radiance_ is None) or overwrite:
                    if (radiance_ is not None) and overwrite:
                        self.delete_model(idnum)
                    else:
                        pass

                    output = Output.restore(outfile)
                    if self.fitted:
                        masking = kwargs['masking'] if 'masking' in kwargs else None
                        radiance_, npackets_, output_ = self.fit_model_to_data(
                                data, output, masking=masking)
                        self.save(data, outfile, radiance_, npackets_, output_)
                        del output_
                    else:
                        radiance_, npackets_ = self.create_model(data,
                                            output, savepackets=savepackets)
                        self.save(data, outfile, radiance_, npackets_)
                    print(f'Completed model {j+1} of {len(self.filenames)}')
                    del output
                else:
                    print(f'Model {j+1} of {len(self.filenames)} '
                          'previously completed.')

            self.radiance += radiance_.values
            self.ninview += npackets_.astype(int).values

        self.radiance *= u.R
        tend = datetime.now()
        print(f'Total time = {tend-tstart}')

    def delete_model(self, idnum):
        with database_connect() as con:
            cur = con.cursor()
            cur.execute('''SELECT idnum, filename FROM uvvsmodels
                           WHERE out_idnum = %s''', (int(idnum), ))
            assert cur.rowcount in (0, 1)
            for mid, mfile in cur.fetchall():
                cur.execute('''DELETE from uvvsmodels
                               WHERE idnum = %s''', (mid, ))
                if os.path.exists(mfile):
                    os.remove(mfile)

    def save(self, data, fname, radiance, packets, output=None):
        # Determine if the model can be saved.
        # Criteria: 1 complete orbit, nothing more.
        orbits = set(data.orbit)
        orb = orbits.pop()

        if len(orbits) != 0:
            print('Model spans more than one orbit. Cannot be saved.')
        else:
            from MESSENGERuvvs import MESSENGERdata
            mdata = MESSENGERdata(self.species, f'orbit = {orb}')
            if len(mdata) != len(data):
                print('Model does not contain the complete orbit. '
                      'Cannot be saved.')
            else:
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

                    tempname = f'temp_{orb}_{str(random.randint(0, 1000000))}'
                    cur.execute(f'''INSERT into uvvsmodels (out_idnum, quantity,
                                    orbit, dphi, mechanism, wavelength,
                                    fitted, filename)
                                    values (%s, %s, %s, %s, %s, %s, %s, %s)''',
                                (idnum, self.quantity, orb, self.dphi,
                                 mech, wave, self.fitted, tempname))

                    # Determine the savefile name
                    idnum_ = pd.read_sql(f'''SELECT idnum
                                             FROM uvvsmodels
                                             WHERE filename='{tempname}';''', con)
                    assert len(idnum_) == 1
                    idnum = int(idnum_.idnum[0])

                    savefile = os.path.join(os.path.dirname(fname),
                                        f'model.orbit{orb:04}.{idnum}.pkl')
                    if self.fitted:
                        with open(savefile, 'wb') as f:
                            pickle.dump((radiance, packets, output), f)
                    else:
                        with open(savefile, 'wb') as f:
                            pickle.dump((radiance, packets), f)
                    cur.execute(f'''UPDATE uvvsmodels
                                    SET filename=%s
                                    WHERE idnum=%s''', (savefile, idnum))

    def restore(self, data, fname):
        # Determine if the model can be restored.
        # Criteria: 1 complete orbit, nothing more.
        orbits = set(data.orbit)
        orb = orbits.pop()

        if len(orbits) != 0:
            print('Model spans more than one orbit. Cannot be saved.')
            radiance, packets, idnum = None, None, None
        else:
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

                result = pd.read_sql(
                    f'''SELECT idnum, filename FROM uvvsmodels
                        WHERE out_idnum={oid} and
                              quantity = '{self.quantity}' and
                              orbit = {orb} and
                              dphi = {self.dphi} and
                              {mech} and
                              {wave} and
                              fitted = {self.fitted}''', con)

            assert len(result) <= 1
            if len(result) == 1:
                savefile = result.filename[0]

                if self.fitted:
                    with open(savefile, 'rb') as f:
                        radiance, packets, output = pickle.load(f)
                else:
                    with open(savefile, 'rb') as f:
                        radiance, packets = pickle.load(f)
                    output = None
                idnum = result.idnum[0]
                
                if len(radiance) != len(data):
                    radiance, packets, idnum, output = None, None, None, None
                else:
                    pass
            else:
                radiance, packets, idnum, output = None, None, None, None

        return radiance, packets, idnum, output
    
    def create_model(self, data, output, savepackets=False, **kwargs):
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
        npack = pd.Series(np.zeros(data.shape[0]), index=data.index)

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
                    losrad_ = losrad[inview]
                    # Determine if any packets are in shadow
                    # Projection of packet onto LOS
                    # Point along LOS the packet represents
                    hit = (x_sc[np.newaxis,:] +
                           bore[np.newaxis,:] * losrad_[:,np.newaxis])
                    rhohit = np.linalg.norm(hit[:,[0,2]], axis=1)
                    out_of_shadow = (rhohit > 1) | (hit[:,1] < 0)
                    wtemp *= out_of_shadow

                    rad.loc[i] = wtemp.sum()
                    npack.loc[i] = sum(inview)
                else:
                    assert False, 'Other quantities not set up.'
            else:
                pass

            if len(data) > 10:
                ind = data.index.get_loc(i)
                if (ind%(len(data)//10)) == 0:
                    print(f'Completed {ind+1} spectra')

        return rad, npack
    
    def fit_model_to_data(self, data, output, masking=None, makeplots=False):

        def should_add_weight(index, saved):
            return index in saved

        def add_weight(x, ratio):
            return np.append(x, ratio)

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
        packets = output.X
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
        # npack = number of packets used to compute radiance
        # saved_packets = list of indicies of the packets used for each spectrum
        # weighting = list of the weights that should be applied
        #   - Final weighting for each packet is mean of weights
        rad = pd.Series(np.zeros(data.shape[0]), index=data.index)
        npack = pd.Series(np.zeros(data.shape[0]), index=data.index)
        saved_packets = pd.Series((np.ndarray((0,), dtype=int)
                                   for _ in range(data.shape[0])),
                                  index=data.index)
        ind0 = packets.Index.unique()
        weighting = pd.Series((np.ndarray((0,)) for _ in range(ind0.shape[0])),
                              index=ind0)
    
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
                    npack.loc[i] = sum(inview)
                    if (wtemp.sum() > 0) and mask[i]:
                        ratio = spectrum.radiance/wtemp.sum()
                
                        # Save which packets are used for each spectrum
                        saved_packets[i] = subset.loc[inview, 'Index'].unique()
                        should = weighting.index.to_series().apply(
                                should_add_weight, args=(saved_packets[i],))
                        weighting.loc[should] = weighting.loc[should].apply(
                                add_weight, args=(ratio,))
                    else:
                        pass
            
                    # lon, lat = output.X0[should].longitude, output.X0[should].latitude
                    # plot_radiance_source(data, i, lon*180/np.pi, lat*180/np.pi)
                else:
                    assert False, 'Other quantities not set up.'
            else:
                pass
        
            if len(data) > 10:
                ind = data.index.get_loc(i)
                if (ind%(len(data)//10)) == 0:
                    print(f'Completed {ind+1} spectra')

        # Determine the proper weightings
        new_weight = weighting.apply(lambda x:x.mean() if x.shape[0] > 0 else 0.)
        assert np.all(np.isfinite(new_weight))
        if np.any(new_weight > 0):
            multiplier = new_weight.loc[packets['Index']].values
            output.X.loc[:, 'frac'] = packets.loc[:, 'frac']*multiplier
            output.X0.loc[:, 'frac'] = output.X0.loc[:, 'frac']*new_weight
            output.totalsource = np.sum(output.X0.frac)
            
            rad, npack = self.create_model(data, output)
        else:
            rad = pd.Series(np.zeros(data.shape[0]), index=data.index)
            npack = pd.Series(np.zeros(data.shape[0]), index=data.index)

        return rad, npack, output

