import os.path
import numpy as np
import pandas as pd
import pickle
import random
import astropy.units as u
from datetime import datetime
from sklearn.neighbors import KDTree
#from scipy.spatial import KDTree, cKDTree
from .ModelResults import ModelResult
from .database_connect import database_connect
from .Input import Input
from .Output import Output


class LOSResult(ModelResult):
    def __init__(self, start_from, data, quantity, dphi=3*u.deg,
                 filenames=None, overwrite=False, savepackets=False,
                 **kwargs):
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
            super().__init__(inputs, format_, filenames=filenames)
        elif isinstance(start_from, Output):
            output = start_from
            inputs = output.inputs
            super().__init__(inputs, format_, output=output)
        else:
            raise TypeError

        tstart = datetime.now()
        self.type = 'LineOfSight'
        self.species = inputs.options.species
        self.origin = inputs.geometry.planet
        self.unit = u.def_unit('R_' + self.origin.object,
                               self.origin.radius)
        self.dphi = dphi.to(u.rad).value

        nspec = len(data)
        self.radiance = np.zeros(nspec)
        self.packets = pd.DataFrame()
        self.ninview = np.zeros(nspec, dtype=int)
        for j,outfile in enumerate(self.filenames):
            if outfile is None:
                radiance_, npackets_, saved_ = self.create_model(data, output,
                                                    savepackets=savepackets)
                print(f'Completed model {j+1} of {len(self.filenames)}')
            else:
                # Search to see if it is already done
                radiance_, npackets_, saved_, idnum = self.restore(data, outfile)

                if (radiance_ is None) or overwrite:
                    if (radiance_ is not None) and overwrite:
                        self.delete_model(idnum)
                    else:
                        pass

                    output = Output.restore(outfile)
                    radiance_, npackets_, saved_ = self.create_model(data,
                                            output, savepackets=savepackets)
                    self.save(data, outfile, radiance_, npackets_, saved_)
                    print(f'Completed model {j+1} of {len(self.filenames)}')
                else:
                    print(f'Model {j+1} of {len(self.filenames)} '
                          'previously completed.')

            self.radiance += radiance_
            self.ninview += npackets_.astype(int)
            self.packets[j] = saved_

        self.radiance *= self.atoms_per_packet
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

    def save(self, data, fname, radiance, packets, saved_packets):
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
                                    orbit, dphi, mechanism, wavelength, filename)
                                    values (%s, %s, %s, %s, %s, %s, %s)''',
                                (idnum, self.quantity, orb, self.dphi,
                                 mech, wave, tempname))

                    # Determine the savefile name
                    idnum_ = pd.read_sql(f'''SELECT idnum
                                             FROM uvvsmodels
                                             WHERE filename='{tempname}';''', con)
                    assert len(idnum_) == 1
                    idnum = int(idnum_.idnum[0])

                    savefile = os.path.join(os.path.dirname(fname),
                                        f'model.orbit{orb:04}.{idnum}.pkl')
                    with open(savefile, 'wb') as f:
                        pickle.dump((radiance, packets, saved_packets), f)
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
                              {wave}''', con)

            assert len(result) <= 1
            if len(result) == 1:
                savefile = result.filename[0]

                with open(savefile, 'rb') as f:
                    radiance, packets, saved_ = pickle.load(f)
                idnum = result.idnum[0]
                if len(radiance) != len(data):
                    radiance, packets, saved_, idnum = None, None, None, None
                else:
                    pass
            else:
                radiance, packets, saved_, idnum = None, None, None, None

        return radiance, packets, saved_, idnum

    def create_model(self, data, output, savepackets, **kwargs):
        # distance of s/c from planet
        dist_from_plan = (np.sqrt(data.x**2 + data.y**2 + data.z**2)).values

        # Angle between look direction and planet.
        ang = np.arccos((-data.x*data.xbore - data.y*data.ybore -
                         data.z*data.zbore)/dist_from_plan)

        # Check to see if look direction intersects the planet anywhere
        asize_plan = np.arcsin(1./dist_from_plan)

        # Don't worry about lines of sight that don't hit the planet
        dist_from_plan[ang > asize_plan] = 1e30

        # Load the outputfile
        packets = output.X
        packets['radvel_sun'] = (packets['vy'] +
                                 output.vrplanet.to(self.unit/u.s).value)

        # Will base shadow on line of sight, not the packets
        out_of_shadow = np.ones(len(packets))
        self.packet_weighting(packets, out_of_shadow, output.aplanet)
        
        xpack = packets[['x', 'y', 'z']].values
        weight = packets['weight'].values
        index = packets['index'].values
        tree = KDTree(xpack)
        # tree = BallTree(xpack)

        # This sets limits on regions where packets might be
        rad, pack = np.zeros(data.shape[0]), np.zeros(data.shape[0])
        if savepackets:
            saved_packets = np.ndarray(data.shape[0], dtype='O')
        else:
            saved_packets = None

        xdata = data[['x', 'y', 'z']].values.astype(float)
        boresight = data[['xbore', 'ybore', 'zbore']].values.astype(float)
        
        # This removes the packets that aren't close to the los
        if 'outeredge' in kwargs:
            oedge = kwargs['outeredge']
        else:
            oedge = output.inputs.options.outeredge*2

        print(f'{data.shape[0]} spectra taken.')
        for i in range(data.shape[0]):
            x_sc = xdata[i,:]
            bore = boresight[i,:]

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
            indicies = np.unique(ind).astype(int)
            subset = xpack[indicies, :]
            wt = weight[indicies]
            subindex = index[indicies]

            xpr = subset-x_sc[np.newaxis, :]
            rpr = np.sqrt(xpr[:, 0]*xpr[:, 0]+xpr[:, 1]*xpr[:, 1]+
                          xpr[:, 2]*xpr[:, 2])
            losrad = np.sum(xpr*bore[np.newaxis, :], axis=1)
            inview = rpr < dist_from_plan[i]
            
            if np.any(inview):
                Apix = np.pi*(rpr[inview]*np.sin(self.dphi))**2 * self.unit**2
                wtemp = wt[inview]/Apix.to(u.cm**2).value
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

                    rad[i] = np.sum(wtemp)
                    pack[i] = len(losrad_)
                    if savepackets:
                        saved_packets[i] = tuple(set(subindex[inview]))
                    else:
                        pass
            else:
                pass

            if len(data) > 10:
                if (i%(len(data)//10)) == 0:
                    print(f'Completed {i+1} spectra')
                    
        del output
        return rad, pack, saved_packets
