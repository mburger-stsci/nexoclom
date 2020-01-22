import os.path
import numpy as np
# from scipy.spatial import distance_matrix
import pandas as pd
import pickle
import random
import astropy.units as u
from datetime import datetime
from sklearn.neighbors import KDTree, BallTree
from .ModelResults import ModelResult
from .database_connect import database_connect
from .Output import Output


class LOSResult(ModelResult):
    def __init__(self, inputs, data, quantity, dphi=3*u.deg,
                 filenames=None, overwrite=False, **kwargs):
        """Determine column or emission along lines of sight.
        This assumes the model has already been run.
        
        Parameters
        ==========
        inputs
            An Inputs object
        
        data
            A Pandas DataFrame object with information on the lines of sight.
            
        quantity
            Quantity to calculate: 'column', 'radiance', 'density'
            
        dphi
            Angular size of the view cone. Default = 3 deg.
            
        filenames
            A filename or list of filenames to use. Default = None is to
            find all files created for the inputs.
            
        overwrite
            If True, deletes any images that have already been computed.
            Default = False
        """
        format_ = {'quantity':quantity}
        super().__init__(inputs, format_, filenames=filenames)
        
        tstart = datetime.now()
        self.type = 'LineOfSight'
        self.species = inputs.options.species
        self.origin = inputs.geometry.planet
        self.unit = u.def_unit('R_' + self.origin.object,
                               self.origin.radius)
        self.dphi = dphi.to(u.rad).value

        nspec = len(data)
        self.radiance = np.zeros(nspec)
        self.packets = np.zeros(nspec)
        self.ninview = np.zeros(nspec, dtype=int)

        for j,outfile in enumerate(self.filenames):
            # Search to see if it is already done
            radiance_, packets_, idnum = self.restore(data, outfile)

            if (radiance_ is None) or overwrite:
                if (radiance_ is not None) and overwrite:
                    self.delete_model(idnum)
                radiance_, packets_, = self.create_model(data, outfile)
                print(f'Completed model {j+1} of {len(self.filenames)}')
            else:
                print(f'Model {j+1} of {len(self.filenames)} '
                      'previously completed.')

            self.radiance += radiance_
            self.packets += packets_

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

    def save(self, data, fname, radiance, packets):
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
                              {wave}''', con)

            assert len(result) <= 1
            if len(result) == 1:
                savefile = result.filename[0]

                with open(savefile, 'rb') as f:
                    radiance, packets = pickle.load(f)
                idnum = result.idnum[0]
                if len(radiance) != len(data):
                    radiance, packets, idnum = None, None, None
                else:
                    pass
            else:
                radiance, packets, idnum = None, None, None

        return radiance, packets, idnum

    def create_model(self, data, outfile, **kwargs):
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
        output = Output.restore(outfile)
        packets = output.X
        packets['radvel_sun'] = (packets['vy'] +
                                 output.vrplanet.to(self.unit/u.s).value)

        # Will base shadow on line of sight, not the packets
        out_of_shadow = np.ones(len(packets))
        self.packet_weighting(packets, out_of_shadow, output.aplanet)
        
        xpack = packets[['x', 'y', 'z']].values
        weight = packets['weight'].values
        ta_ = datetime.today()
        tree = BallTree(xpack)
        tb_ = datetime.today()
        print((tb_-ta_).microseconds)
        
        from IPython import embed
        embed() # drop into an IPython session.

        assert 0

        # This sets limits on regions where packets might be
        rad, pack = np.zeros(len(data)), np.zeros(len(data))
        rad_, pack_ = np.zeros(len(data)), np.zeros(len(data))

        xdata = data[['x', 'y', 'z']].values.astype(float)
        boresight = data[['xbore', 'ybore', 'zbore']].values.astype(float)
        
        # This removes the packets that aren't close to the los
        if 'outeredge' in kwargs:
            oedge = kwargs['outeredge']
        else:
            oedge = output.inputs.options.outeredge*2

        print(f'{len(data)} spectra taken.')
        t0, t1, t2, t3 = 0, 0, 0, 0
        for i in range(len(data)):
            x_sc = xdata[i,:]
            bore = boresight[i,:]

            dd = 30  # Furthest distance we need to look
            x_far = x_sc+bore*dd
            while np.linalg.norm(x_far) > oedge:
                dd -= 0.1
                x_far = x_sc+bore*dd
                
            ###########################
            t0_ = datetime.today()
            t = [0.05]
            while t[-1] < dd:
                t.append(t[-1] + t[-1]*np.sin(self.dphi))
            t =  np.array(t)
            Xbore = x_sc[np.newaxis, :]+bore[np.newaxis, :]*t[:, np.newaxis]

            wid = t*np.sin(self.dphi)
            ind = np.concatenate(tree.query_radius(Xbore, wid))
            indicies = np.unique(ind)
            subset = xpack[indicies, :]
            wt_ = weight[indicies]
            
            xpr = subset-x_sc[np.newaxis, :]
            # rpr_ = np.sqrt(xpr[:,0]*xpr[:,0] + xpr[:,1]*xpr[:,1] +
            #               xpr[:,2]*xpr[:,2])
            losrad_ = np.sum(xpr*bore[np.newaxis, :], axis=1)

            #############################
            t1_ = datetime.today()
            x_far = x_sc + bore*oedge
            b_min = np.minimum(x_sc, x_far)
            b_max = np.maximum(x_sc, x_far)
            
            wid = dd*np.sin(self.dphi + 0.5*np.pi/180)

            mask = ((xpack[:,0] >= b_min[0]-wid) &
                    (xpack[:,0] <= b_max[0]+wid) &
                    (xpack[:,1] >= b_min[1]-wid) &
                    (xpack[:,1] <= b_max[1]+wid) &
                    (xpack[:,2] >= b_min[2]-wid) &
                    (xpack[:,2] <= b_max[2]+wid)).nonzero()[0]
            subset = xpack[mask]
            wt = weight[mask]

            # Distance of packets from spacecraft
            xpr = subset - x_sc[np.newaxis,:]
            rpr = np.sqrt(xpr[:,0]*xpr[:,0] + xpr[:,1]*xpr[:,1] +
                          xpr[:,2]*xpr[:,2])

            # Packet-s/c boresight angle
            losrad = np.sum(xpr * bore[np.newaxis,:], axis=1)
            costheta = losrad/rpr
            costheta[costheta > 1] = 1.
            costheta[costheta < -1] = -1.

            inview = ((costheta >= np.cos(self.dphi)) &
                      (rpr < dist_from_plan[i]))
            t1a_ = datetime.today()
            #############################
            if np.any(inview):
                t2_ = datetime.today()
                Apix_ = np.pi*(losrad_*np.sin(self.dphi))**2 * self.unit**2
                wtemp_ = wt_/Apix_.to(u.cm**2).value
                if self.quantity == 'radiance':
                    # Determine if any packets are in shadow
                    # Projection of packet onto LOS
                    # Point along LOS the packet represents
                    hit = (x_sc[np.newaxis,:] +
                           bore[np.newaxis,:] * losrad_[:,np.newaxis])
                    rhohit_ = np.linalg.norm(hit[:,[0,2]], axis=1)
                    out_of_shadow_ = (rhohit_ > 1) | (hit[:,1] < 0)
                    wtemp_ *= out_of_shadow_

                    rad_[i] = np.sum(wtemp_)
                    pack_[i] = len(losrad_)
                #########################
                t3_ = datetime.now()
                losrad_ = losrad[inview]
                Apix = np.pi*(losrad_*np.sin(self.dphi))**2 * self.unit**2
                wtemp = wt[inview]/Apix.to(u.cm**2).value
                if self.quantity == 'radiance':
                    # Determine if any packets are in shadow
                    # Point along LOS the packet represents
                    hit = (x_sc[np.newaxis,:] +
                           bore[np.newaxis,:] * losrad_[:,np.newaxis])
                    rhohit = np.linalg.norm(hit[:,[0,2]], axis=1)
                    out_of_shadow = (rhohit > 1) | (hit[:,1] < 0)
                    wtemp *= out_of_shadow
    
                    rad[i] = np.sum(wtemp)
                    pack[i] = np.sum(inview)
                t4_ = datetime.now()
            else:
                pass

            if len(data) > 10:
                if (i % (len(data)//10)) == 0:
                    print(f'Completed {i+1} spectra')
            t0 += (t1_-t0_).microseconds
            t1 += (t1a_-t1_).microseconds
            t2 += (t3_-t2_).microseconds
            t3 += (t4_-t3_).microseconds

        print(t0/1e6, t1/1e6, t2/1e6, t3/1e6)
        import matplotlib.pyplot as plt
        plt.plot(rad_)
        plt.plot(rad)
        plt.show()
        assert 0
        # del output
        self.save(data, outfile, rad, pack)

        return rad, pack
