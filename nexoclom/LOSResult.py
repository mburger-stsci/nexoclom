import os.path
import numpy as np
import pandas as pd
import pickle
import astropy.units as u
from MESSENGERuvvs import MESSENGERdata
from .ModelResults import ModelResult
from .database_connect import database_connect
from .Output import Output


class LOSResult(ModelResult):
    def __init__(self, inputs, data, quantity, dphi=3*u.deg,
                 filenames=None, overwrite=False):
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
        format = {'quantity':quantity}
        super().__init__(inputs, format, filenames=filenames)
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

    @staticmethod
    def delete_model(self, idnum):
        with database_connect() as con:
            cur = con.cursor()
            cur.execute('''SELECT idnum, filename FROM uvvsmodels
                           WHERE out_idnum = %s''', (idnum, ))
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
            mdata = MESSENGERdata(self.species, f'orbit = {orb}')
            if len(mdata) != len(data):
                print('Model does not contain the complete orbit. '
                      'Cannot be saved.')
            else:
                con = database_connect()
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

                cur.execute(f'''INSERT into uvvsmodels (out_idnum, quantity,
                                orbit, dphi, mechanism, wavelength, filename)
                                values (%s, %s, %s, %s, %s, %s, 'temp')''',
                            (idnum, self.quantity, orb, self.dphi,
                             mech, wave))

                # Determine the savefile name
                idnum_ = pd.read_sql('''SELECT idnum
                                        FROM uvvsmodels
                                        WHERE filename='temp';''', con)
                assert len(idnum_) == 1
                idnum = int(idnum_.idnum[0])

                savefile = os.path.join(os.path.dirname(fname),
                                        f'model.orbit{orb:04}.{idnum}.pkl')
                with open(savefile, 'wb') as f:
                    pickle.dump((radiance, packets), f)
                cur.execute(f'''UPDATE uvvsmodels
                                SET filename=%s
                                WHERE idnum=%s''', (savefile, idnum))
                con.close()

    def restore(self, data, fname):
        # Determine if the model can be restored.
        # Criteria: 1 complete orbit, nothing more.
        orbits = set(data.orbit)
        orb = orbits.pop()

        if len(orbits) != 0:
            print('Model spans more than one orbit. Cannot be saved.')
            radiance, packets = None, None
        else:
            mdata = MESSENGERdata(self.species, f'orbit = {orb}')
            if len(mdata) != len(data):
                print('Model does not contain the complete orbit. '
                      'Cannot be saved.')
                radiance, packets, idnum = None, None, None
            else:
                con = database_connect()
                con.autocommit = True

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
                else:
                    radiance, packets, idnum = None, None, None

            return radiance, packets, idnum

    def create_model(self, data, outfile):
        # distance of s/c from planet
        
        dist_from_plan = np.sqrt(data.x**2 + data.y**2 + data.z**2)

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

        # This sets limits on regions where packets might be
        xx_, yy_, zz_ = (np.zeros((2,len(data))), np.zeros((2,len(data))),
                         np.zeros((2,len(data))))
        oedge = output.inputs.options.outeredge
        xx_[1,:], yy_[1,:], zz_[1,:] = (data.xbore*oedge*2,
                                        data.ybore*oedge*2,
                                        data.zbore*oedge*2)
        xx = (data.x[np.newaxis,:] + xx_)
        yy = (data.y[np.newaxis,:] + yy_)
        zz = (data.z[np.newaxis,:] + zz_)
        
        box_size = output.inputs.options.outeredge*2*np.sin(self.dphi) + 0.1
        
        xx_min = np.min(xx-box_size, axis=0)
        yy_min = np.min(yy-box_size, axis=0)
        zz_min = np.min(zz-box_size, axis=0)
        xx_max = np.max(xx+box_size, axis=0)
        yy_max = np.max(yy+box_size, axis=0)
        zz_max = np.max(zz+box_size, axis=0)

        rad, pack = np.zeros(len(data)), np.zeros(len(data))

        index = np.arange(len(data))
        for j, row in zip(index,data.iterrows()):
            i, row = row
            # This removes the packets that aren't close to the los
            mask = ((packets.x >= xx_min[j]) &
                    (packets.x <= xx_max[j]) &
                    (packets.y >= yy_min[j]) &
                    (packets.y <= yy_max[j]) &
                    (packets.z >= zz_min[j]) &
                    (packets.z <= zz_max[j]))
            subset = packets[mask]
            
            # Distance of packets from spacecraft
            row_pr = row[['x', 'y', 'z']].values.astype(float)
            xpr = subset[['x', 'y', 'z']].values - row_pr[np.newaxis,:]
            bore = row[['xbore', 'ybore', 'zbore']].values.astype(float)
            rpr = np.linalg.norm(xpr, axis=1)
            
            # Packet-s/c boresight angle
            costheta = np.sum(xpr * bore[np.newaxis,:], axis=1) / rpr
            costheta[costheta > 1] = 1.
            costheta[costheta < -1] = -1.
            
            inview = ((costheta >= np.cos(self.dphi)) &
                      (subset['weight'] > 0) &
                      (rpr < dist_from_plan[i]))

            if np.any(inview):
                Apix = np.pi*(rpr[inview]*np.sin(self.dphi))**2 * self.unit**2
                wtemp = subset.weight[inview]/Apix.to(u.cm**2).value
                if self.quantity == 'radiance':
                    # Determine if any packets are in shadow
                    # Projection of packet onto LOS
                    losrad = rpr[inview] * costheta[inview]

                    # Point along LOS the packet represents
                    hit = (row_pr[np.newaxis,:] +
                           bore[np.newaxis,:] * losrad[:,np.newaxis])
                    rhohit = np.linalg.norm(hit[:,[0,1]], axis=1)
                    out_of_shadow = (rhohit > 1) | (hit[:,1] < 0)
                    wtemp *= out_of_shadow

                    rad[j] = np.sum(wtemp)
                    pack[j] = np.sum(inview)

        del output
        self.save(data, outfile, rad, pack)

        return rad, pack
