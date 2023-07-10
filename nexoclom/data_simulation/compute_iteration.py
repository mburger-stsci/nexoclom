import os
import numpy as np
import pandas as pd
import pickle
import astropy.units as u
from sklearn.neighbors import KDTree
import sqlalchemy as sqla
import sqlalchemy.dialects.postgresql as pg
from datetime import datetime
import time
from nexoclom.particle_tracking.Output import Output
from nexoclom import engine


class IterationResult:
    def __init__(self, iteration, losresult):
        self.radiance = iteration['radiance']
        self.npackets = iteration['npackets']
        self.totalsource = iteration['totalsource']
        self.outputfile = iteration['outputfile']
        self.out_idnum = iteration['out_idnum']
        self.included = iteration['included']
        
        self.modelfile = None
        self.model_idnum = None
        self.fitted = False
        self.used_packets = iteration.get('used', None)
        self.used_packets0 = iteration.get('used0', None)

        self.quantity = losresult.quantity
        self.query = losresult.query
        self.dphi = losresult.dphi
        self.mechanism = losresult.mechanism
        self.wavelength = losresult.wavelength
        self.fitted = losresult.fitted

    def save_iteration(self):
        '''
        Insert the result of a LOS iteration into the database
        :param iteration_result: LOS result from a single outputfile
        :return: name of saved file
        '''
        # Insert the result into the database
        metadata_obj = sqla.MetaData()
        table = sqla.Table("uvvsmodels", metadata_obj, autoload_with=engine)
        
        ufit_id = (self.unfit_outid
                   if isinstance(self, IterationResultFitted)
                   else None)
        
        insert_stmt = pg.insert(table).values(
            out_idnum=self.out_idnum,
            unfit_idnum=ufit_id,
            quantity=self.quantity,
            query=self.query,
            dphi=self.dphi,
            mechanism=self.mechanism,
            wavelength=[w.value for w in self.wavelength],
            fitted=self.fitted)
        
        with engine.connect() as con:
            result = con.execute(insert_stmt)
            con.commit()
            
        self.idnum = result.inserted_primary_key[0]
        savefile = os.path.join(os.path.dirname(self.outputfile),
                                f'model.{self.idnum}.pkl')
        self.modelfile = savefile
        print(f'Saving modelfile: {self.modelfile}')
        update = sqla.update(table).where(table.columns.idnum == self.idnum).values(
            filename=savefile)
        with engine.connect() as con:
            con.execute(update)
            con.commit()

        with open(savefile, 'wb') as f:
            pickle.dump(self, f)


class IterationResultFitted(IterationResult):
    def __init__(self, iteration, losresult):
        super().__init__(iteration, losresult)
        
        self.unfit_outputfile = iteration['unfit_outputfile']
        self.unfit_outid = iteration['unfit_outid']
        self.unfit_modelfile = iteration['unfit_modelfile']
        self.fitted = True


def compute_iteration(self, outputfile, scdata, delay=False):
    t0 = datetime.today()
    
    if delay:
        time.sleep(np.random.random()*10)
    else:
        pass
    
    xcols= ['x', 'y', 'z']
    borecols = ['xbore', 'ybore', 'zbore']
    
    # distance of s/c from planet
    # This is used to determine if the line of sight needs to be cut
    # short because it intersects the planet.
    data = scdata.data
    dist_from_plan = np.sqrt(data.x**2 + data.y**2 + data.z**2)

    # Angle between look direction and planet.
    ang = np.arccos((-data.x * data.xbore - data.y * data.ybore -
                     data.z * data.zbore) / dist_from_plan)

    # Check to see if look direction intersects the planet anywhere
    asize_plan = np.arcsin(1. / dist_from_plan)

    # Don't worry about lines of sight that don't hit the planet
    dist_from_plan.loc[ang > asize_plan] = 1e30

    # simulate the data
    output = Output.restore(outputfile)

    X0_index = output.X0.index
    packets = output.X
    if 'Index' not in packets.columns:
        packets['Index'] = list(packets.index)
    if isinstance(output.vrplanet, float):
        output.vrplanet *= output.unit/u.s
    packets['radvel_sun'] = (packets['vy'] +
                             output.vrplanet.to(self.unit / u.s).value)
    aplanet = output.aplanet
    totalsource = output.totalsource
    idnum = output.idnum
    del output

    # Note: A packet is in shadow if the line-of-sight it is on is
    #       in shadow. This is because the cone used is larger than
    #       the slit.

    # This sets limits on regions where packets might be
    tree = KDTree(packets[xcols].values)

    rad = pd.Series(np.zeros(data.shape[0]), index=data.index)
    npack = pd.Series(np.zeros(data.shape[0]), index=data.index,
                      dtype=int)
    used = pd.Series([set() for _ in range(data.shape[0])], index=data.index)
    used0 = pd.Series([set() for _ in range(data.shape[0])], index=data.index)
    included = pd.Series(index=X0_index, dtype=bool)
    included[:] = False

    print(f'{data.shape[0]} spectra taken.')
    for i, spectrum in data.iterrows():
        x_sc = spectrum[xcols].values.astype(float)
        bore = spectrum[borecols].values.astype(float)
    
        # Distance from spacecraft to edge of field of view
        a = 1
        b = 2*np.sum(x_sc*bore)
        c = np.linalg.norm(x_sc)**2 - self.inputs.options.outeredge**2
        dd = (-b + np.sqrt(b**2 - 4*a*c))/2
    
        # Compute coordinates of the LOS spaced farther apart the farther out
        t = [np.sin(self.dphi)]
        while t[-1] < dd:
            t.append(t[-1] + t[-1] * np.sin(self.dphi))
        t = np.array(t)
        Xbore = x_sc[np.newaxis, :] + bore[np.newaxis, :] * t[:, np.newaxis]
    
        # Narrow down number of packets
        wid = t * np.sin(self.dphi*2)
        ind = np.concatenate(tree.query_radius(Xbore, wid))
        ilocs = np.unique(ind).astype(int)
    
        subset = packets.iloc[ilocs]
        subset_rel_sc = subset[xcols].values - x_sc[np.newaxis, :]
        subset_dist_sc = np.linalg.norm(subset_rel_sc, axis=1)
        losrad = np.sum(subset_rel_sc * bore[np.newaxis, :], axis=1)
        cosang = np.sum(subset_rel_sc * bore[np.newaxis, :], axis=1)/subset_dist_sc
        cosang[cosang > 1] = 1
        ang = np.arccos(cosang)
        assert np.all(np.isfinite(ang))
    
        # Projection of packet onto line of sight
        inview = (losrad < dist_from_plan.loc[i]) & (ang <= self.dphi)
    
        if np.any(inview):
            subset = subset.loc[inview]
            subset_dist_sc = subset_dist_sc[inview]
            losrad = losrad[inview]
            included[subset.Index] = True

            self.packet_weighting(subset, aplanet)
            Apix = np.pi * (subset_dist_sc * np.sin(self.dphi))**2 * (
                self.unit.to(u.cm))**2
            wtemp = subset['weight'] / Apix
        
            if self.quantity == 'radiance':
                # Determine if any packets are in shadow
                # Projection of packet onto LOS
                # Point along LOS the packet represents
                hit = (x_sc[np.newaxis, :] +
                       bore[np.newaxis, :] * losrad[:, np.newaxis])
                rhohit = np.linalg.norm(hit[:, [0, 2]], axis=1)
                out_of_shadow = (rhohit > 1) | (hit[:, 1] < 0)
                wtemp *= out_of_shadow
            
                rad_ = wtemp.sum()
                npack_ = np.sum(inview)
                used.loc[i] = set(subset.loc[wtemp > 0].index)
                used0.loc[i] = set(subset.loc[wtemp > 0, 'Index'])
            else:
                assert False, 'Other quantities not set up.'
            assert np.isfinite(rad_)
        
            rad.loc[i] = rad_
            npack.loc[i] = npack_
    
        if len(data) > 10:
            ind = data.index.get_loc(i)
            if (ind % (len(data) // 10)) == 0:
                print(f'Completed {ind + 1} spectra')

    iteration_ = {'radiance': rad,
                  'npackets': npack,
                  'totalsource': totalsource,
                  'outputfile': outputfile,
                  'out_idnum': idnum,
                  'query': scdata.query,
                  'used': used,
                  'used0': used0,
                  'included': included}
    iteration_result = IterationResult(iteration_, self)
    iteration_result.save_iteration()
    
    t1 = datetime.today()
    print(f'Iteration time: {(t1-t0).seconds/60} minutes')

    del packets
    return iteration_result
