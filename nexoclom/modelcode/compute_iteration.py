import numpy as np
import pandas as pd
import astropy.units as u
from sklearn.neighbors import KDTree
from nexoclom.modelcode.Output import Output
from nexoclom.modelcode.ModelResult import IterationResult


def compute_iteration(self, outputfile, scdata):
    xcols = ['x', 'y', 'z']
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

    packets = output.X
    packets['radvel_sun'] = (packets['vy'] +
                             output.vrplanet.to(self.unit / u.s).value)
    aplanet = output.aplanet
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
        
            rad.loc[i] = rad_
            npack.loc[i] = npack_
    
        if len(data) > 10:
            ind = data.index.get_loc(i)
            if (ind % (len(data) // 10)) == 0:
                print(f'Completed {ind + 1} spectra')

    iteration_ = {'radiance': rad,
                  'npackets': npack,
                  'totalsource': output.totalsource,
                  'outputfile': outputfile,
                  'out_idnum': output.idnum,
                  'query': scdata.query,
                  'used': used,
                  'used0': used0}
    iteration_result = IterationResult(iteration_)
    modelfile = self.save(iteration_result)
    iteration_result.modelfile = modelfile

    del packets
    return iteration_result
