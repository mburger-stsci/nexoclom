import numpy as np
import pandas as pd
import astropy.units as u
from sklearn.neighbors import BallTree
from nexoclom.modelcode.Output import Output
import nexoclom.math as mathMB


def make_source_map(outputfile, params, todo=None):
    """
    At each point in lon/lat grid want:
        * Source flux (atoms/cm2/s
        * Speed distribution (f_v vs v)
        * Azimuthal distribution (f_az vs az) -> measured CCW from north
        * Altitude distribution (f_alt vs alt) -> tangent = 0, normal = 90
        :param smear_radius:
        :param nlonbins:
        :param nlatbins:
        :param nvelbins:
        :param nazbins:
        :param naltbins:
        :param use_condor:
        :param normalize:
        :param do_source:
        :param do_available:
        :return:
    """
    
    if todo == 'source':
        print('Determining modeled source')
    elif todo == 'available':
        print('Determining available source')
    else:
        return None
    
    if params is None:
        params = {}
    else:
        pass
    
    smear_radius = params.get('smear_radius', np.radians(10))
    nlonbins = params.get('nlonbins', 180)
    nlatbins = params.get('nlatbins', 90)
    nvelbins = params.get('nvelbins', 100)
    nazbins = params.get('nazbins', 45)
    naltbins = params.get('naltbins', 23)
    smear_abundance = params.get('smear_abundance', True)
    
    print(outputfile)
    output = Output.restore(outputfile)
    X0 = output.X0
    R_planet = output.inputs.geometry.planet.radius.to(u.km)
    del output
    
    vmax = np.ceil(X0['v'].max() * R_planet.value)
    
    included = X0.frac > 0
    if todo == 'source':
        weight = X0.frac
    elif todo == 'available':
        weight = pd.Series(index=X0.index, dtype=float)
        weight[:] = 1
    else:
        assert False

    # List of all X0 near the point on the surface
    # Calculate the histograms and available fraction
    
    # Need surface distribution longitude, latitude for BallTree
    abundance = mathMB.Histogram2d(X0.loc[included, 'longitude'],
                                   X0.loc[included, 'latitude'],
                                   weights=weight[included],
                                   range=[[0, 2 * np.pi],
                                          [-np.pi / 2, np.pi / 2]],
                                   bins=(nlonbins, nlatbins))
    gridlatitude, gridlongitude = np.meshgrid(abundance.y,
                                              abundance.x)
    # If not normalized, leave unitless. Hard to know what units
    # should be
    distribution = {}
    distribution['abundance_uncor'] = abundance.histogram
    distribution['longitude'] = abundance.x * u.rad
    distribution['latitude'] = abundance.y * u.rad
    
    # Speed distribution - full planet
    velocity = mathMB.Histogram(X0.loc[included, 'v'] * R_planet.value,
                                bins=nvelbins, range=[0, vmax],
                                weights=weight[included])
    distribution['speed_dist'] = velocity.histogram
    distribution['speed'] = velocity.x * u.km / u.s
    
    # Altitude distribution
    altitude = mathMB.Histogram(X0.loc[included, 'altitude'],
                                bins=naltbins,
                                range=[0, np.pi / 2],
                                weights=weight[included])
    distribution['altitude_dist'] = altitude.histogram
    distribution['altitude'] = altitude.x * u.rad
    
    # Azimuth distribution
    azimuth = mathMB.Histogram(X0.loc[included, 'azimuth'],
                               bins=nazbins,
                               range=[0, 2 * np.pi],
                               weights=weight[included])
    distribution['azimuth_dist'] = azimuth.histogram
    distribution['azimuth'] = azimuth.x * u.rad
    
    # Lon, lat as list of points
    points = np.array([gridlatitude.flatten(),
                       gridlongitude.flatten()]).T
    tree = BallTree(X0[['latitude', 'longitude']], metric='haversine')
    ind = tree.query_radius(points, smear_radius*np.cos(points[:,0]))

    # Determine spatial variations over the surface
    n_included = np.zeros((points.shape[0],))  # X0 seen by UVVS
    n_total = np.zeros((points.shape[0],))  # all X0
    v_point = np.zeros((points.shape[0], nvelbins))  # spd dist
    alt_point = np.zeros((points.shape[0], naltbins))  # alt dist
    az_point = np.zeros((points.shape[0], nazbins))  # az dist
    if smear_abundance:
        abundance = np.zeros((points.shape[0], ))
    else:
        pass
    for index in range(points.shape[0]):
        if len(ind[index]) > 0:
            subset = X0.loc[ind[index]]
            sub_incl = included[ind[index]]
            sub_weight = weight[ind[index]]
            
            n_included[index] = sub_incl.sum()
            n_total[index] = len(subset)
            
            if smear_abundance:
                abundance[index] = sub_weight.sum()
            else:
                pass
            
            vpoint_ = mathMB.Histogram(subset.loc[sub_incl, 'v'] * R_planet.value,
                                       bins=nvelbins,
                                       range=[0, vmax],
                                       weights=sub_weight[sub_incl])
            v_point[index, :] = vpoint_.histogram
            
            altpoint_ = mathMB.Histogram(subset.loc[sub_incl, 'altitude'],
                                         bins=naltbins,
                                         range=[0, np.pi / 2],
                                         weights=sub_weight[sub_incl])
            alt_point[index, :] = altpoint_.histogram
            
            azpoint_ = mathMB.Histogram(subset.loc[sub_incl, 'azimuth'],
                                        bins=nazbins,
                                        range=[0, 2 * np.pi],
                                        weights=sub_weight[sub_incl])
            az_point[index, :] = azpoint_.histogram
        else:
            pass
        
    if smear_abundance:
        distribution['abundance_uncor'] = abundance.reshape(gridlongitude.shape)
    else:
        pass
    
    distribution['n_included'] = n_included.reshape(gridlongitude.shape)
    distribution['n_total'] = n_total.reshape(gridlongitude.shape)
    distribution['speed_dist_map'] = v_point.reshape(
        gridlongitude.shape + (nvelbins,))
    distribution['altitude_dist_map'] = alt_point.reshape(
        gridlongitude.shape + (naltbins,))
    distribution['azimuth_dist_map'] = az_point.reshape(
        gridlongitude.shape + (nazbins,))
    
    return distribution
