"""Tests for source_distribution.py
Functions to test:
    gaussiandist (Not tested - this isn't used)
    sputdist
    MaxwellianDist
    xyz_from_lonlat
    surface_distribution
    surface_spot
    idlversion
    speed_distribution
    angular_distribution

Distributions to be compared with results from IDL version of code"""
import numpy as np
import os.path
from scipy.io import readsav
import astropy.units as u
import bokeh.plotting as plt
try:
    from ..source_distribution import *
except ValueError:
    from nexoclom.source_distribution import *

TOL = 1e-5


def test_maxwellian_dist():
    idl = readsav(os.path.join(os.path.dirname(__file__), 'data',
                               'maxdist.sav'))
    vv =  idl['vv']*u.km/u.s
    temperature = [t*u.K for t in idl['temp']]
    atom = [str(a)[2:-1] for a in idl['atom']]
    f_idl = idl['max_dist']

    for i,t in enumerate(temperature):
        for j,a in enumerate(atom):
            f_idl_ = np.array(f_idl[j,i,:], dtype=float)
            f_idl_ /= np.max(f_idl_)
            f_py = MaxwellianDist(vv, t, a)

            print(np.max(np.abs(f_idl_-f_py)))
            assert np.all(np.abs(f_idl_-f_py) < TOL), (
                f'atom={a}, temperature={t}')

def test_sputdist():
    idl = readsav(os.path.join(os.path.dirname(__file__), 'data',
                               'sputdist.sav'))
    vv =  idl['vv']*u.km/u.s
    atom = [str(a)[2:-1] for a in idl['atom']]
    alpha = [a for a in idl['alpha']]
    bet = [b for b in idl['bet']]
    U = [u_*u.eV for u_ in idl['U']]
    f_idl = idl['sput_dist']

    for i,al in enumerate(alpha):
        for j,b in enumerate(bet):
            for k,u_ in enumerate(U):
                for l,a in enumerate(atom):
                    f_idl_ = np.array(f_idl[l,k,j,i,:], dtype=float)
                    f_idl_ /= np.max(f_idl_)
                    f_py = sputdist(vv, u_, al, b, a)

            print(np.max(np.abs(f_idl_-f_py)))
            assert np.all(np.abs(f_idl_-f_py) < TOL), (
                f'atom={a}, alpha={al}, beta={b}, U={u_}')


if __name__ == '__main__':
    test_maxwellian_dist()
    test_sputdist()
