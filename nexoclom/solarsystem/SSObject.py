"""Create an object for each Solar System body containing basic information.
Information stored:
* object: Object name
* orbits: Body that the object orbits
* radius: in km
* mass: in kg
* a: orbital semi major axis. In AU if orbits the Sun; km
    if orbits a planet
* e: orbital eccentricity
* tilt: tilt of planetary axis in degrees
* rotperiod: rotation period in hours
* orbperiod: orbital period in days
* GM: mass * G in m**3/s**2
* moons: returned as a list of SSObjects

Values are astropy units quantities when appropriate.
"""
import os
import pandas as pd
from astropy import constants as c
from astropy import units as u
import spiceypy as spice
from astroquery.jplhorizons import Horizons
from nexoclom import __file__ as basefile
from nexoclom.solarsystem.spice_routines import load_kernels


basepath = os.path.dirname(basefile)
pklfile = os.path.join(basepath, 'data', 'PlanetaryConstants.pkl')

class SSObject:
    """Creates Solar System object."""
    def __init__(self, obj):
        if not os.path.exists(pklfile):
            set_up_planetary_constants()
        else:
            pass

        constants = pd.read_pickle(pklfile)
        
        if obj.title() in constants.index:
            row = constants.loc[obj.title()]

            self.object = obj.title()
            self.orbits = row.orbits
            self.radius = row.radius * u.km
            self.mass = row.mass * u.kg
            self.a = row.a
            self.e = row.e
            self.naif_id = row.naif_id
            # self.tilt = row.tilt * u.deg
            # self.rotperiod = row.rot_period * u.h
            # self.orbperiod = row.orb_period * u.d
            self.GM = row.GM * u.m**3/u.s**2

            self.moons = [SSObject(moon) for moon in
                constants.loc[constants.orbits == self.object].index]
            if len(self.moons) == 0:
                self.moons = None
            else:
                pass

            if self.orbits == 'Milky Way':
                self.type = 'Star'
                self.a *= u.km
            elif self.orbits == 'Sun':
                self.type = 'Planet'
                self.a *= u.au
            else:
                self.type = 'Moon'
                self.a *= u.km
        else:
            print(f'Object {obj} does not exist in table.')
            self.object = None

    def __len__(self):
        # Returns number of objects (e.g. Planet + moons) in the SSObeject
        return 1 if self.moons is None else len(self.moons)+1

    def __eq__(self, other):
        return self.object == other.object

    def __hash__(self):
        return hash((self.object, ))
    
    def __repr__(self):
        return self.__str__()

    def __str__(self):
        out = (f'Object: {self.object}\n'
               f'NAIF ID: {self.naif_id}\n'
               f'Type = {self.type}\n'
               f'Orbits {self.orbits}\n'
               f'Equatorial Radius = {self.radius:0.2f}\n'
               f'Mass = {self.mass:0.2e}\n'
               f'a = {self.a:0.2f}\n'
               f'Eccentricity = {self.e:0.2f}\n'
               f'GM = {self.GM:0.2e}')
               # f'Tilt = {self.tilt:0.2f}\n'
               # f'Rotation Period = {self.rotperiod:0.2f}\n'
               # f'Orbital Period = {self.orbperiod:0.2f}\n'
        return out
    
def set_up_planetary_constants():
    load_kernels()
    objfile = os.path.join(basepath, 'data', 'SolarSystemContents.csv')
    objects = pd.read_csv(objfile).set_index('object')
    
    constants = pd.DataFrame(columns=['orbits', 'radius', 'mass', 'a', 'e',
                                      'GM', 'naif_id'],
                             index=objects.index)
    constants.loc[constants.index, 'orbits'] = objects.loc[constants.index, 'orbits']
    
    naif_file = os.path.join(basepath, 'data', 'naif_ids.dat')
    naif_ids = pd.read_csv(naif_file, sep=':')
    naif_ids.columns = [x.strip() for x in naif_ids.columns]
    naif_ids.NAME = naif_ids.NAME.apply(
        lambda name: name.strip().split(' ')[0].strip().replace("'", ""))
    naif_ids.drop_duplicates('NAME', inplace=True)
    naif_ids.set_index('NAME', inplace=True)
    for obj in constants.index:
        print(obj)
        if obj.upper() in naif_ids.index:
            constants.loc[obj, 'naif_id'] = int(naif_ids.loc[obj.upper(),
                                                             'NAIF ID'])
        else:
            pass
        
        radius = spice.bodvrd(obj.upper(), 'RADII', 3)
        constants.loc[obj, 'radius'] = radius[1][:1].mean()
        gm = spice.bodvrd(obj.upper(), 'GM', 1)[1][0] * u.km**3/u.s**2
        constants.loc[obj, 'GM'] = -1*gm.to(u.m**3/u.s**2).value
        constants.loc[obj, 'mass'] = (gm/c.G).to(u.kg).value
        
        if obj == 'Sun':
            constants.loc[obj, ['a', 'e']] = 0
        else:
            hobj = Horizons(id=constants.loc[obj, 'naif_id'],
                            location=0,
                            epochs={'start': '2023-01-01',
                                    'stop': '2024-01-01',
                                    'step': '1d'})
            elements = hobj.elements()
            constants.loc[obj, 'a'] = elements['a'].mean()
            constants.loc[obj, 'e'] = elements['e'].mean()
            # constants.loc[obj, 'orbperiod'] = elements['P'].mean()
            # constants.loc[obj, 'rotperiod'] = 0.
            # constants.loc[obj, 'tilt'] = 0
        
    spice.kclear()
    constants.to_pickle(pklfile)
    constants.to_csv(pklfile.replace('.pkl', '.csv'), index=False)
    
if __name__ == '__main__':
    set_up_planetary_constants()
