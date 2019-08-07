"""Classes used by the Inputs class"""
import os.path
import numpy as np
import pandas as pd
from astropy.time import Time
import astropy.units as u
from solarsystemMB import SSObject
from .database_connect import database_connect

class InputError(Exception):
    """Raised when a required parameter is not included in the inputfile."""
    def __init__(self, expression, message):
        self.expression = expression
        self.message = message


def isNone(x):
    try:
        q = x.value
    except:
        q = x

    if type(q) == str:
        return f'is NULL' if q is None else f"= '{q}'"
    else:
        return f'is NULL' if q is None else f"= {q}"


def inRange(field, x, delta):
    return f'ABS({field} - {x}) <= {delta/2}'

dtor = np.pi/180.


class Geometry:
    def __init__(self, gparam):
        """Define a Geometry object.
        
        See :doc:`inputfiles#Geometry` for more information.
        """
        # Planet
        if 'planet' in gparam:
            planet = gparam['planet'].title()
            self.planet = SSObject(planet)
        else:
            raise InputError('Geometry.__init__',
                             'Planet not defined in inputfile.')

        # All possible objects
        objlist = [self.planet.object]
        if self.planet.moons is not None:
            objlist.extend([m.object for m in self.planet.moons])
        else:
            pass

        # Choose the starting point
        self.startpoint = (gparam['startpoint'].title()
                           if 'startpoint' in gparam
                           else self.planet.object)
        if self.startpoint not in objlist:
            print(f'{self.startpoint} is not a valid starting point.')
            olist = '\n\t'.join(objlist)
            print(f'Valid choices are:\n\t{olist}')
            raise ValueError
        else:
            pass

        # Choose which objects to include
        # This is given as a list of names
        # Default = geometry.planet and geometry.startpoint
        if 'objects' in gparam:
            inc = set(i.strip().title()
                      for i in gparam['objects'].split(','))
        else:
            inc = {self.planet.object, self.startpoint}

        for i in inc:
            raise InputError('Geometry.__init__',
                             f'Invalid object {i} in geometry.include')
        self.objects = set(SSObject(o) for o in inc)

        # Different objects are created for geometry_with_starttime and
        # geometry_without_starttime
        if 'starttime' in gparam:
            self.type = 'geometry with starttime'
            try:
                self.time = Time(gparam['starttime'].upper())
            except:
                raise InputError('Geometry.__init__',
                         f'Invalid starttime format: {gparam["starttime"]}')
        else:
            self.type = 'geometry without starttime'
            if len(planet) == 1:
                self.phi = None
            elif 'phi' in gparam:
                phi_ = gparam['phi'].split(',')
                phi = tuple(float(p)*u.rad for p in phi_)
                nmoons = len(self.objects - {self.planet})
                if len(phi) == nmoons:
                    self.phi = phi
                else:
                    raise InputError('Geometry.__init__',
                          'The wrong number of orbital positions was given.')
            else:
                raise InputError('Geometry.__init__',
                                 'geometry.phi was not specified.')

            # Subsolar longitude and latitude
            if 'subsolarpoint' in gparam:
                subs = gparam['subsolarpoint'].split(',')
                try:
                    self.subsolarpoint = (float(subs[0])*u.rad,
                                          float(subs[1])*u.rad)
                except:
                    raise InputError('Geometry.__init__',
                         'The format for geometry.subsolarpoint is wrong.')
            else:
                self.subsolarpoint = (0*u.rad, 0*u.rad)

            # True Anomaly Angle
            self.taa = (float(gparam['taa'])*u.rad
                        if 'taa' in gparam
                        else 0.*u.rad)

    def __str__(self):
        result = ''
        for key,value in self.__dict__.items():
            result += (f'geometry.{key} = {value}\n')
        return result.strip()

    def search(self, startlist=None):
        # Make list of objects in planet system
        objs = [obj.object for obj in self.objects]
        objs.sort()
        objs2 = ','.join(objs)

        if startlist is None:
            startstr = ''
        else:
            start_ = [str(s) for s in startlist]
            startstr = f"and geo_idnum in ({', '.join(start_)})"

        if self.time is None:
            # Fields to query:
            #   planet, startpoint, objects, phi, subsolarpoint, TAA
            dtaa = (5.*u.deg).to(u.rad)
            taa = [self.taa-dtaa/2., self.taa+dtaa/2.]
            taa = [taa[0].value, taa[1].value]
            if taa[0] < 0.:
                taabit = '(taa>={} or taa<{})'.format(2*np.pi+taa[0], taa[1])
            elif taa[1] > 2*np.pi:
                taabit = '(taa>={} or taa<{})'.format(taa[0],
                                                      taa[1] % (2*np.pi))
            else:
                taabit = inRange('taa', self.taa.value, dtaa.value)

            phi = [p.value for p in self.phi]
            assert phi[0] == 0., 'phi for planet should be zero.'
            ptxt = [inRange('phi[{}]'.format(i+1), p, 5.*dtor) for
                    i,p in enumerate(phi)]
            ptxt2 = ' and '.join(ptxt)

            sspt0 = inRange('subsolarpt[0]', self.subsolarpoint[0].value,
                            5*dtor)
            sspt1 = inRange('subsolarpt[1]', self.subsolarpoint[1].value,
                            5*dtor)

            query = f'''SELECT geo_idnum FROM geometry
                        WHERE planet = '{self.planet.object}' and
                              startpoint = '{self.startpoint}' and
                              objects = ARRAY['{objs2}']::SSObject[] and
                              {ptxt2} and
                              {sspt0} and
                              {sspt1} and
                              {taabit} {startstr}'''
        else:
            # Fields to query
            # planet, StartPoint, objects, time
            # query =
            assert 0, 'Not working yet.'

        with database_connect() as con:
            result = pd.read_sql(query, con)
        if len(result) == 0:
            return None
        else:
            return result.geo_idnum.to_list()
###############################################################


class SurfaceInteraction:
    def __init__(self, sparam):
        """Define a SurfaceInteraction object.

        See :doc:`inputfiles#SurfaceInteraction` for more information.
        """
        sticktype = (sparam['sticktype'].lower()
                     if 'sticktype' in sparam
                     else None)
        if sticktype == 'temperature dependent':
            self.sticktype = sticktype
            
            if 'accomfactor' in sparam:
                self.accomfactor = sparam['accomfactor']
            else:
                raise InputError('SurfaceInteraction.__init__',
                                 'surface_interaction.accomfactor not given.')
        elif sticktype == 'surface map':
            self.sticktype = sticktype
            if 'stick_mapfile' in sparam:
                self.stick_mapfile = sparam['stick_mapfile']
            else:
                raise InputError('SurfaceInteraction.__init__',
                                 'surface_interaction.stick_mapfile not given.')
            
            if not os.path.exists(self.stick_mapfile):
                raise InputError('SurfaceInteraction.__init__',
                                 f'File Not Found: {self.stick_mapfile}')
            else:
                pass
            
            if 'accomfactor' in sparam:
                self.accomfactor = sparam['accomfactor']
            else:
                raise InputError('SurfaceInteraction.__init__',
                                 'surface_interaction.accomfactor not given.')
        elif 'stickcoef' in sparam:
            # Constant sticking
            self.sticktype = 'constant'
            self.stickcoef = sparam['stickcoef']
            if self.stickcoef < 0:
                self.stickcoef = 0
            elif self.stickcoef > 1:
                self.stickcoef = 1
            else:
                pass
        else:
            self.sticktype = 'constant'
            self.stickcoef = 1.
        
    def __str__(self):
        result = ''
        for key,value in self.__dict__.items():
            result += (f'surface_interaction.{key} = {value}\n')
        return result.strip()

    def search(self, startlist=None):
        if startlist is None:
            startstr = ''
        else:
            start_ = [str(s) for s in startlist]
            startstr = f"and st_idnum in ({', '.join(start_)})"

        if self.stickcoef == 1:
            query = f'''SELECT st_idnum FROM sticking_info
                        WHERE stickcoef=1 {startstr}'''
        else:
            query = f'''SELECT st_idnum FROM sticking_info
                        WHERE stickcoef={self.stickcoef}
                             tsurf {self.tsurf} and
                             stickfn {self.stickfn} and
                             stick_mapfile {self.stick_mapfile} and
                             epsilon {self.epsilon} and
                             n {self.n} and
                             tmin {self.tmin} and
                             emitfn {self.emitfn} and
                             accom_mapfile {self.accom_mapfile}
                             {startstr}'''

        with database_connect() as con:
            result = pd.read_sql(query, con)
        if len(result) == 0:
            return None
        else:
            return result.st_idnum.to_list()


class Forces:
    def __init__(self, fparam):
        """Define a Forces object.

        See :doc:`inputfiles#Forces` for more information.
        """
    
        self.gravity = (eval(fparam['gravity'])
                        if 'gravity' in fparam
                        else True)
        self.radpres = (eval(fparam['radpres'])
                        if 'radpres' in fparam
                        else True)

    def __str__(self):
        result = ''
        for key,value in self.__dict__.items():
            result += (f'forces.{key} = {value}\n')
        return result.strip()

    def search(self, startlist=None):
        if startlist is None:
            startstr = ''
        else:
            start_ = [str(s) for s in startlist]
            startstr = f"and f_idnum in ({', '.join(start_)})"

        query = f'''SELECT f_idnum FROM forces
                    WHERE gravity={self.gravity} and
                          radpres={self.radpres} {startstr}'''
        with database_connect() as con:
            result = pd.read_sql(query, con)

        if len(result) == 0:
            return None
        else:
            return result.f_idnum.to_list()

class SpatialDist:
    def __init__(self, sparam):
        """Define a SpatialDist object.

        See :doc:`inputfiles#SpatialDist` for more information.
        """
        if 'type' in sparam:
            self.type = sparam['type']
        else:
            raise InputError('SpatialDist.__init__',
                             'spatial_dist.type not given')
        
        if self.type == 'uniform':
            self.exobase = (float(sparam['exobase'])
                            if 'exobase' in sparam
                            else 1.)  # Unit gets set later
            if 'longitude' in sparam:
                lon0, lon1 = (float(l.strip())
                              for l in sparam['longitude'].split(','))
                lon0 = max(lon0, 0.)
                lon0 = min(lon0, 2*np.pi)
                lon1 = max(lon1, 0.)
                lon1 = min(lon1, 2*np.pi)
                self.longitude = (lon0*u.rad, lon1*u.rad)
            else:
                self.longitude = (0.*u.rad, 2*np.pi*u.rad)
                
            if 'latitude' in sparam:
                lat0, lat1 = (float(l.strip())*u.rad
                              for l in sparam['latitude'].split(','))
                lat0 = max(lat0, -np.pi/2)
                lat0 = min(lat0, np.pi/2)
                lat1 = max(lat1, -np.pi/2)
                lat1 = min(lat1, np.pi/2)
                if lat0 > lat1:
                    raise InputError('SpatialDist.__init__',
                         'spatial_dist.latitude[0] > spatial_dist.latitude[1]')
                self.latitude = (lat0*u.rad, lat1*u.rad)
            else:
                self.latitude = (-np.pi/2*u.rad, np.pi/2*u.rad)
        elif self.type == 'surface map':
            self.exobase = (float(sparam['exobase'])
                            if 'exobase' in sparam
                            else 1.)  # Unit gets set later
            
            if 'mapfile' in sparam:
                self.mapfile = sparam['mapfile']
            else:
                raise InputError('SpatialDist.__init__',
                                 'spatial_dist.mapfile not given.')
            
            if not os.path.exists(self.mapfile):
                raise InputError('SpatialDist.__init__',
                                 f'File Not Found: {self.mapfile}')
            else:
                pass
            
            if 'coordsystem' in sparam:
                coordsystems = {'solar-fixed', 'planet-fixed', 'moon-fixed'}
                if sparam['coordsystem'] in coordsystems:
                    self.coordsystem = sparam['coordsystem']
                else:
                    raise InputError('SpatialDist.__init__',
                     f'{sparam["coordsystem"]} not a valid coordinat system.')
            else:
                self.coordsystem = 'default'
        elif self.type == 'surface spot':
            self.exobase = (float(sparam['exobase'])
                            if 'exobase' in sparam
                            else 1.)  # Unit gets set later
            if 'longitude' in sparam:
                self.longitude = float(sparam['longitude'])*u.rad
            else:
                raise InputError('SpatialDist.__init__',
                                 'spatial_dist.longitude not given.')
            
            if 'latitude' in sparam:
                self.latitude = float(sparam['latitude'])*u.rad
            else:
                raise InputError('SpatialDist.__init__',
                                 'spatial_dist.latitude not given.')

            if 'sigma' in sparam:
                self.sigma = float(sparam['sigma'])*u.rad
            else:
                raise InputError('SpatialDist.__init__',
                                 'spatial_dist.sigma not given.')
        else:
            raise InputError('SpatialDist.__init__',
                             f'spatial_dist.type = {self.type} not defined.')

    def __str__(self):
        result = ''
        for key,value in self.__dict__.items():
            result += (f'spatial_dist.{key} = {value}\n')
        return result.strip()

    def search(self, startlist=None):
        if startlist is None:
            startstr = ''
        else:
            start_ = [str(s) for s in startlist]
            startstr = f"and spat_idnum in ({', '.join(start_)})"

        if self.longitude is None:
            long0 = 'longitude[1] = 0.'
            long1 = 'longitude[2] = 0.'
        else:
            long0 = inRange('longitude[1]', self.longitude[0].value, 5*dtor)
            long1 = inRange('longitude[2]', self.longitude[1].value, 5*dtor)

        if self.latitude is None:
            lat0 = 'latitude[1] = 0.'
            lat1 = 'latitude[2] = 0.'
        else:
            lat0 = inRange('latitude[1]', self.latitude[0].value, 5*dtor)
            lat1 = inRange('latitude[2]', self.latitude[1].value, 5*dtor)

        query = f'''SELECT spat_idnum FROM spatialdist
                    WHERE type = '{self.type}' and
                         {inRange('exobase', self.exobase, 0.05)} and
                         use_map {isNone(self.use_map)} and
                         mapfile {isNone(self.mapfile)} and
                         {long0} and
                         {long1} and
                         {lat0} and
                         {lat1} {startstr}'''

        with database_connect() as con:
            result = pd.read_sql(query, con)
        if len(result) == 0:
            return None
        else:
            return result.spat_idnum.to_list()
###############################################################


class SpeedDist:
    """Define a SpeedDist object.

    See :doc:`inputfiles#SpeedDist` for more information.
    """
    
    def __init__(self, sparam):
        self.type = sparam['type']

        if self.type == 'gaussian':
            if 'vprob' in sparam:
                self.vprob = float(sparam['vprob'])*u.km/u.s
            else:
                raise InputError('SpatialDist.__init__',
                                 'speed_dist.vprob not given.')
            if 'sigma' in sparam:
                self.sigma = float(sparam['sigma'])*u.km/u.s
            else:
                raise InputError('SpatialDist.__init__',
                                 'speed_dist.sigma not given.')
        elif self.type == 'sputtering':
            if 'alpha' in sparam:
                self.alpha = float(sparam['alpha'])
            else:
                raise InputError('SpatialDist.__init__',
                                 'speed_dist.alpha not given.')
            if 'beta' in sparam:
                self.beta = float(sparam['beta'])
            else:
                raise InputError('SpatialDist.__init__',
                                 'speed_dist.beta not given.')
            if 'u' in sparam:
                self.U = float(sparam['u'])*u.eV
            else:
                raise InputError('SpatialDist.__init__',
                                 'speed_dist.U not given.')
        elif self.type == 'maxwellian':
            if 'temperature' in sparam:
                self.alpha = float(sparam['temperature'])*u.K
            else:
                raise InputError('SpatialDist.__init__',
                                 'speed_dist.temperature not given.')
        elif self.type == 'flat':
            if 'vprob' in sparam:
                self.alpha = float(sparam['vprob'])*u.km/u.s
            else:
                raise InputError('SpatialDist.__init__',
                                 'speed_dist.vprob not given.')
            
            if 'delv' in sparam:
                self.alpha = float(sparam['delv'])*u.km/u.s
            else:
                raise InputError('SpatialDist.__init__',
                                 'speed_dist.delv not given.')

        else:
            assert 0, f'SpeedDist.type = {self.type} not available'

    def __str__(self):
        result = ''
        for key,value in self.__dict__.items():
            result += (f'speed_dist{key} = {value}\n')
        return result.strip()

    def search(self, startlist=None):
        if startlist is None:
            startstr = ''
        else:
            start_ = [str(s) for s in startlist]
            startstr = f"and spd_idnum in ({', '.join(start_)})"

        if self.vprob is None:
            vstr = 'vprob is NULL'
        else:
            vstr = inRange('vprob', self.vprob.value,
                           self.vprob.value*0.05)

        if self.temperature is None:
            Tstr = 'temperature is NULL'
        else:
            Tstr = inRange('temperature', self.temperature.value,
                           self.temperature.value*0.05)

        query = f'''SELECT spd_idnum FROM speeddist
                    WHERE type = '{self.type}' and
                          {vstr} and
                          sigma {isNone(self.sigma)} and
                          U  {isNone(self.U)} and
                          alpha  {isNone(self.alpha)} and
                          beta  {isNone(self.beta)} and
                          {Tstr} and
                          delv  {isNone(self.delv)} {startstr}'''

        with database_connect() as con:
            result = pd.read_sql(query, con)
        if len(result) == 0:
            return None
        else:
            return result.spd_idnum.to_list()


class AngularDist:
    def __init__(self, aparam):
        """Define a AngularDist object.

        See :doc:`inputfiles#AngularDist` for more information.
        """
        if 'type' in aparam:
            self.type = aparam['type']
            if self.type == 'radial':
                pass
            elif self.type == 'isotropic':
                if 'azimuth' in aparam:
                    az0, az1 = (float(l.strip())
                                for l in aparam['azimuth'].split(','))
                    az0 = max(az0, 0.)
                    az0 = min(az0, 2*np.pi)
                    az1 = max(az1, 0.)
                    az1 = min(az1, 2*np.pi)
                    self.aziumth = (az0*u.rad, az1*u.rad)
                else:
                    self.azimuth = (0*u.rad, 2*np.pi*u.rad)

                if 'altitude' in aparam:
                    alt0, alt1 = (float(l.strip())*u.rad
                                  for l in aparam['altitude'].split(','))
                    alt0 = max(alt0, -np.pi/2)
                    alt0 = min(alt0, np.pi/2)
                    alt1 = max(alt1, -np.pi/2)
                    alt1 = min(alt1, np.pi/2)
                    if alt0 > alt1:
                        raise InputError('AngularDist.__init__',
                         'angular_dist.altitude[0] > angular_dist.altitude[1]')
                    self.altitude = (alt0*u.rad, alt1*u.rad)
                else:
                    self.altitude = (-np.pi/2*u.rad, np.pi/2*u.rad)
            else:
                raise InputError('AngularDist.__init__',
                             f'angular_dist.type = {self.type} not defined.')
        else:
            self.type = 'isotropic'
            self.azimuth = (0*u.rad, 2*np.pi*u.rad)
            self.altitude = (-np.pi/2*u.rad, np.pi/2*u.rad)

    def __str__(self):
        result = ''
        for key,value in self.__dict__.items():
            result += (f'angular_dist.{key} = {value}\n')
        return result.strip()

    def search(self, startlist=None):
        if startlist is None:
            startstr = ''
        else:
            start_ = [str(s) for s in startlist]
            startstr = f"and ang_idnum in ({', '.join(start_)})"

        if self.azimuth is None:
            az0 = 'azimuth[1] is NULL'
            az1 = 'azimuth[2] is NULL'
        else:
            az0 = inRange('azimuth[1]', self.azimuth[0].value, 5*dtor)
            az1 = inRange('azimuth[2]', self.azimuth[1].value, 5*dtor)

        if self.altitude is None:
            alt0 = 'altitude[1] is NULL'
            alt1 = 'altitude[2] is NULL'
        else:
            alt0 = inRange('altitude[1]', self.altitude[0].value, 5*dtor)
            alt1 = inRange('altitude[2]', self.altitude[1].value, 5*dtor)
        n = isNone(self.n)

        query = f'''SELECT ang_idnum from angulardist
                    WHERE type = '{self.type}' and
                          {az0} and {az1} and
                          {alt0} and {alt1} and
                          n {n} {startstr}'''
        with database_connect() as con:
            result = pd.read_sql(query, con)

        if len(result) == 0:
            return None
        else:
            return result.ang_idnum.to_list()
###############################################################


class Options:
    def __init__(self, oparam):
        """Define a Options object.

        See :doc:`inputfiles#Options` for more information.
        """
        if 'endtime' in oparam:
            self.endtime = float(oparam['endtime'])*u.s
        else:
            raise InputError('Options.__init__',
                             'options.endtime not specified.')

        if 'species' in oparam:
            self.species = oparam['species']
        else:
            raise InputError('Options.__init__',
                             'options.species not specified.')

        self.lifetime = (float(oparam['lifetime'])*u.s
                         if 'lifetime' in oparam
                         else 0.*u.s)
        
        self.outeredge = (float(oparam['outer_edge'])
                          if 'outeredge' in oparam
                          else 1e30)
                          
        self.step_size = (float(oparam['step_size'])
                          if 'step_size' in oparam
                          else 0.)

        if self.step_size == 0:
            self.resolution = (float(oparam['resolution'])
                               if 'resolution' in oparam else 1e-4)
        else:
            self.resolution = None

    def __str__(self):
        result = ''
        for key,value in self.__dict__.items():
            result += (f'options.{key} = {value}\n')
        return result.strip()

    def search(self, startlist=None):
        if startlist is None:
            startstr = ''
        else:
            start_ = [str(s) for s in startlist]
            startstr = f"and opt_idnum in ({', '.join(start_)})"

        endtime = inRange('endtime', self.endtime.value,
                          self.endtime.value*0.05)
        outeredge = isNone(self.outeredge)
        nsteps = isNone(self.nsteps)
        res = isNone(self.resolution)

        query = f'''SELECT opt_idnum from options
                    WHERE {endtime} and
                          resolution {res} and
                          at_once = {self.at_once} and
                          atom = '{self.atom}' and
                          lifetime = {self.lifetime.value} and
                          fullsystem = {self.fullsystem} and
                          outeredge {outeredge} and
                          motion = {self.motion} and
                          streamlines = {self.streamlines} and
                          nsteps {nsteps} {startstr}'''
        with database_connect() as con:
            result = pd.read_sql(query, con)

        if len(result) == 0:
            return None
        else:
            return result.opt_idnum.to_list()
