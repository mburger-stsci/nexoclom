import pickle
from scipy.io import readsav
import astropy.units as u
from nexoclom.utilities.exceptions import InputError


class SourceMap:
    def __init__(self, sourcemap):
        self.abundance = None
        self.longitude = None
        self.latitude = None
        self.speed = None
        self.speed_dist = None
        self.azimuth = None
        self.azimuth_dist = None
        self.altitude = None
        self.altitude_dist = None
        self.fraction_observed = None
        self.coordinate_system = 'solar-fixed'

        if isinstance(sourcemap, dict):
            self.load_dict(sourcemap)
        elif isinstance(sourcemap, str) and sourcemap.endswith('.pkl'):
            with open(sourcemap, 'rb') as file:
                sourcemap_ = pickle.load(file)
            if isinstance(sourcemap_, SourceMap):
                self.load_dict(sourcemap_.__dict__)
            elif isinstance(sourcemap_, dict):
                self.load_dict(sourcemap_)
            else:
                raise InputError('SourceMap', 'problem with mapfile')
        elif isinstance(sourcemap, str) and sourcemap.endswith('.sav'):
            sourcemap_ = readsav(sourcemap)
            self.abundance = sourcemap_.get('abundance', None)
            
            self.longitude = sourcemap_.get('longitude', None)
            if self.longitude is not None:
                self.longitude *= u.rad
                
            self.latitude = sourcemap_.get('latitude', None)
            if self.latitude is not None:
                self.latitude *= u.rad
                
            self.speed = sourcemap_.get('speed', None)
            if self.speed is not None:
                self.speed *= u.km/u.s

            self.speed_dist = sourcemap_.get('speed_dist', None)
            if self.speed_dist is not None:
                self.speed_dist *= u.km/u.s

            self.azimuth = sourcemap_.get('azimuth', None)
            if self.azimuth is not None:
                self.azimuth *= u.rad

            self.azimuth_dist = sourcemap_.get('azimuth_dist', None)
            if self.azimuth_dist is not None:
                self.azimuth_dist *= u.rad

            self.altitude = sourcemap_.get('altitude', None)
            if self.altitude is not None:
                self.altitude *= u.rad

            self.altitude_dist = sourcemap_.get('altitude_dist', None)
            if self.altitude_dist is not None:
                self.altitude_dist *= u.rad
                
            self.fraction_observed= sourcemap_.get('fraction_observed', None)
            self.coordinate_system = str(sourcemap_.get('coordinate_system',
                                                        'solar-fixed'))
        else:
            raise InputError('SourceMap', 'problem')
            
    def load_dict(self, sourcemap):
        self.abundance = sourcemap.get('abundance', None)
        self.longitude = sourcemap.get('longitude', None)
        self.latitude = sourcemap.get('latitude', None)
        self.speed = sourcemap.get('speed', None)
        self.speed_dist = sourcemap.get('speed_dist', None)
        self.azimuth = sourcemap.get('azimuth', None)
        self.azimuth_dist = sourcemap.get('azimuth_dist', None)
        self.altitude = sourcemap.get('altitude', None)
        self.altitude_dist = sourcemap.get('altitude_dist', None)
        self.fraction_observed = sourcemap.get('fraction_observed', None)
        self.coordinate_system = sourcemap.get('coordinate_system', 'solar-fixed')
        
    # def show_source_map:
    
    # def transform_reference_frame(self, output):
    #     """If the image center is not the planet, transform to a
    #        moon-centric reference frame."""
    #     assert 0, 'Not ready yet.'
    #
    #     # Load output
    #
    #     # # Transform to moon-centric frame if necessary
    #     # if result.origin != result.inputs.geometry.planet:
    #     #     assert 0, 'Need to do transparamsion for a moon.'
    #     # else:
    #     #     origin = np.array([0., 0., 0.])*output.x.unit
    #     #     sc = 1.
    #
    #     # Choose which packets to use
    #     # touse = output.frac >= 0 if keepall else output.frac > 0
    #
    #     # packet positions relative to origin -- not rotated
    #     # pts_sun = np.array((output.x[touse]-origin[0],
    #     #                     output.y[touse]-origin[1],
    #     #                     output.z[touse]-origin[2]))*output.x.unit
    #     #
    #     # # Velocities relative to sun
    #     # vels_sun = np.array((output.vx[touse],
    #     #                      output.vy[touse],
    #     #                      output.vz[touse]))*output.vx.unit
    #
    #     # Fractional content
    #     # frac = output.frac[touse]
    #
    #     return output #, pts_sun, vels_sun, frac
