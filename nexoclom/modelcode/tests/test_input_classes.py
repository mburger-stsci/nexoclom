""" Tests Input.__init__() and classes in input_classes.py
Note that this is more of a regression test than a unit test.
Compares with previously computed results."""
import os
import pickle
import pytest
import astropy.units as u
from astropy.time import Time
from nexoclom.solarsystem import SSObject
from nexoclom import Input, __file__ as basefile

basepath = os.path.dirname(basefile)
if __name__ == '__main__':
    inputpath = os.path.join('test_data', 'inputfiles')
else:
    inputpath = os.path.join(basepath, 'modelcode', 'tests', 'test_data', 'inputfiles')

@pytest.mark.modelcode
@pytest.mark.skip
def test_input_classes():
    saved_data = os.path.join(basepath, 'modelcode', 'tests', 'test_data',
                              'input_classes_data.pkl')
    inputfiles, results = pickle.load(open(saved_data, 'rb'))

    for inputfile_, result in zip(inputfiles, results):
        inputfile = os.path.join(os.path.dirname(__file__), inputfile_)
        inputs = Input(inputfile)
        if inputs != result:
            print(inputfile)
            print('geometry', inputs.geometry == result.geometry)
            print('surfaceinteraction', inputs.surfaceinteraction == result.surfaceinteraction)
            print('forces', inputs.forces == result.forces)
            print('spatialdist', inputs.spatialdist == result.spatialdist)
            print(inputs.spatialdist)
            print(result.spatialdist)
            print('speeddist', inputs.speeddist == result.speeddist)
            print('angulardist', inputs.angulardist == result.angulardist)
            print('options', inputs.options == result.options)

        assert Input(inputfile) == result

@pytest.mark.modelcode
def test_geometry():
    inputfile0 = os.path.join(inputpath, 'Geometry.01.input')
    geometry0 = Input(inputfile0).geometry
    result = {'planet': SSObject('Jupiter'),
              'startpoint': 'Io',
              'objects': {SSObject('Jupiter'), SSObject('Io'), SSObject('Europa')},
              'type': 'geometry without starttime',
              'phi': (1*u.rad, 2*u.rad),
              'subsolarpoint': (3.14*u.rad, 0*u.rad),
              'taa': 1.57*u.rad}
    assert geometry0.__dict__ == result
    
    inputfile1 = os.path.join(inputpath, 'Geometry.02.input')
    geometry1 = Input(inputfile1).geometry
    result = {'planet': SSObject('Jupiter'),
              'startpoint': 'Io',
              'objects': {SSObject('Jupiter'), SSObject('Io')},
              'type': 'geometry with starttime',
              'time': Time('2022-03-08T19:53:21')}
    assert geometry1.__dict__ == result

    inputfile2 = os.path.join(inputpath, 'Geometry.03.input')
    geometry2 = Input(inputfile2).geometry
    result = {'planet':SSObject('Mercury'),
              'startpoint': 'Mercury',
              'objects': {SSObject('Mercury')},
              'type': 'geometry without starttime',
              'subsolarpoint': (0 * u.rad, 0 * u.rad),
              'phi': None,
              'taa': 3.14 * u.rad}
    assert geometry2.__dict__ == result
    
    assert geometry0 == geometry0
    assert geometry0 != geometry1
    assert geometry0 != geometry2

@pytest.mark.modelcode
def test_SurfaceInteraction():
    # sticktype = 'constant'
    inputfile0 = os.path.join(inputpath, 'SurfaceInteraction.01.input')
    interaction0 = Input(inputfile0).surfaceinteraction
    result = {'sticktype': 'constant',
              'stickcoef': 1.,
              'accomfactor': None}
    assert interaction0.__dict__ == result

    inputfile1 = os.path.join(inputpath, 'SurfaceInteraction.02.input')
    interaction1 = Input(inputfile1).surfaceinteraction
    result = {'sticktype': 'constant',
              'stickcoef': 0.5,
              'accomfactor': 0.2}
    assert interaction1.__dict__ == result

    assert interaction0 == interaction0
    assert interaction0 != interaction1

    # sticktype = 'temperature dependent

if __name__ == '__main__':
    test_geometry()
    test_SurfaceInteraction()
