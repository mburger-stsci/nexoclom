"""Test that all types of inputfiles can be read successfully.

To Test:
    (g.1) Geometry With Time Stamp
    (g.2) Geometry Without Time Stamp
    (g.3) Startpoint given
    (g.4) Startpoint not given
    (g.5) Specifiy objects for planet with moon
    (g.6) Don't specify objects for planet with moon
    (g.7) phi for planets with moon
    
    (si.1) Complete sticking
    (si.2) Constant sticking with some bouncing
    (si.3) Temperature dependent sticking
    (si.4) Sticking coefficient from a surface map
    
    (f.1) Gravity on
    (f.2) Gravity off
    (f.3) Radiation pressure on
    (f.4) Radiation pressure off
    
"""
import os.path
try:
    from ..Input import Input
except:
    from nexoclom import Input

def test_Input():
    inputs = Input(os.path.join(os.path.dirname(__file__),
                                'inputfiles',
                                'Ca.isotropic.maxwellian.50000.input'))
    inputs.run(1e5)

    print(inputs)
    assert True
    
    
if __name__ == '__main__':
    test_Input()
