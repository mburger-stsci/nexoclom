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
from MESSENGERuvvs import MESSENGERdata

inputfiles = [#'Na.Gaussian.3_1.no_accom_tempsticking.input']
              #'Na.Gaussian.3_1.no_accom.input',
              # 'Na.maxwellian.1200.accom.input',
              #'Ca.surfacemap.maxwellian.50000.input',
              'Ca.spot.maxwellian.input']

orbit = 36
overwrite = False

def test_Input():
    for infile in inputfiles:
        data = MESSENGERdata(infile[0:2], f'orbit={orbit}')
        inputs = Input(os.path.join(os.path.dirname(__file__),
                                    'inputfiles', infile))
        # inputs.run(1e5, overwrite=overwrite)
        # image = inputs.produce_image('inputfiles/MercuryEmission.format',
        #                              overwrite=overwrite)
        # sfile = os.path.join(os.path.dirname(__file__),
        #                      'outputs', infile.replace('.input', '.html'))
        # image.display(show=True, savefile=sfile)

        data.model(inputs, 1e5, overwrite=overwrite, fit_method='middle50',
                   label='Middle 50')
        data.model(inputs, 1e5, overwrite=overwrite, fit_method='middle10',
                   label='Middle 10')
        data.model(inputs, 1e5, overwrite=overwrite, fit_method='middle90',
                   label='Middle 90')

        sfile = os.path.join(os.path.dirname(__file__), 'figures',
                             f'{infile}_Orbit{orbit:04d}.html')
        data.plot(filename=sfile)


    assert True
    return data

if __name__ == '__main__':
    data = test_Input()
