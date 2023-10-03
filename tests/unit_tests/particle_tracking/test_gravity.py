""" Test that gravity is working correctly.
For particles ejected from the surface subject only to gravity
    (1) Conserve energy
        E/m = KE/m + PE/m = v**2/2 + GM/r
    (2) Escape when they are supposed to or return to the surface
    (3) Bound atoms reach proper altitude
    (4) Bound particles hit surface with initial velocity
"""
import os
import numpy as np
import pytest
import astropy.units  as u
from nexoclom import Input, Output, __path__
from nexoclom.solarsystem import SSObject


basepath = os.path.dirname(__path__[0])
inputpath = os.path.join(basepath, 'tests', 'test_data', 'inputfiles')


@pytest.mark.particle_tracking
def test_gravity():
    """
    Test cases to run
        (1) constant step size
        (2) variable stepsize
    """
    Mercury = SSObject('Mercury')
    
    outputs = []
    for run in ('constant', 'variable'):
        # set up the inputs
        inputfile = os.path.join(inputpath, 'Gravity.input')
        inputs = Input(inputfile)
        if run == 'variable':
            inputs.options.step_size = 0
            inputs.options.resolution = 0.0001
        else:
            inputs.options.step_size = 30
        
        # Run the model
        inputs.run(1e3, 1e3, overwrite=False, compress=False)
        _, outputfiles, _, _ = inputs.search()
        output = Output.restore(outputfiles[0])
    
        # Check for energy conservation
        radius = np.sqrt(output.X.x.values**2 + output.X.y.values**2 +
                         output.X.z.values**2) * output.unit
        velocity = np.sqrt(output.X.vx.values**2 + output.X.vy.values**2 +
                           output.X.vz.values**2) * output.unit/u.s
        energy = 0.5*velocity**2 + Mercury.GM/radius
        
        for i in output.X.Index.unique():
            e = energy[(output.X.Index == i) & np.isfinite(energy)]
            assert np.all(np.isclose(e, e.mean()))
            
        outputs.append(output)
        
    from inspect import currentframe, getframeinfo
    frameinfo = getframeinfo(currentframe())
    print(frameinfo.filename, frameinfo.lineno)
    from IPython import embed; embed()
    import sys; sys.exit()
    
    # v0 = np.sqrt(output.X0.vx**2 + output.X0.vy**2 + output.X0.vz**2)
    # v0 = (v0.values * output.unit/u.s).to(u.km/u.s)
    # output.X['r'] = np.sqrt(output.X.x**2 + output.X.y**2 + output.X.z**2)
    #
    # time = inputs.options.endtime - output.X.time.values * u.s
    # rmax = np.zeros(len(output.X0))*output.unit
    # bound = np.array([False for _ in range(len(rmax))])
    # t_end = np.zeros(len(output.X0))
    #
    # # Find the last point for each particle
    # for ind in output.X.Index.unique():
    #     sub = output.X[output.X.Index == ind]
    #     t0 = np.where(sub.frac == 0)[0]
    #     rmax[ind] = sub.r.max()*output.unit
    #
    #     # Chck to see if things start to come back down
    #     bound[ind] = sub.iloc[-1, -1] < sub.iloc[-2, -1]
    #     if len(t0) > 0:
    #         t_end[ind] = sub.iloc[t0[0]-1].time + inputs.options.step_size
    #     else:
    #         pass
    #
    # # Verify those with < escape velocity hit surface
    # # v_escape = np.sqrt(-2*Mercury.GM/Mercury.radius).to(u.km/u.s)
    # peak_altitude = (Mercury.GM/(.5*v0**2 + Mercury.GM/Mercury.radius)).to(output.unit)
    #
    # assert np.all(np.isclose(rmax[bound], peak_altitude[bound]))
    #
    # # Verify position is correct
    #
if __name__ == '__main__':
    test_gravity()
