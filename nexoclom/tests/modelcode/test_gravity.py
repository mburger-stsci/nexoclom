""" Test that atoms ejected radially outward return to the surface at the
proper rate and .
"""
import os
import numpy as np
import pytest
import astropy.units  as u
from nexoclom import Input, Output, __file__ as basefile
from nexoclom.solarsystem import SSObject

basepath = os.path.dirname(basefile)
if __name__ == '__main__':
    inputpath = '/Users/mburger/Work/Research/NeutralCloudModel/nexoclom/nexoclom/'
    inputpath = os.path.join(inputpath, 'tests', 'test_data', 'inputfiles')
else:
    inputpath = os.path.join(basepath, 'tests', 'test_data', 'inputfiles')


@pytest.mark.modelcode
def test_gravity_with_constant_stepsize():
    Mercury = SSObject('Mercury')
    
    inputfile = os.path.join(inputpath, 'Gravity.input')
    inputs = Input(inputfile)
    inputs.run(1e4, 1e4, overwrite=False, compress=False)
    _, outputfiles, _, _ = inputs.search()
    output = Output.restore(outputfiles[0])
    v0 = np.sqrt(output.X0.vx**2 + output.X0.vy**2 + output.X0.vz**2)
    v0 = (v0.values * output.unit/u.s).to(u.km/u.s)
    output.X['r'] = np.sqrt(output.X.x**2 + output.X.y**2 + output.X.z**2)
    
    time = inputs.options.endtime - output.X.time.values * u.s
    rmax = np.zeros(len(output.X0))*output.unit
    bound = np.array([False for _ in range(len(rmax))])
    t_end = np.zeros(len(output.X0))
    
    # Find the last point for each particle
    for ind in output.X.Index.unique():
        sub = output.X[output.X.Index == ind]
        t0 = np.where(sub.frac == 0)[0]
        rmax[ind] = sub.r.max()*output.unit
        
        # Chck to see if things start to come back down
        bound[ind] = sub.iloc[-1, -1] < sub.iloc[-2, -1]
        if len(t0) > 0:
            t_end[ind] = sub.iloc[t0[0]-1].time + inputs.options.step_size
        else:
            pass
    
    # Verify those with < escape velocity hit surface
    # v_escape = np.sqrt(-2*Mercury.GM/Mercury.radius).to(u.km/u.s)
    peak_altitude = (Mercury.GM/(.5*v0**2 + Mercury.GM/Mercury.radius)).to(output.unit)
    
    assert np.all(np.isclose(rmax[bound], peak_altitude[bound]))

    # Verify position is correct
    
if __name__ == '__main__':
    test_gravity_with_constant_stepsize()
