"""Computes acceleration and ionization on a packet due to specified forces

Gravitational acceleration

Equations of motion:
    dvxdt = sum_objects (GM * (x-x_obj))/(r_obj)^3
    dvydt = sum_objects (GM * (y-y_obj))/(r_obj)^3
    dvzdt = sum_objects (GM * (z-z_obj))/(r_obj)^3
        -- r_obj = sqrt( (x-x_obj)^2 + (y-y_obj)^2 + (z-z_obj)^2 )
    dndt = instantaneous change in density

Current version: Assumes there is only a planet -- does not do moons yet
"""
import numpy as np


def State(t, x, v, X, output):
    xx = ['x', 'y', 'z']
    rh = ['x', 'z']
    # compute gravitational acceleration
    if output.inputs.forces.gravity:
        r3 = (x[0,:]**2 + x[1,:]**2 + x[2,:]**2)**1.5
        agrav = output.GM * x/r3
        
        R3 = (np.sum(X[xx].values**2, axis=1))**1.5
        Agrav = output.GM * X[xx].values/R3[:,np.newaxis]
    else:
        agrav = np.zeros_like(x)
        Agrav = np.zeros(X[xx].shape)
        
    assert np.all(r3 == R3)
    assert np.all(agrav.transpose() == Agrav)

    # compute radiation acceleration
    arad = np.zeros_like(x)
    Arad = np.zeros(X[xx].shape)
    if output.inputs.forces.radpres:
        rho = x[0,:]**2 + x[2,:]**2
        Rho = np.sum(X[rh]**2, axis=1)
        assert np.all(rho == Rho)
        
        out_of_shadow = (rho > 1) | (x[1,:] < 0)
        Out_of_shadow = (Rho > 1) | (X['y'] < 0)
        assert np.all(out_of_shadow == Out_of_shadow)

        # radial velocity of each packet realtive to the Sun
        vv = v[1,:] + output.vrplanet
        VV = X['vy'] + output.vrplanet

        # Compute radiation acceleration
        arad[1,:] = (np.interp(vv, output.radpres.velocity,
                               output.radpres.accel) * out_of_shadow)
        Arad[:,1] = np.interp(VV, output.radpres.velocity,
                              output.radpres.accel) * Out_of_shadow
    else:
        pass
    
    assert np.all(Arad == arad.transpose())

    # Compute total acceleration
    accel = agrav + arad
    Accel = Agrav + Arad
    assert np.all(accel.transpose() == Accel)
    assert np.all(np.isfinite(accel))

    # Compute ionization rate
    if output.inputs.options.lifetime > 0:
        # Explicitly set lifetime
        ionizerate = np.ones_like(t)/output.inputs.options.lifetime.value
        Ionizerate = np.ones(X['t'].shape)/output.inputs.options.lifetime.value
    else:
        if output.loss_info.photo is not None:
            # Compute photoionization rate
            rho = x[0,:]**2 + x[2,:]**2
            out_of_shadow = (rho > 1) | (x[1,:] < 0)
            photorate = output.loss_info.photo * out_of_shadow
            
            Rho = np.sum(X[rh]**2, axis=1)
            Out_of_shadow = (Rho > 1) | (X['y'] < 0)
            Photorate = output.loss_info.photo * out_of_shadow
        else:
            photorate = np.zeros(X['t'].shape)
         
        assert np.all(photorate == Photorate)

        '''
        magcoord = xyz_to_magcoord(t, x, output.inputs, output.planet)

        if output.loss_info.eimp:
            # Compute electron impact rate
            assert 0, 'Electron impacts not set up yet'
        else:
            eimprate = 0.

        if output.loss_info.chX:
            # Compute charge exchange rate
            assert 0, 'Charge exchange not set up yet'
        else:
            chxrate = 0.
        '''

        ionizerate = photorate  #+ eimprate + chxrate

    return accel, ionizerate
