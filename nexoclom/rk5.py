import numpy as np
from astropy.time import Time
from .state import state
from .State_ import State

# RK coefficients
c = np.array([0, 0.2, 0.3, 0.8, 8./9., 1., 1.])
b = np.array([35./384., 0., 500./1113., 125./192., -2187./6784.,
              11./84., 0.])
bs = np.array([5179./57600., 0., 7571./16695., 393./640., -92097./339200.,
      187./2100., 1./40.])
bd = b - bs

a = np.zeros((7,7))
a[1,0] = 0.2
a[2,:2] = [3./40., 9./40.]
a[3,:3] = [44./45., -56./15., 32./9.]
a[4,:4] = [19372./6561., -25360./2187., 64448./6561., -212./729.]
a[5,:5] = [9017./3168., -355./33., 46732./5247., 49./176., -5103./18656.]
a[6,:] =  b

c2, c3, c4, c5, c6, c7 = 0.2, 0.3, 0.8, 8./9., 1., 1.
b1, b2, b3, b4, b5, b6, b7 = (35./384., 0., 500./1113., 125./192.,
                              -2187./6784., 11./84., 0.)
b1s, b2s, b3s, b4s, b5s, b6s, b7s = (5179./57600., 0., 7571./16695.,
    393./640., -92097./339200., 187./2100., 1./40.)
b1d, b2d, b3d, b4d, b5d, b6d, b7d = (b1-b1s, b2-b2s, b3-b3s, b4-b4s,
                                     b5-b5s, b6-b6s, b7-b7s)

a21 = 0.2
a31, a32 = 3./40., 9./40.
a41, a42, a43 = 44./45., -56./15., 32./9.
a51, a52, a53, a54 = (19372./6561., -25360./2187., 64448./6561.,
                      -212./729.)
a61, a62, a63, a64, a65 = (9017./3168., -355./33., 46732./5247.,
                           49./176., -5103./18656.)
a71, a72, a73, a74, a75, a76 = b1, b2, b3, b4, b5, b6


def rk5(output, X0):
    ''' Perform a single rk5 step '''
    cols = ['time', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'frac']

    t1 = X0.time.values
    x1 = X0[['x', 'y', 'z']].values.transpose()
    v1 = X0[['vx', 'vy', 'vz']].values.transpose()
    f1 = X0.frac.values
    h = X0.h.values
    
    f1 = np.log(f1)
    X0.frac = X0.frac.apply(np.log)
    assert np.all(f1 == X0.frac)
    
    t0_ = Time.now()
    x = np.zeros((len(X0), len(cols), 8))
    x[:,:,0] = X0[cols].values
    accel = np.zeros((len(X0), 3, 8))
    ioniz = np.zeros((len(X0), 8))

    for n in range(6):
        accel[:,:,n], ioniz[:,n] = state(x[:,:,n], output)
        x[:,0,n+1] = -h*c[n+1]
        for i in np.arange(n+1):
            x[:,1:4,n+1] += h[:,np.newaxis]*a[n+1,i]*x[:,4:7,i]
            x[:,4:7,n+1] += h[:,np.newaxis]*a[n+1,i]*accel[:,:,i]
            x[:,7,n+1] -= h*a[n+1,i]*ioniz[:,i]
        x[:,:,n+1] += x[:,:,0]

    delta = np.zeros((len(X0), len(cols)))
    for i in range(6):
        delta[:,1:4] += bd[i]*x[:,4:7,i]
        delta[:,4:7] += bd[i]*accel[:,:,i]
        delta[:,7] += bd[i]*ioniz[:,i]
    delta = np.abs(h[:,np.newaxis]*delta)
    t1_ = Time.now()

    # Step 1: Determine initial acceleration and ionization rate
    a1, i1 = State(t1, x1, v1, output)
    
    # Step 2
    t2 = t1 - c2*h
    x2 = x1 + h*a21*v1
    v2 = v1 + h*a21*a1
    f2 = f1 - h*a21*i1
    a2, i2 = State(t2, x2, v2, output)

    # Step 3
    t3 = t1 - c3*h
    x3 = x1 + h*(a31*v1 + a32*v2)
    v3 = v1 + h*(a31*a1 + a32*a2)
    f3 = f1 - h*(a31*i1 + a32*i2)
    a3, i3 = State(t3, x3, v3, output)

    # Step 4
    t4 = t1 - c4*h
    x4 = x1 + h*(a41*v1 + a42*v2 + a43*v3)
    v4 = v1 + h*(a41*a1 + a42*a2 + a43*a3)
    f4 = f1 - h*(a41*i1 + a42*i2 + a43*i3)
    a4, i4 = State(t4, x4, v4, output)

    # Step 5
    t5 = t1 - c5*h
    x5 = x1 + h*(a51*v1 + a52*v2 + a53*v3 + a54*v4)
    v5 = v1 + h*(a51*a1 + a52*a2 + a53*a3 + a54*a4)
    f5 = f1 - h*(a51*i1 + a52*i2 + a53*i3 + a54*i4)
    a5, i5 = State(t5, x5, v5, output)

    # Step 6
    t6 = t1 - c6*h
    x6 = x1 + h*(a61*v1 + a62*v2 + a63*v3 + a64*v4 + a65*v5)
    v6 = v1 + h*(a61*a1 + a62*a2 + a63*a3 + a64*a4 + a65*a5)
    f6 = f1 - h*(a61*i1 + a62*i2 + a63*i3 + a64*i4 + a65*i5)
    a6, i6 = State(t6, x6, v6, output)

    # Step 7: Compute the final result
    t7 = t1 - c7*h
    x7 = x1 + h*(a71*v1 + a72*v2 + a73*v3 + a74*v4 + a75*v5 + a76*v6)
    v7 = v1 + h*(a71*a1 + a72*a2 + a73*a3 + a74*a4 + a75*a5 + a76*a6)
    f7 = f1 - h*(a71*i1 + a72*i2 + a73*i3 + a74*i4 + a75*i5 + a76*i6)

    # Step 8: Estimate the error
    # Note - difference from IDL version. b1d*loc.v where loc.v = final
    deltax = abs(h*(b1d*v1+b2d*v2+b3d*v3+b4d*v4+b5d*v5+b6d*v6))  # +
    # b7d*v7))
    deltav = abs(h*(b1d*a1+b2d*a2+b3d*a3+b4d*a4+b5d*a5+b6d*a6))  # +
    # b7d*a7))
    deltaf = abs(h*(b1d*i1+b2d*i2+b3d*i3+b4d*i4+b5d*i5+b6d*i6))  # +
    # b7d*i7))
    t2_ = Time.now()

    assert np.all(np.isclose(accel[:,:,0], a1.transpose()))
    assert np.all(np.isclose(ioniz[:,0], i1))

    assert np.all(np.isclose(x[:, 1:4, 1], x2.transpose()))
    assert np.all(np.isclose(x[:, 4:7, 1], v2.transpose()))
    assert np.all(np.isclose(x[:, 7, 1], f2))
    assert np.all(np.isclose(accel[:,:,1], a2.transpose()))
    assert np.all(np.isclose(ioniz[:,1], i2))
    assert np.all(np.isclose(x[:, 0, 1], t2))

    assert np.all(np.isclose(x[:, 0, 2], t3))
    assert np.all(np.isclose(x[:, 1:4, 2], x3.transpose()))
    assert np.all(np.isclose(x[:, 4:7, 2], v3.transpose()))
    assert np.all(np.isclose(x[:, 7, 2], f3))
    assert np.all(np.isclose(accel[:,:,2], a3.transpose()))
    assert np.all(np.isclose(ioniz[:,2], i3))

    assert np.all(np.isclose(x[:, 0, 3], t4))
    assert np.all(np.isclose(x[:, 1:4, 3], x4.transpose()))
    assert np.all(np.isclose(x[:, 4:7, 3], v4.transpose()))
    assert np.all(np.isclose(x[:, 7, 3], f4))
    assert np.all(np.isclose(accel[:,:,3], a4.transpose()))
    assert np.all(np.isclose(ioniz[:,3], i4))

    assert np.all(np.isclose(x[:, 0, 4], t5))
    assert np.all(np.isclose(x[:, 1:4, 4], x5.transpose()))
    assert np.all(np.isclose(x[:, 4:7, 4], v5.transpose()))
    assert np.all(np.isclose(x[:, 7, 4], f5))
    assert np.all(np.isclose(accel[:,:,4], a5.transpose()))
    assert np.all(np.isclose(ioniz[:,4], i5))

    assert np.all(np.isclose(x[:, 0, 5], t6))
    assert np.all(np.isclose(x[:, 1:4, 5], x6.transpose()))
    assert np.all(np.isclose(x[:, 4:7, 5], v6.transpose()))
    assert np.all(np.isclose(x[:, 7, 5], f6))
    assert np.all(np.isclose(accel[:,:,5], a6.transpose()))
    assert np.all(np.isclose(ioniz[:,5], i6))

    assert np.all(np.isclose(x[:, 0, 6], t7))
    assert np.all(np.isclose(x[:, 1:4, 6], x7.transpose()))
    assert np.all(np.isclose(x[:, 4:7, 6], v7.transpose()))
    assert np.all(np.isclose(x[:, 7, 6], f7))

    assert np.all(np.isclose(delta[:,1:4], deltax.transpose()))
    assert np.all(np.isclose(delta[:,4:7], deltav.transpose()))
    assert np.all(np.isclose(delta[:,7], deltaf))
    
    print(t1_-t0_)
    print(t2_-t1_)
    assert 0, 'Done'

    # Put frac back the way it should be
    result = x[:,:,6]
    result[:,:,6] = np.exp(result[:,:,6])
    
    f1 = np.exp(f1)
    f7 = np.exp(f7)

    return t7, x7, v7, f7, deltax, deltav, deltaf
