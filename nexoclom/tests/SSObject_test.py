import numpy as np
import astropy.units as u
from solarsystemMB import SSObject, planet_dist
from mymath import minmaxmean
import matplotlib.pyplot as plt

# Create existing objects
mercury = SSObject('Mercury')
io = SSObject('iO')

# Create a non-existing object
#bad = SSObject('asdf')

print('Mercury')
print(mercury)
print()
print('Io')
print(io)
print()
#print('asdf')
#print(bad)

# Find some planet distances
taa = np.linspace(0, 2*np.pi, 73)*u.rad
taa2 = np.linspace(0, 360, 73)*u.deg
taa3 = np.linspace(0, 2*np.pi, 73)

r, v_r = planet_dist('Mercury', taa=taa)
r2, v_r2 = planet_dist(mercury, taa=taa2)
r3, v_r3 = planet_dist(mercury, taa=taa3)

print(minmaxmean(r-r2))
print(minmaxmean(v_r-v_r2))
print(minmaxmean(r-r3))
print(minmaxmean(v_r-v_r3))

fix, ax = plt.subplots(ncols=2)
ax[0].plot(taa, r)
ax[1].plot(taa, v_r)
plt.show()
