import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

fcon = 'control_L40.cam.i.0009-01-01-00000.nc'
control = nc.Dataset(fcon,'r')

fp5 = 'ic_P0.5bar_L40_ic.nc'
datp5 = nc.Dataset(fp5,'r')

f1 = 'ic_P1bar_L40_ic.nc'
dat1 = nc.Dataset(f1,'r')

f2 = 'ic_P2bar_L40_ic.nc'
dat2 = nc.Dataset(f2,'r')

#fig, axes = plt.subplots(ncols=2,nrows=2,figsize=(12,10))
fig = plt.figure(figsize=(12,10))

fields = ['T', 'US', 'VS', 'Q']

for i in np.arange(len(fields)):
    ax = plt.subplot(2,2,i+1)
    ax.semilogy(control[fields[i]][0,:,22,0], control['lev'][:]*control['P0'][0]/1e5,  'k--', lw=2)
    ax.plot(dat2[fields[i]][0,:,22,0], dat2['lev'][:]*dat2['P0'][0]/1e5, 'b-')
    ax.plot(dat1[fields[i]][0,:,22,0], dat1['lev'][:]*dat1['P0'][0]/1e5,  'g-')
    ax.plot(datp5[fields[i]][0,:,22,0], datp5['lev'][:]*datp5['P0'][0]/1e5, 'r-')
    ax.invert_yaxis()
    ax.set_xlabel(fields[i])
    ax.set_ylabel('P')

plt.show()
