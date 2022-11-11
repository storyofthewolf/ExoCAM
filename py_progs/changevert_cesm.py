#====================================================================
# Author: Eric T. Wolf
# Translator: Russell Deitrick
# Date:  some time long ago
# needs fixing?
#====================================================================
# Change number of vertical levels in a CAM initial condition file.
# This is achieved starting with a 66 level WACCM initial condition
# file and cutting it down to size.
# Adapted for Python by R. Deitrick from
# idl_progs/changevert_cesm.pro (E.T. Wolf)
#====================================================================

import netCDF4 as nc
import numpy as np
from scipy.interpolate import interp1d
import argparse
import pathlib
import exocampy_tools as exo

#user needs to edit these!
exocam_path = '/home/deitrr/projects/def-czg/deitrr/ExoCAM/'
ccsm_inputdata = '/home/deitrr/projects/def-czg/deitrr/ccsm/inputdata/'

ccsm_inputdata_paths = [ccsm_inputdata + 'atm/cam/inic/fv/',
                        ccsm_inputdata + 'atm/waccm/ic/',
			exocam_path + '/cesm1.2.1/initial_files/cam_aqua_fv/',
                        exocam_path + '/cesm1.2.1/initial_files/cam_aqua_se/',
                        exocam_path + '/cesm1.2.1/initial_files/cam_land_fv/',
                        exocam_path + '/cesm1.2.1/initial_files/cam_mixed_fv/',
                        exocam_path + '/cesm1.2.1/initial_files/mars/atm/']

parser = argparse.ArgumentParser()
parser.add_argument('fname_out',nargs=1,help='Name of new IC file')
parser.add_argument('-n','--num_lev', nargs=1, default=[40], help='Number of levels in new IC file')
parser.add_argument('-ic','--input_IC_file', nargs=1, default=['cami-mam3_0000-01-01_1.9x2.5_L30_c090306.nc'],
                     help='input climate file')
parser.add_argument('-w','--overwrite', action='store_true', help = 'force overwrite of output files')
args = parser.parse_args()


#================================
#=======  name of new file ======
#================================
fname_out = args.fname_out[0]
nlev_out = int(args.num_lev[0])
nilev_out = nlev_out + 1


#read in appropriate lev, hyai for N>26
lev_fname_new = exocam_path + 'cesm1.2.1/initial_files/other/oxygen_CE.cam2.avg.nc'
ncid = nc.Dataset(lev_fname_new,'r')
lev_new = ncid['lev'][:]
ilev_new = ncid['ilev'][:]
lat_new = ncid['lat'][:]
lon_new = ncid['lon'][:]
hyai_new = ncid['hyai'][:]
hybi_new = ncid['hybi'][:]
hyam_new = ncid['hyam'][:]
hybm_new = ncid['hybm'][:]
ncid.close()

nlev_new = len(lev_new)
nilev_new = len(ilev_new)
nlat_new = len(lat_new)
nlon_new = len(lon_new)


#===========================================
#=======  name file with climate data ======
#===========================================
clim_fname_in = args.input_IC_file[0]

filefound = False
for fdir in ccsm_inputdata_paths:
  if pathlib.Path(fdir+clim_fname_in).exists():
    ncid = nc.Dataset(fdir + clim_fname_in, 'r')
    filefound = True
    break
if filefound == False:
  raise IOError('File %s not found!'%clim_fname_in)

lev_old = ncid['lev'][:]
ilev_old = ncid['ilev'][:]
lat = ncid['lat'][:]
lon = ncid['lon'][:]
hyai_old = ncid['hyai'][:]
hybi_old = ncid['hybi'][:]
hyam_old = ncid['hyam'][:]
hybm_old = ncid['hybm'][:]

nlev_old = len(lev_old)
nilev_old = len(ilev_old)
nlat = len(lat)
nlon = len(lon)
nslat = nlat-1
nslon = nlon

CLDICE_old = ncid['CLDICE'][:]
CLDLIQ_old = ncid['CLDLIQ'][:]
Q_old = ncid['Q'][:]
T_old = ncid['T'][:]
US_old = ncid['US'][:]
VS_old = ncid['VS'][:]

#stuff that doesn't need to be changed
P0 = ncid['P0'][:]
slat = ncid['slat'][:]
slon = ncid['slon'][:]
w_stag = ncid['w_stag'][:]
time = ncid['time'][:]
time_bnds = ncid['time_bnds'][:]
date_written = ncid['date_written'][:]
time_written = ncid['time_written'][:]
ntrm = ncid['ntrm'][:]
ntrn = ncid['ntrn'][:]
ntrk = ncid['ntrk'][:]
ndbase = ncid['ndbase'][:]
nsbase = ncid['nsbase'][:]
nbdate = ncid['nbdate'][:]
nbsec = ncid['nbsec'][:]
mdt = ncid['mdt'][:]
gw = ncid['gw'][:]
ndcur = ncid['ndcur'][:]
nscur = ncid['nscur'][:]
date = ncid['date'][:]
datesec = ncid['datesec'][:]
nsteph = ncid['nsteph'][:]
ICEFRAC = ncid['ICEFRAC'][:]
PS = ncid['PS'][:]
SICTHK = ncid['SICTHK'][:]
SNOWHICE = ncid['SNOWHICE'][:]
TS1 = ncid['TS1'][:]
TS2 = ncid['TS2'][:]
TS3 = ncid['TS3'][:]
TS4 = ncid['TS4'][:]
TSICE = ncid['TSICE'][:]
ncid.close()

# create pressure grids from hybrid sigma coordinates
lev_P_new, ilev_P_new = exo.hybrid2pressure(nlon,nlat,nlev_new,PS.squeeze(),P0,hyam_new,hybm_new,hyai_new,hybi_new)

lev_P_old, ilev_P_old =  exo.hybrid2pressure(nlon,nlat,nlev_old,PS.squeeze(),P0,hyam_old,hybm_old,hyai_old,hybi_old)


#-----------------------------------------
#----- Interpolate relevant fields -------
#-----------------------------------------

#out variables
T_out = np.zeros((nlev_out,nlat,nlon))
CLDICE_out = np.zeros_like(T_out)
CLDLIQ_out = np.zeros_like(T_out)
CLOUD_out = np.zeros_like(T_out)
Q_out = np.zeros_like(T_out)
US_out = np.zeros((nlev_out,nslat,nlon))
VS_out = np.zeros_like(T_out)

lev_out = np.zeros(nlev_out)
hyam_out = np.zeros(nlev_out)
hybm_out = np.zeros(nlev_out)
ilev_out = np.zeros(nilev_out)
hyai_out = np.zeros(nilev_out)
hybi_out = np.zeros(nilev_out)

n = 66 - nlev_out  #index of lowest pressure (highest altitude) of new vertical grid
for x in np.arange(nlon):
  for y in np.arange(nlat):
    T_out[:,y,x] = interp1d(lev_P_old[:,y,x].data, T_old[0,:,y,x].data, fill_value='extrapolate')(lev_P_new[n:,y,x].data)
    CLDICE_out[:,y,x] = interp1d(lev_P_old[:,y,x].data, CLDICE_old[0,:,y,x].data, fill_value='extrapolate')(lev_P_new[n:,y,x].data)
    CLDLIQ_out[:,y,x] = interp1d(lev_P_old[:,y,x].data, CLDLIQ_old[0,:,y,x].data, fill_value='extrapolate')(lev_P_new[n:,y,x].data)
    Q_out[:,y,x] = interp1d(lev_P_old[:,y,x].data, Q_old[0,:,y,x].data, fill_value='extrapolate')(lev_P_new[n:,y,x].data)
    VS_out[:,y,x] = interp1d(lev_P_old[:,y,x].data, VS_old[0,:,y,x].data, fill_value='extrapolate')(lev_P_new[n:,y,x].data)
    if y < nslat:
      US_out[:,y,x] = interp1d(lev_P_old[:,y,x].data, US_old[0,:,y,x].data, fill_value='extrapolate')(lev_P_new[n:,y,x].data)

for i in np.arange(nilev_out-1,-1,-1):
  j = nilev_new - nilev_out + i
  hyai_out[i] = hyai_new[j]
  hybi_out[i] = hybi_new[j]
  ilev_out[i] = ilev_new[j]

for i in np.arange(0,nlev_out):
  j = nlev_new - nlev_out + i
  hyam_out[i] = hyam_new[j]
  hybm_out[i] = hybi_new[j]
  lev_out[i] = lev_new[j]


if not pathlib.Path(fname_out).exists() or args.overwrite:
  print("Creating file '%s'..."%fname_out)
  id = nc.Dataset(fname_out,'w')
  dim1 = id.createDimension("lat",nlat)
  dim2 = id.createDimension("lon",nlon)
  dim3 = id.createDimension("slat",nslat)
  dim4 = id.createDimension("slon",nlon)  #same as lon, since it wraps around
  dim5 = id.createDimension("lev",nlev_out)
  dim6 = id.createDimension("ilev",nilev_out)
  dim10 = id.createDimension("time",1)
  dim11 = id.createDimension("nbnd",2)
  dim12 = id.createDimension("chars",8)

  varid1 = id.createVariable("P0","f8")
  varid1.long_name = "reference pressure"
  varid1.units = "Pa"
  varid1[:] = P0

  varid2 = id.createVariable("lat","f8",(dim1.name,))
  varid2.long_name = "latitude"
  varid2.units = "degrees_north"
  varid2[:] = lat

  varid3 = id.createVariable("lon","f8",(dim2.name,))
  varid3.long_name = "longitude"
  varid3.units = "degrees_east"
  varid3[:] = lon

  varid4 = id.createVariable("slat","f8",(dim3.name,))
  varid4.long_name = "staggered latitude"
  varid4.units = "degrees_north"
  varid4[:] = slat

  varid5 = id.createVariable("slon","f8",(dim4.name,))
  varid5.long_name = "staggered longitude"
  varid5.units = "degrees_east"
  varid5[:] = slon

  varid6 = id.createVariable("w_stag","f8",(dim3.name,))
  varid6.long_name = "staggered latitude weights"
  varid6[:] = w_stag

  varid7 = id.createVariable("lev","f8",(dim5.name,))
  varid7.long_name = "hybrid level at midpoints (1000*(A+B))"
  varid7.units = "level"
  varid7.positive = "down"
  varid7.standard_name = "atmosphere_hybrid_sigma_pressure_coordinate"
  varid7.formula = "a: hyam b: hybm p0: P0 ps: PS"
  varid7[:] = lev_out

  varid8 = id.createVariable("ilev","f8",(dim6.name,))
  varid8.long_name = "hybrid level at interfaces (1000*(A+B))"
  varid8.units = "level"
  varid8.positive = "down"
  varid8.standard_name = "atmosphere_hybrid_sigma_pressure_coordinate"
  varid8.formula = "a: hyai b: hybi p0: P0 ps: PS"
  varid8[:] = ilev_out

  varid12 = id.createVariable("time","f8",(dim10.name,))
  varid12.long_name = "time"
  varid12.units = "days since 1995-01-01 00:00:00"
  varid12.calendar = "noleap"
  varid12.bounds = "time_bnds"
  varid12[:] = time

  varid13 = id.createVariable("time_bnds","f8",(dim11.name,))
  varid13.long_name = "time interval endpoints"
  varid13[:] = time_bnds[:2]

  varid14 = id.createVariable("date_written","S1",(dim10.name,dim12.name))
  varid14[:] = date_written

  varid15 = id.createVariable("time_written","S1",(dim10.name,dim12.name))
  varid15[:] = time_written

  varid16 = id.createVariable("ntrm","i4")
  varid16.long_name = "spectral truncation parameter M"
  varid16[:] = ntrm

  varid17 = id.createVariable("ntrn","i4")
  varid17.long_name = "spectral truncation parameter N"
  varid17[:] = ntrn

  varid18 = id.createVariable("ntrk","i4")
  varid18.long_name = "spectral truncation parameter K"
  varid18[:] = ntrk

  varid19 = id.createVariable("ndbase","i4")
  varid19.long_name = "base day"
  varid19[:] = ndbase

  varid20 = id.createVariable("nsbase","i4")
  varid20.long_name = "seconds of base day"
  varid20[:] = nsbase

  varid21 = id.createVariable("nbdate","i4")
  varid21.long_name = "base date (YYYYMMDD)"
  varid21[:] = nbdate

  varid22 = id.createVariable("nbsec","i4")
  varid22.long_name = "seconds of base date"
  varid22[:] = nbsec

  varid23 = id.createVariable("mdt","i4")
  varid23.long_name = "timestep"
  varid23.units = "s"
  varid23[:] = mdt

  varid24 = id.createVariable("hyai","f8",(dim6.name,))
  varid24.long_name = "hybrid A coefficient at layer interfaces"
  varid24[:] = hyai_out

  varid25 = id.createVariable("hybi","f8",(dim6.name,))
  varid25.long_name = "hybrid B coefficient at layer interfaces"
  varid25[:] = hybi_out

  varid26 = id.createVariable("hyam","f8",(dim5.name,))
  varid26.long_name = "hybrid A coefficient at layer midpoints"
  varid26[:] = hyam_out

  varid27 = id.createVariable("hybm","f8",(dim5.name,))
  varid27.long_name = "hybrid B coefficient at layer midpoints"
  varid27[:] = hybm_out

  varid28 = id.createVariable("gw","f8",(dim1.name,))
  varid28.long_name = "gauss_weights"
  varid28[:] = gw

  varid29 = id.createVariable("ndcur","i4",(dim10.name,))
  varid29.long_name = "current day (from base day)"
  varid29[:] = ndcur

  varid30 = id.createVariable("nscur","i4",(dim10.name,))
  varid30.long_name = "current seconds of current day"
  varid30[:] = nscur

  varid31 = id.createVariable("date","i4",(dim10.name,))
  varid31.long_name = "current date (YYYYMMDD)"
  varid31[:] = date

  varid33 = id.createVariable("datesec","i4",(dim10.name,))
  varid33.long_name = "current seconds of current date"
  varid33[:] = datesec

  varid34 = id.createVariable("nsteph","i4",(dim10.name,))
  varid34.long_name = "current timestep"
  varid34[:] = nsteph

  #order dimensions (time,level,latitude,longitude)
  varid57 = id.createVariable("CLDICE","f8",(dim10.name,dim5.name,dim1.name,dim2.name))
  varid57.units = "kg/kg"
  varid57.long_name = "Grid box averaged ice condensate amount"
  varid57[:] = CLDICE_out

  varid58 = id.createVariable("CLDLIQ","f8",(dim10.name,dim5.name,dim1.name,dim2.name))
  varid58.units = "kg/kg"
  varid58.long_name = "Grid box averaged liquid condensate amount"
  varid58[:] = CLDLIQ_out

  varid76 = id.createVariable("ICEFRAC","f8",(dim10.name,dim1.name,dim2.name))
  varid76.units = "fraction"
  varid76.long_name = "Fraction of sfc area covered by sea-ice"
  varid76[:] = ICEFRAC

  varid99 = id.createVariable("PS","f8",(dim10.name,dim1.name,dim2.name))
  varid99.units = "Pa"
  varid99.long_name = "surface pressure"
  varid99[:] = PS

  varid100 = id.createVariable("Q","f8",(dim10.name,dim5.name,dim1.name,dim2.name))
  varid100.units = "kg/kg"
  varid100.long_name = "Specific humidity"
  varid100[:] = Q_out

  varid103 = id.createVariable("SICTHK","f8",(dim10.name,dim1.name,dim2.name))
  varid103.units = "m"
  varid103.long_name = "Sea ice thickness"
  varid103[:] = SICTHK

  varid104 = id.createVariable("SNOWHICE","f8",(dim10.name,dim1.name,dim2.name))
  varid104.units = "m"
  varid104.long_name = "Water equivalent snow depth"
  varid104[:] = SNOWHICE

  varid105 = id.createVariable("T","f8",(dim10.name,dim5.name,dim1.name,dim2.name))
  varid105.units = "K"
  varid105.long_name = "Temperature"
  varid105[:] = T_out

  varid110 = id.createVariable("TS1","f8",(dim10.name,dim1.name,dim2.name))
  varid110.units = "K"
  varid110.long_name = "TS1 subsoil temperature"
  varid110[:] = TS1

  varid111 = id.createVariable("TS2","f8",(dim10.name,dim1.name,dim2.name))
  varid111.units = "K"
  varid111.long_name = "TS2 subsoil temperature"
  varid111[:] = TS2

  varid112 = id.createVariable("TS3","f8",(dim10.name,dim1.name,dim2.name))
  varid112.units = "K"
  varid112.long_name = "TS3 subsoil temperature"
  varid112[:] = TS3

  varid113 = id.createVariable("TS4","f8",(dim10.name,dim1.name,dim2.name))
  varid113.units = "K"
  varid113.long_name = "TS4 subsoil temperature"
  varid113[:] = TS4

  varid114 = id.createVariable("TSICE","f8",(dim10.name,dim1.name,dim2.name))
  varid114.units = "K"
  varid114.long_name = "Ice temperature"
  varid114[:] = TSICE 

  varid117 = id.createVariable("US","f8",(dim10.name,dim5.name,dim3.name,dim2.name))
  varid117.units = "m/s"
  varid117.units = "Zonal wind staggered"
  varid117[:] = US_out

  varid119 = id.createVariable("VS","f8",(dim10.name,dim5.name,dim1.name,dim4.name))
  varid119.units = "m/s"
  varid119.units = "Meridional wind staggered"
  varid119[:] = VS_out

  id.close()

