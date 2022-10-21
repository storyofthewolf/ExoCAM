#====================================================================
# Author: Eric T. Wolf
# Translator: Russell Deitrick
# Date:  some time long ago
# needs fixing
#====================================================================
# Change number of vertical levels in a CAM initial condition file.
# This is achieved starting with a 66 level WACCM initial condition
# file and cutting it down to size.
# Adapted for Python by R. Deitrick from
# idl_progs/changevert_cesm.pro (E.T. Wolf)
#====================================================================

import netCDF4 as nc
import numpy as np
import argparse
import pathlib
import exocampy_tools as exo

exocam_path = '/home/deitrr/ExoCAM/'
ccsm_inputdata_path = '/home/deitrr/scratch/ccsm/inputdata/atm/cam/inic/fv/'

parser = argparse.ArgumentParser()
parser.add_argument('fname_out',nargs=1,help='Name of new IC file')
parser.add_argument('-n','--num_lev', nargs=1, default=[36], help='Number of levels in new IC file')
parser.add_argument('-ic','--input_IC_file', nargs=1, default=['cami-mam3_0000-01-01_1.9x2.5_L30_c090306.nc'],
                     help='input climate file')
parser.add_argument('-w','--overwrite', action='store_true', help = 'force overwrite of output files')
args = parser.parse_args()


#================================
#=======  name of new file ======
#================================
fname_out = args.fname_out[0]
nlev_out = args.num_lev[0]
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
ncid = nc.Dataset(ccsm_inputdata_path + clim_fname_in, 'r')

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

import pdb; pdb.set_trace()

