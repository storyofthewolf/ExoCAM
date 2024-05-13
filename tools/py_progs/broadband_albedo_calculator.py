"""
Created on Sat June 5 21:25:32 2021
@author: Vidya Venkatesan (vidyav1@uci.edu)

Project description: This code is a python version of ice_gcm.pro which calculates the albedo from stellar 
spectrum and spectrum of reflectance. 
Note: any spectrum file can be used along with any
reflectance spectrum provided the minGrid and maxGrid parameters are changed to reflect the appropriate units. 
Units are not important (as long as the wavelength units are the same for each file) for the 
final albedo calculations since it is normalized by both the wavelength and the stellar spectrum.
Required Inputs:
    stellarfile = file containing two columns:
                         1 - wavelength, microns
                         2 - flux, in any units (this is normalized
                             for the albedo calculation)
    albedofile  = file containing three columns (as generate by the USGS):
                         1 - wavelength, in microns
                         2 - reflectance at that wavelength (0 - 1)
                         3 - error, if any
"""
#Importing all the libraries
import numpy as np
import scipy.interpolate as ip
import matplotlib.pyplot as plt
import scipy.interpolate as interpol

#Set number of lines of header information to ignore
stel_nh = 2
alb_nh = 0

#Conversion factor
#wavelengths must be converted to microns
#check header information for given SED
#1.0e-4 angstroms to microns
stel_convert_wavl = 1.0e-4
alb_convert_wavl = 1.0

#Reading in the stellar file 
stellar_file = "/discover/nobackup/etwolf/models/ExoRT/data/solar/raw/TOI700_SED_HST.txt"
with open(stellar_file, 'r') as f: 
    lines = f.readlines()
    header = lines[0:stel_nh]
    lamda= [float(line.split()[0]) * stel_convert_wavl for line in lines[stel_nh:]]
    flux= [float(line.split()[1]) for line in lines[stel_nh:]] 
    
#Reading in the Albedo file for surface
#albedo_file = "../spectral_albedos/50%Mixture.txt"
#albedo_file = "../spectral_albedos/snow100um.txt"
albedo_file = "../spectral_albedos/Bluemarineice.txt"
with open(albedo_file, 'r') as f:
    lines = f.readlines()
    header = lines[0:alb_nh]
    wave = [float(line.split()[0]) * alb_convert_wavl for line in lines[alb_nh:]]
    albedo = [float(line.split()[1]) for line in lines[alb_nh:]]

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Using stellar spectral file...")
print("  ", stellar_file)
print("Using albedo file...")
print("  ", albedo_file)
    
#find the minimum and maximum in both stellar and albedo files and match them for interpolation purpose
a = np.min(wave)
b = np.max(wave)
c = np.min(lamda)
d = np.max(lamda)
start_g = max(a,c)
end_g = min(b,d)

# Checking maximum and minimum grid is important else code result is NAN. Default = 2.5, Warm ice 2.3
mingrid = start_g
maxgrid = end_g


#Set up the grid
ngrid = 10000
n = np.arange(10000)
k = maxgrid - mingrid
wavelengthgrid = (mingrid + k * n / (ngrid-1))
dlamda = (np.max(wavelengthgrid)- np.min(wavelengthgrid))/(len(wavelengthgrid)-1)

# Interpolation
stellarInterpolate4 = np.interp(wavelengthgrid,lamda,flux)
temp = interpol.CubicSpline(wave,albedo)
albedoInterpolate = temp(wavelengthgrid)

#Normalization 
flux4 = stellarInterpolate4/np.max(stellarInterpolate4)

#More calculations for scaling 
Snorm = 1360
total_sed4 = sum(stellarInterpolate4*dlamda)
ScaleFac = Snorm/total_sed4
scaledflux4 = total_sed4*ScaleFac
new_total_sed4 = scaledflux4
new_flux4 = stellarInterpolate4*ScaleFac
total_new_flux_4 = sum(new_flux4 * dlamda)

#Setting up cuts between IR and Vis
#cutoff wavelengths in microns
cutoff = 0.76923 
ind_v = np.where(wavelengthgrid[:]<cutoff)
ind_i = np.where(wavelengthgrid[:]> cutoff )

print("Using cuttoff wavelength (microns) ...")
print("  ", cutoff)

#use this to calculate bond albedo
ind_w = np.where(wavelengthgrid[:]>0)

#Multiplying reflectance spectrum to normalized stellar specrtrum
albedo_sed_M_v = (albedoInterpolate[ind_v]) * (new_flux4[ind_v]) *dlamda
albedo_sed_M_i = (albedoInterpolate[ind_i]) * (new_flux4[ind_i]) *dlamda
albedo_sed = (albedoInterpolate[ind_w])*(new_flux4[ind_w]) *dlamda

#Calculate the sum
total_albedo_sed_M_v = sum(albedo_sed_M_v)
total_albedo_sed_M_i = sum(albedo_sed_M_i)
total_albedo_sed = sum(albedo_sed)

total_sed_M_v = sum(new_flux4[ind_v]*dlamda)
total_sed_M_i = sum(new_flux4[ind_i]*dlamda)
total_sed = sum(new_flux4[ind_w]*dlamda)

#Calculate albedo for IR, Vis, and bond albedo for a given host star SED
albedo_M_v = total_albedo_sed_M_v/total_sed_M_v
albedo_M_i = total_albedo_sed_M_i/total_sed_M_i
Bond_albedo = total_albedo_sed/total_sed

#Print the numbers
print('Albedo Vis =', albedo_M_v)
print('Albedo IR =', albedo_M_i)
print('Broadband Albedo =', Bond_albedo)
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
