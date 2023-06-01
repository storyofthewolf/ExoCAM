
pro haze_emission
;--------------------------------------------------------------------
;AUTHOR:  WOLF, E. T.
;
;PURPOSE:  Create source file in netcdf format describing the
;production of organic haze aerosols in the atmosphere of the Early
;Earth.  A lognormal distribution is assumed.
;Haze aersol density and production rate determined by TRAINER, 2006
;-------------------------------------------------------------------
;
;Confirmed to meteor_smoke_kalashnikova.nc file within a factor of
;1.2.  I believe this small discrepency most likely has to do with my
;level thickness being grossly approximated using a constant scale
;height at all altitudes
;
;Note that any deviations in the haze aerosol production function seem
;to be corrected by setting carma_emis_total [KT yr-1] to the expected
;production rate.  This is required to be set in WACCM/CARMA, carma.F90



print, ""
print, "__________________________________________________________"
print, ""

print, "_______________________PARAMETERS_________________________"

;TEST PARAMETERS to duplicate meteor_smoke_kalashnikova.nc
;G=1.6E10
;rho=2.0
;initial_radius=1.3e-7

G=1.0E14                                          ;g yr-1    annual production rate
rho=0.64                                         ;g cm-3    effective density
initial_radius=1.0e-7                             ;cm        radius of spherical monomers 
initial_mass=4.0*!pi/3*(initial_radius)^3*rho     ;g         mass per particle
N=G/initial_mass/(365.0*24.0*60.0*60.0)           ;# s-1    number of particles emitted per second
R_earth=6.378E8                                   ;cm       radius of the earth
G_KT=G/1000.0/1000.0/1000.0

print, "G [g yr-1] = ",G
print, "G_KT [kt] = ",G_KT
print, "rho [g cm-3] = ",rho
print, "initial_radius [cm] = ",initial_radius
print, "initial_mass [g] = ",initial_mass
print, "N [# s-1] = ",N 

;Define lognormal distribution for haze aerosol prodction.
;Production of haze aersol is the largest where the photolysis of
;methane is the greatest.  This coincides with the largest
;concentrations of C2H4, C2H2, which are precursors to haze aerosol
;formation.  PAVLOV, 2001

;haze aerosol production function in height coordinates

;TEST PARAMETERS to duplicate meteor_smoke_kalashnikova.nc
;z_bar=88.0
;sigma_Z=5.0
;z=findgen(120)

z_bar=65.0                              ;km         height of maximum haze aerosol production
sigma_z=6.0                             ;km           standard deviation
z=findgen(100)                          ;km         height vector (for plotting purposes)             
  
norm_z=N/(sqrt(2*!pi)*sigma_z)*exp(-(z-z_bar)^2/(2*sigma_z^2))
;number s-1 km-1, normal distribution


;plot, norm_z, z

;check consistency
print, "consistency check: ", total(norm_z)  
print, ""

print,"________________________BIN_PARAMETERS_____________________"

rmrat=2.0                                ;Ratio of particle mass between successive bins
NBIN=32                                  ;Number of size bins 
m_bin=initial_mass

FOR i=2,NBIN DO BEGIN
   IF (i eq 2) THEN print, "Mass Bin ",i-1," = ",m_bin," g :  Radius Bin ",i-1," = ",initial_radius," cm"
   m_bin=m_bin*rmrat
   r_bin=(3.0*m_bin/(4.0*!pi*rho))^(1.0/3)
 ;  print, "Mass Bin ",i," = ",m_bin," g :  Radius Bin ",i," = ",r_bin," cm"
ENDFOR

print,""
print,"____________________DISTRIBUTION_PARAMETERS_________________"
print,""

P_0=101300.0                              ;Pa           reference pressure, SLP
H=7.0                                    ;km           scale height 
z_prime=findgen(51)+40.0                  ;km           height indexing (50 levels used)
p_prime=P_0*exp(-z_prime/H)              ;Pa           pressure indexing
delta_z=1.0e5                            ;cm km-1           differential height level (1 km)
max_lev=max(z_prime)                     ;km   
min_lev=min(z_prime)                     ;km
delta_atm=(max_lev-min_lev)*1.0e5        ;cm
print, "delta_atm = ",delta_atm/1.0e5
    
N_C2=N/(4.0*!pi*R_earth^2)      ;[# cm-2 s-1]     column number of particles emitted per cm2

                                ;recall from rules of normal
                                ;distribution that N_C2 =int(norm dz)
                                ;thus the integral of the normal
                                ;distribution function over all
                                ;heights should equal N_C2
   
print, "N_C2 [# s-1 cm-2] = ",N_C2

;normal distribution
;sigma_z is divided by 1.0e5 to convert from km > cm. N_C2 ~ [# cm-2 s-1] and therefore
;norm ~ [# cm-3 s-1] binned per kilometer

norm=N_C2/(sqrt(2*!pi)*sigma_z*delta_z)*exp(-(z_prime-z_bar)^2/(2*sigma_z^2))  ;number per cm-3 binned per kilo



FOR j=0,50 DO BEGIN
   print, j," ","Z = ",z_prime(j)," km  | P(Z) = ",p_prime(j)," Pa  |  F = ",norm(j)," # cm-3 s-1"
ENDFOR

print, "consistency check: ", total(norm*delta_z)

;plot, norm, p_prime, /ylog, yrange=[400,0.1]

N_CHECK=0.0
FOR j=0,50 DO BEGIN
   N_CHECK=N_CHECK+norm(j)*(4.0*!pi*R_Earth^2.0*delta_Z)
;   print, norm(j),N_CHECK
ENDFOR

print, "N_CHECK = ", N_CHECK

print, "___________________WROTE_NETCDF____________________________"
print, ""

filename="early_earth_haze.nc"
;filename="smoke_test.nc"
id = NCDF_CREATE(filename, /CLOBBER)
dim1 = NCDF_DIMDEF(id,'lev',51)

varid_1 = NCDF_VARDEF(id, 'MHAZE',[dim1], /float)
varid_2 = NCDF_VARDEF(id, 'lev',[dim1], /float)

NCDF_ATTPUT, id, varid_1, "title", "haze production function"
NCDF_ATTPUT, id, varid_1, "units", "# cm-3 s-1"
NCDF_ATTPUT, id, varid_1, "_fillValue", 0.0
NCDF_ATTPUT, id, varid_2, "title", "pressure"
NCDF_ATTPUT, id, varid_2, "units", "Pa"
NCDF_ATTPUT, id, varid_2, "_fillValue",0.0
NCDF_CONTROL, id, /ENDEF

NCDF_VARPUT,id,varid_1, norm
NCDF_VARPUT,id,varid_2, p_prime

NCDF_CLOSE, id


print, ""
print, "______________________DONE__________________________________"
print, ""

end
