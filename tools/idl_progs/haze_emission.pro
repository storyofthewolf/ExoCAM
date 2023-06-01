pro haze_emission
;--------------------------------------------------------------------
;AUTHOR:  WOLF, E. T.
;
;PURPOSE:  Create production rate profile
;


do_write  = 1
do_zenith = 1

;input files
nfiles=18
filenames=strarr(nfiles)

nlev = 200
nzen = 18
path = "/gpfsm/dnb53/etwolf/atmos_haze_prod"
folder = "/Wolf_hazmat_prox_mCO2_1e-2_mCH4_3.70E-03"
filenames(17) = path + folder + "/hazmat_prox_mCH4_3.70E-03_zenith_015_haze_rate.dat"
filenames(16) = path + folder + "/hazmat_prox_mCH4_3.70E-03_zenith_060_haze_rate.dat"
filenames(15) = path + folder + "/hazmat_prox_mCH4_3.70E-03_zenith_105_haze_rate.dat"
filenames(14) = path + folder + "/hazmat_prox_mCH4_3.70E-03_zenith_150_haze_rate.dat"
filenames(13) = path + folder + "/hazmat_prox_mCH4_3.70E-03_zenith_195_haze_rate.dat"
filenames(12) = path + folder + "/hazmat_prox_mCH4_3.70E-03_zenith_240_haze_rate.dat"
filenames(11) = path + folder + "/hazmat_prox_mCH4_3.70E-03_zenith_285_haze_rate.dat"
filenames(10) = path + folder + "/hazmat_prox_mCH4_3.70E-03_zenith_330_haze_rate.dat"
filenames(9) = path + folder + "/hazmat_prox_mCH4_3.70E-03_zenith_420_haze_rate.dat"
filenames(8) = path + folder + "/hazmat_prox_mCH4_3.70E-03_zenith_465_haze_rate.dat"
filenames(7) = path + folder + "/hazmat_prox_mCH4_3.70E-03_zenith_510_haze_rate.dat"
filenames(6) = path + folder + "/hazmat_prox_mCH4_3.70E-03_zenith_555_haze_rate.dat"
filenames(5) = path + folder + "/hazmat_prox_mCH4_3.70E-03_zenith_600_haze_rate.dat"
filenames(4) = path + folder + "/hazmat_prox_mCH4_3.70E-03_zenith_645_haze_rate.dat"
filenames(3) = path + folder + "/hazmat_prox_mCH4_3.70E-03_zenith_690_haze_rate.dat"
filenames(2) = path + folder + "/hazmat_prox_mCH4_3.70E-03_zenith_780_haze_rate.dat"
filenames(1) = path + folder + "/hazmat_prox_mCH4_3.70E-03_zenith_825_haze_rate.dat"
filenames(0) = path + folder + "/hazmat_prox_mCH4_3.70E-03_zenith_870_haze_rate.dat"


pressure = fltarr(nlev, nfiles)
prod     = fltarr(nlev, nfiles)
header = strarr(1)
press_temp   = fltarr(nlev)
prod_temp    = fltarr(nlev)

for z=0, nfiles-1 do begin

  data = fltarr(3,nlev)
  OPENR,lun,filenames(z),/GET_LUN
  READF,lun,header
  READF,lun,data
  FREE_LUN,lun

  press_temp(*) = data(0,*) * 1.0e5 ; convert bars to pascals
  prod_temp(*) = data(2,*)         ; g/cm3/s 
  press_temp   = reverse(press_temp)
  prod_temp    = reverse(prod_temp)

  pressure(*,z) = press_temp(*)
  prod(*,z)     = prod_temp(*)

print, prod(*,z)

endfor




;for i=0,nlev-1 do begin 
;  print, i, pressure(i), prod(i)
;endfor

if (do_zenith eq 1) then begin

;  zen = [0, 10, 20, 30, 40, 50, 60, 70, 89]
;  zen = reverse(zen)
;  nz = n_elements(zen)
;  cosz = fltarr(nz)
;  prod_zen = fltarr(nlev, nzen)
;  for z=0, nzen-1 do begin 
;    prod_zen(*,z) = prod(*) * cos(zen(z)*!pi/180)
;    cosz(z) = cos(zen(z)*!pi/180)
;  endfor    

zen = [ 1.5,  6.0, 10.5, 15.0, 19.5, 24.0, 28.5, 33.0, $
       37.5, 42.0, 46.5, 51.0, 55.5, 60.0, 64.5, 69.0, $
       73.5, 78.0, 82.5, 87.0 ]
zen = reverse(zen)
cosz = fltarr(nfiles)
for z=0, nfiles-1 do cosz(z) = cos(zen(z)*!pi/180) 

endif



print, "___________________WROTE_NETCDF____________________________"
print, ""

filename="early_earth_haze_zen_atmos.nc"

id = NCDF_CREATE(filename, /CLOBBER)
dim1 = NCDF_DIMDEF(id,'lev',nlev)
if (do_zenith) then dim2 = NCDF_DIMDEF(id,'zenith',nzen)

if (do_zenith) then begin
  varid_1 = NCDF_VARDEF(id, 'MHAZE',[dim1, dim2], /float)
  varid_3 = NCDF_VARDEF(id, 'cosz',[dim2], /float)
endif else begin
  varid_1 = NCDF_VARDEF(id, 'MHAZE',[dim1], /float)
endelse
varid_2 = NCDF_VARDEF(id, 'lev',[dim1], /float)


NCDF_ATTPUT, id, varid_1, "title", "haze production function"
NCDF_ATTPUT, id, varid_1, "units", "g cm-3 s-1"
NCDF_ATTPUT, id, varid_1, "_fillValue", 0.0

NCDF_ATTPUT, id, varid_2, "title", "pressure"
NCDF_ATTPUT, id, varid_2, "units", "Pa"
NCDF_ATTPUT, id, varid_2, "_fillValue",0.0

NCDF_ATTPUT, id, varid_2, "title", "zenith angle"
NCDF_ATTPUT, id, varid_3, "units", "degrees"
NCDF_ATTPUT, id, varid_3, "_fillValue",0.0

NCDF_CONTROL, id, /ENDEF

;if (do_zenith) then begin
;  NCDF_VARPUT,id,varid_1, prod_zen
;endif else begin
  NCDF_VARPUT,id,varid_1, prod
;endelse

NCDF_VARPUT,id,varid_2, pressure(*,0)
NCDF_VARPUT,id,varid_3, cosz

NCDF_CLOSE, id


print, ""
print, "______________________DONE__________________________________"
print, ""

end
