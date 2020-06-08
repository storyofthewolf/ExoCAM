pro carma_optics
;---------------------------------------------------------------
;AUTHOR: WOLF, E.T.
;4/06/09  
;
;PURPOSE: Create initial condition file for haze optical 
; properties, namelist variable (carma_optics_file)
;
;METHOD: 
;
;NOTES:  add functionality for fractal particles
;        make sure to check haze bin parameters are correct.
;---------------------------------------------------------------



print, "begin carmaoptics . . . "

;-------------------parameters-----------------------------------
outputfile="khare.n40.nc"      ;name of CARMA optics file
wvl_cam = 19                 ;number of cam spectral intervals
nelem = 1                    ;number of CARMA elements
nbin_param = 40                    ;number of CARMA size bins
rmrat_param = 2.5            ;ratio between masses of sucessive bins
;------------------some constants-------------------------------
c = 2.998E8     ;[m s-1]
h = 6.626E-34 ;[J s-1]
kb = 1.38E-23   ;[J K-1]

;----------------Read titan optical constants--------------------
inputfile='/Users/wolf/IDLWorkspace/titan.optical.constants.txt'
rows = 90
num_wvlopt = rows
wvlopt = fltarr(rows)     ;[um]
a = fltarr(rows)
b = fltarr(rows)
n = fltarr(rows)       ; real part of index of refraction
k = fltarr(rows)       ; imaginary part of the index of refraction, computed from a,b
data = fltarr(4,rows)
header = strarr(4)
  
OPENR,lun,inputfile,/GET_LUN

READF,lun,header
READF,lun, data

wvlopt(*) = data(0,*)
a(*) = data(1,*)
b(*) = data(2,*)
n(*) = data(3,*)

;compute imaginary index of refraction
k(*) = a(*)*10.0^(-b(*))

;reverse order of arrays
wvlopt = reverse(wvlopt)
n = reverse(n)
k = reverse(k)

;vectors need for future use
;wvlopt, n,k = this will be interpolated to better smaple the
;wavelenth bins

FREE_LUN,lun

;----------------------Read CAM WAV parameters---------------------
cam_wavfile = '/Users/wolf/IDLWorkspace/cam.wavbins.txt'
rows = 19
cam_wvlmin = fltarr(rows)
cam_wvlmax = fltarr(rows)
cam_wvlmid = fltarr(rows)
datawav = fltarr(3,rows)
headerwav = strarr(2)

OPENR,lun1,cam_wavfile,/GET_LUN

READF,lun1,headerwav
READF,lun1,datawav

cam_wvlmin(*) = datawav(1,*)
cam_wvlmax(*) = datawav(2,*)
cam_wvlmid(*) = (cam_wvlmin(*)+cam_wvlmax(*))/2.0

FREE_LUN,lun1

;------------- interpolate optical properties to size bins-----------
;NOTES: Khare optical constants are interpolated to better sample
;the carma size bins.  This will allow a more accurate Chandrasekhar 
;flux weighting 

rows = 1004
int_num_wvlopt = rows
;interpolate arrays of khare optical constants

int_wvlopt = fltarr(rows)
int_n = fltarr(rows)
int_k = fltarr(rows)

int_wvlopt(*) = 0.0
int_n(*) = 0.0
int_k(*) = 0.0

;CAM sw wavelength bins stretch from 0.200 to 5.000 um
delta_wvl = 5.0-0.2
grid_space=delta_wvl/(rows-5)

int_wvlopt(0) = cam_wvlmin(0)-2.0*grid_space
FOR l = 1, int_num_wvlopt-1 DO BEGIN
  int_wvlopt(l)=int_wvlopt(l-1)+grid_space 
ENDFOR


;interpolate n,k to fit new grid
;NOTES: I should be safe with error check to prevent array seg fault
;since I ensure the range of int_num_wvlopt is sufficiently large.
;The range of khare optical constants also exceeds need range.

FOR u=2, int_num_wvlopt-1 DO BEGIN
  FOR v=0, num_wvlopt-2 DO BEGIN
    ;find encompassing data

    if((wvlopt(v) lt int_wvlopt(u)) && (wvlopt(v+1) gt int_wvlopt(u))) then begin

      ;interpolate n
      I_X = [wvlopt(v-1), wvlopt(v),wvlopt(v+1),wvlopt(v+2)]
      I_Y = [n(v-1),n(v),n(v+1),n(v+2)]
      I_U = [int_wvlopt(u)]
      I_INT = interpol(I_Y,I_X,I_U)
      int_n(u) = I_INT

      ;interpolate k
      I_X = [wvlopt(v-1), wvlopt(v),wvlopt(v+1),wvlopt(v+2)]
      I_Y = [k(v-1),k(v),k(v+1),k(v+2)]
      I_U = [int_wvlopt(u)]
      I_INT = interpol(I_Y,I_X,I_U)
      int_k(u) = I_INT
 
    endif
  ENDFOR
ENDFOR

;-----------------------haze bin information-----------------------
;haze particle parametes
rho_haze = 0.64       ; [g cm-3] haze particle density        
nbin = nbin_param
rmin = 1.0E-7         ; [cm] 
rmrat = rmrat_param
rmassmin = 4.0/3.0*!pi*(rmin^3.0)*rho_haze    ;[g]

dr = fltarr(nbin)    ;[um] vector containing the width of each bin in um
bin_radius = fltarr(nbin)  ;[cm]
bin_mass = fltarr(nbin)    ;[g]
vrfact = ((3.0/2.0/!pi/(rmrat+1.0))^(1.0/3.0))*(rmrat^(1.0/3.0)-1.0)

FOR bin_idx=0, nbin-1 DO BEGIN
  bin_mass(bin_idx)=rmassmin*rmrat^(bin_idx)
  bin_radius(bin_idx)=((3.0*bin_mass(bin_idx))/(4.0*!pi*rho_haze))^(1.0/3)     ;[cm]
  dr(bin_idx)=vrfact*(bin_mass(bin_idx)/rho_haze)^(1.0/3.0)*(1.0E4) ;[um] bin width, 1.0e4 converts cm to um 
ENDFOR

bin_radius(*)=bin_radius(*)*1.0E4   ; [um]

;-------------------------print parameters---------------------------
print, "----------------------------------"
print, "filename = ",outputfile
print, "wvl_cam = ",wvl_cam
print, "nelem = ",nelem
print, "nbin = ",nbin
print, "rmrat = ",rmrat
print, "inputfile = ",inputfile
print, "----------------------------------"


;------------------------computation of optical properties---------------
kcarma_dat = fltarr(nelem,nbin,wvl_cam)
wcarma_dat = fltarr(nelem,nbin,wvl_cam)
gcarma_dat = fltarr(nelem,nbin,wvl_cam)

FOR elem_idx = 0, nelem-1 DO BEGIN
 
  FOR bin_idx = 0, nbin-1 DO BEGIN   ;  r bin loop

    FOR wvlbin_idx = 0, wvl_cam-1 DO BEGIN  ; wavlength bin loop

      Kext = 0.0
      W = 0.0
      G = 0.0
      planck_tot = 0.0

      FOR wvlopt_idx = 0, int_num_wvlopt-1 DO BEGIN  ;loop over all interpolated fine grid points

        if ((int_wvlopt(wvlopt_idx) GE cam_wvlmin(wvlbin_idx)) && (int_wvlopt(wvlopt_idx) LE  cam_wvlmax(wvlbin_idx))) then begin

         
          cm = complex(int_n(wvlopt_idx),-int_k(wvlopt_idx))                   ; set complex indices of refraction

          ;calculate optical parameters as average of bin end points and center
          ;should i sample more heavily for averaging?
                            
          alpha_left = 2.0 * !pi * (bin_radius(bin_idx) - dr(bin_idx) / 2.0) / int_wvlopt(wvlopt_idx) 
          alpha_center = 2.0 * !pi * bin_radius(bin_idx) / int_wvlopt(wvlopt_idx)  
          alpha_right = 2.0 * !pi * (bin_radius(bin_idx) + dr(bin_idx) / 2.0) / int_wvlopt(wvlopt_idx)     
          
          mie_single, alpha_left,cm,dqext,dqscat,dqbk,dg
          Qext_left = dqext
          Qscat_left = dqscat
          G_left = dg
          W_left = dqscat/dqext  

          mie_single, alpha_center,cm,dqext,dqscat,dqbk,dg
          Qext_center = dqext
          Qscat_center = dqscat
          G_center = dg
          W_center = dqscat/dqext  

          mie_single, alpha_right,cm,dqext,dqscat,dqbk,dg
          Qext_right = dqext
          Qscat_right = dqscat
          G_right = dg
          W_right = dqscat/dqext  

          Qext = (Qext_left + Qext_center + Qext_right) / 3.0
          Qscat = (Qscat_left + Qscat_center + Qscat_right) / 3.0
          Gavg = (G_left + G_center + G_right) / 3.0
          
          delta_wvl = ((int_wvlopt(wvlopt_idx+1) - int_wvlopt(wvlopt_idx-1))/2.0)*1.0E-6
          
          ;FIX HERE!!!!  EW - 7/18/2011
          ;Below I had implemented a planck weighting of the optical
          ;constants within each spectral interval.  This is incorrect. For shortwave need to use
          ;weight via the solar insolation at wach wavelength.
          ;planck = 2*h*c^2.0/((int_wvlopt(wvlopt_idx)*1.0E-6)^5.0*exp(h*c/kb/(int_wvlopt(wvlopt_idx)*1.0E-6)/5777.0)-1)
          ;planck_tot = planck_tot + planck * delta_wvl       

          ;Should be weighted sum of Kext * solarflux in  each fine subinterval
          Kext = Kext + 3.0*Qext/4.0/rho_haze/(bin_radius(bin_idx)*1.0E-4) *1.0E-4 * delta_wvl * planck
          W = W + Qscat/Qext * planck*delta_wvl
          G = G + Gavg * planck*delta_wvl

        endif 
         
      ENDFOR

      ;FIX: divide out by total "solar" irradiance in the spectral interval
      kcarma_dat(elem_idx,bin_idx,wvlbin_idx) = Kext/planck_tot
      wcarma_dat(elem_idx,bin_idx,wvlbin_idx) = W/planck_tot
      gcarma_dat(elem_idx,bin_idx,wvlbin_idx) = G/planck_tot
    ENDFOR    ; end cam wav loop   

   ENDFOR     ; end r bin loop

ENDFOR

;--------------netcdf creation-----------------------------
print, "building netcdf file . . . "

id = NCDF_CREATE(outputfile, /CLOBBER)
wvl_cam_id = NCDF_DIMDEF(id,'wvl_cam',wvl_cam)
nelem_id = NCDF_DIMDEF(id,'elements',nelem)
nbin_id = NCDF_DIMDEF(id,'bins',nbin)

kcarma_id = NCDF_VARDEF(id, 'ext_carma',[nelem_id,nbin_id,wvl_cam_id])
wcarma_id = NCDF_VARDEF(id, 'ssa_carma',[nelem_id,nbin_id,wvl_cam_id])
gcarma_id = NCDF_VARDEF(id, 'asm_carma',[nelem_id,nbin_id,wvl_cam_id])

NCDF_ATTPUT, id, kcarma_id, "title", "specific extinction - early earth haze"
NCDF_ATTPUT, id, kcarma_id, "units", "m2 g-1"
NCDF_ATTPUT, id, wcarma_id, "title", "single scattering albedo - early earth haze"
NCDF_ATTPUT, id, wcarma_id, "units", "fraction"
NCDF_ATTPUT, id, gcarma_id, "title", "asymmetry parameter - early earth haze"
NCDF_ATTPUT, id, gcarma_id, "units", "fraction"
NCDF_CONTROL, id, /ENDEF

NCDF_VARPUT, id, kcarma_id, kcarma_dat
NCDF_VARPUT, id, wcarma_id, wcarma_dat
NCDF_VARPUT, id, gcarma_id, gcarma_dat

NCDF_CLOSE, id
end
