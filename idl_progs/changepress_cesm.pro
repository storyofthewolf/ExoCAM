pro changepress_cesm
;-----------------------------------------------------------------------
;AUTHOR: WOLF, E.T.
; 9/3/2013
; revised and simplified 10/21/2021
;-----------------------------------------------------------------------
; Description:  Changes mean surface temperature in atmosphere initial
; condition file (ncdata). In ExoCAM, the atmospheric composition is set in
; exoplanet_mod.F90 by setting the partial pressures of consituent gases.  The
; ncdata data file must then be changed to match.  Note, that small
; (~percentage) deviations between exoplanet_mod.F90 total pressure and PS in
; the initial condition file is ok, and the value in exoplanet_mod.F90 will be
; used in the run.  However larger (~factors) deviations between the two will
; result in model crashes, and thus new ics will need to be made.
;
; do_dry option which sets initial Q, CLDICE, and CLDLIQ to zero
;
; NOTES:  presently this doesn't work for the cubed sphere dynamical core, 
;         must generalize.
;
;-----------------------------------------------------------------------


do_write = 1  ; if 1 write output file
do_dry = 1    ; remove water vapor

;choose one or the other
do_flat_psfield = 0 ; write output pressures equal to sum of partial pressures
                                ; pressure field is assumed flat across
                                ; surface, regardless of topography

do_scale_psfield = 1 ; write output pressures by scaling press field from 
                                ; the input file in being used (ncdata).

if (do_flat_psfield eq 1 and do_scale_psfield eq 1) then begin
  print, " cannot have both do_scale_psfield and do_flat_psfield equal 1 "
  print, " ending without calculation "
  GOTO, FINISH
endif

if (do_scale_psfield eq 1) then print, "scaling input surface pressure field"
if (do_flat_psfield eq 1)  then print, "applying flat surface pressure field"

;++++++++++++++++++++++++++++++++++++++++++++++++++++++
;------------------------------------------------------
;---  Choose your new mean surface pressure in BARS ---
new_PS_mean = 7.0   ; [bars]
;------------------------------------------------------
;++++++++++++++++++++++++++++++++++++++++++++++++++++++

;=======================================================================================================================
;=======================================================================================================================
;  FILES
;=======================================================================================================================
;=======================================================================================================================
;choose one from each

;---  assorted template initial files ---
;ncdata_in = '/discover/nobackup/etwolf/models/CESM_Mars/marsfiles/atm/mars.cam.i.0002-01-01-00000.nc'
;ncdata_in = '/gpfsm/dnb53/etwolf/models/ExoCAM/cesm1.2.1/initial_files/cam_aqua_fv/ic_1bar_L51_300Kiso_Q0.01_ic.nc'
;ncdata_in = '/gpfsm/dnb53/etwolf/models/ExoCAM/cesm1.2.1/initial_files/cam_aqua_fv/ic_P0.25bar_L40_ic.nc'
;ncdata_in = '/gpfsm/dnb53/etwolf/cesm_scratch/rundir/mars_2barCO2/run/mars_2barCO2.cam.i.0012-02-01-00000.nc'
;ncdata_in = '/gpfsm/dnb53/etwolf/models/ExoCAM/cesm1.2.1/initial_files/cam_aqua_fv/ic_P0.1bar_L40_dry_ic.nc'
;ncdata_in = '/gpfsm/dnb53/etwolf/cesm_scratch/archive/mars_0.1barCO2/rest/0031-01-01-00000/mars_0.1barCO2.cam.i.0031-01-01-00000.nc'
;ncdata_in = '/gpfsm/dnb53/etwolf/cesm_scratch/archive/mars_dev2/rest/0001-01-21-00000/mars_dev2.cam.i.0001-01-21-00000.nc'
;ncdata_in = '/gpfsm/dnb53/etwolf/models/ExoCAM/cesm1.2.1/initial_files/cam4_aqua_se/cami_0000-01-01_ne16np4_L48_c120525_aquaplanet.nc'
ncdata_in = '/gpfsm/dnb53/etwolf/models/ExoCAM/cesm1.2.1/initial_files/mars/atm/paleo_mars_6.0bar_dry.cam.i.0021-01-01-00000.nc'

;--- assorted output file names ---
;file_out = 'ic_100bar_L51_iso300_ic.nc'   ; writes to local directory
;file_out =  '/gpfsm/dnb53/etwolf/models/ExoCAM/cesm1.2.1/initial_files/cam_land_fv/ic_P1.5bar_L51_ic.nc'
;file_out =  '/gpfsm/dnb53/etwolf/models/ExoCAM/cesm1.2.1/initial_files/cam_aqua_fv/ic_P0.01bar_L40_ic.nc'
;file_out =  '/discover/nobackup/etwolf/models/CESM_Mars/marsfiles/atm/mars_4bar.cami.4x5.dry.nc'
;file_out =  '/gpfsm/dnb53/etwolf/models/ExoCAM/cesm1.2.1/initial_files/cam_aqua_fv/ic_P8bar_L40_dry_ic.nc'
;file_out =  '/gpfsm/dnb53/etwolf/models/CESM_Mars/marsfiles/atm/mars_dev2_0.01bar.cam.i.0001-01-21-00000.nc'
;file_out = '/gpfsm/dnb53/etwolf/models/ExoCAM/cesm1.2.1/initial_files/cam_aqua_se/cami_0000-01-01_ne16np4_L48_10bar_c120525_aquaplanet.nc'
file_out = '/gpfsm/dnb53/etwolf/models/ExoCAM/cesm1.2.1/initial_files/mars/atm/paleo_mars_7.0bar_dry.cam.i.0021-01-01-00000.nc'

;=======================================================================================================================
;=======================================================================================================================
;=======================================================================================================================

; copy template file to new file locations
outstring = "cp " + ncdata_in + " "+ file_out
if (do_write eq 1) then  spawn, outstring

;------load metrics--------------------------
  ncid=ncdf_open(ncdata_in, /nowrite)
  ncdf_varget,ncid,'PS',PS_in              ;read surface pressure array, PS[lon,lat]
  ncdf_varget,ncid,'P0',P0_in              ;read reference surface pressure
  ncdf_varget,ncid,'T',T_in              ;read reference surface pressure
;  ncdf_varget,ncid,'TS',TS_in              ;read reference surface pressure

  if (do_dry eq 1) then begin
    ; set water vapor and cloud fields to zero
    ncdf_varget,ncid,'Q',Q   & Q(*,*,*) = 0.0
    ncdf_varget,ncid,'CLDLIQ',CLDLIQ  & CLDLIQ(*,*,*) = 0.0
    ncdf_varget,ncid,'CLDICE',CLDICE  & CLDICE(*,*,*) = 0.0
  endif

  ; read hybrid sigma levels
  ncdf_varget,ncid,'hyai',hyai
  ncdf_varget,ncid,'hybi',hybi
  ncdf_varget,ncid,'hyam',hyam
  ncdf_varget,ncid,'hybm',hybm

  ; read lon, lat, lev fields
  ncdf_varget,ncid,'lon',lon
  ncdf_varget,ncid,'lat',lat
  ncdf_varget,ncid,'lev',lev
  ncdf_close,ncid

  nlon=n_elements(lon) 
  nlat=n_elements(lat) 
  nlev=n_elements(lev) 
  nilev = nlev + 1
  PS=fltarr(nlon,nlat)

  area_weighted_avg_gen, lon, lat, PS_in, PS_avg
  print, "initial file pressure: ", PS_avg, P0_in
  scalefac = new_PS_mean/(P0_in/1.0e5)
  print, "scalefac: ", scalefac

  for x=0,nlon-1 do begin
    for y=0,nlat-1 do begin
      if (do_flat_psfield) then  PS(x,y) = new_PS_mean*1.0e5
      if (do_scale_psfield) then  PS(x,y) = PS_in(x,y)*scalefac
    endfor
  endfor
  if (do_flat_psfield) then P0_out = new_PS_mean*1.0e5
  if (do_scale_psfield) then P0_out = P0_in*scalefac

  print, "P0 new", P0_out

if (do_write eq 1) then begin
  ;--------- update ic file ----------
  ncid = ncdf_open(file_out, /WRITE)
  print, "updating atmosphere initial conditions file"
  print, file_out

  ;straight scaling
  ncdf_varput, ncid, 'P0',P0_out
  ncdf_varput, ncid, 'PS',PS
  if (do_dry eq 1) then begin
    ncdf_varput, ncid, 'Q',Q
    ncdf_varput, ncid, 'CLDICE',CLDICE
    ncdf_varput, ncid, 'CLDLIQ',CLDLIQ
  endif
  ncdf_close, ncid
endif

FINISH:
end
