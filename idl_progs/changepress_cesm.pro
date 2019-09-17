pro changepress_cesm
;----------------------------------------------
;AUTHOR: WOLF, E.T.
;9/3/2013
;-----------------------------------------------
; NOTES:
; note do_dry option which sets initial Q, CLDICE, and CLDLIQ to zero
;
; Change pressure for initial conditions files in CESM
; Known issue:  In CESM's coordinate system, the layers below
; 100 mb are dialted or contracted to follow the terrain (i.e. surface
; pressure).  When changing the total pressure to much less than or
; greater than 1 bar, the bottom model layers may become sub-optimal
; and cause problems. 
;
; at high-P the layers become
; at low-P the layers bunch together and can cross, causing crash
;
; Future work should actualy adjust/create new hybrid sigma
; coefficient when changing pressure.
;-----------------------------------------------


do_write = 1  ;if 1 write output file
do_dry = 1

mwn2 = 28.
mwar = 40.
mwco2 = 44.
mwch4 = 16.

cpn2 = 1.039e3
cpar = 0.520e3
cpco2 = 0.846e3
cpch4 = 2.226e3

co2bar = 0.0 ;bar
ch4bar = 0.0 ;bar
n2bar = 10.0  ; - co2bar - ch4bar      ;bar

psbar = co2bar + ch4bar  + n2bar

co2vmr = co2bar / psbar ;vmr
ch4vmr = ch4bar / psbar ;vmr
n2vmr = n2bar / psbar   ;vmr

mwdry = co2vmr*mwco2 + ch4vmr*mwch4 + n2vmr*mwn2 
cpdry = co2vmr*cpco2 + ch4vmr*cpch4 + n2vmr*cpn2 


print, "============================="
print, "N2BAR: ", N2BAR
print, "CO2BAR: ", CO2BAR
print, "CH4BAR: ", CH4BAR
print, "------------------------------"
print, "TOTAL: ",N2BAR+CO2BAR+CH4BAR
print, "============================="
print, "N2VMR: ", N2VMR
print, "CO2VMR: ", CO2VMR
print, "CH4VMR: ", CH4VMR
print, "------------------------------"
print, "TOTAL: ",N2VMR+CO2VMR+CH4VMR
print, "------------------------------"
print, "mwdry: ", mwdry
print, "cpdry: ", cpdry
print, "------------------------------"


;==============================================================
;  file_out =  '/projects/btoon/wolfet/exofiles/atm/CO2_0.00000976562bar_L45_ic.nc'
;  file_out =   '/projects/btoon/wolfet/exofiles/atm/ic_1barN2_0.2barCO2_L40_ic.nc'
  file_out =  '/gpfsm/dnb53/etwolf/models/ExoCAM/cesm1.2.1/initial_files/cam4_aqua_fv/ic_10bar_L51_zmean_dry_ic.nc'
;  file_out = '/projects/btoon/wolfet/exofiles/atm/ohz_t2700K_temp_ic.nc'
;  file_out = '/projects/btoon/wolfet/exofiles/atm/archean_360ppmCO2_1bN2_L40_ic.nc'
;  ncdata_in = '/projects/btoon/wolfet/exofiles/atm/cami_0001-01-01_4x5_L26_c060608.nc'

; standard present day initial conditions

 ncdata_in =  '/gpfsm/dnb53/etwolf/models/ExoCAM/cesm1.2.1/initial_files/cam4_aqua_fv/ic_1bar_L51_zmean_ic.nc'
;ncdata_in = '/projects/wolfet/EXO_RESTART/trappist1e_1barN2_0.1barCO2_aqua_0061-01-01/trappist1e_1barN2_0.1barCO2_aqua.cam.i.0061-01-01-00000.nc'
;ncdata_in = '/gpfs/summit/scratch/wolfet/archive/t3300_s1400_p22.1392_0.25bar/rest/0017-01-01-00000/t3300_s1400_p22.1392_0.25bar.cam.i.0017-01-01-00000.nc'
  ;ncdata_in = '/projects/btoon/wolfet/exofiles/atm/ic_1bar_L51_ic.nc'
;  ncdata_in = '/projects/btoon/wolfet/exofiles/atm/control_L48.cam.i.0048-01-01-00000.nc'

   ;standard ohz start
;           ncdata_in = '/home/wolfet/ExoCAM_cesm1.2.1/initial_files/cam4_fv_4x5_aquaplanet/control_L40.cam.i.0009-01-01-00000.nc'
        ;  ncdata_in = '/lustre/janus_scratch/wolfet/t2700K_s304.836_p14.8548_ohz/run/t2700K_s304.836_p14.8548_ohz.cam.i.0002-01-01-00000.nc'
;hot climate
;ncdata_in = '/lustre/janus_scratch/wolfet/archive/CO2_184320ppm/rest/0049-01-01-00000/CO2_184320ppm.cam.i.0049-01-01-00000.nc'
  outstring = "cp " + ncdata_in + " "+ file_out


if (do_write eq 1) then  spawn, outstring

;------load metrics--------------------------
  ncid=ncdf_open(ncdata_in, /nowrite)
  ncdf_varget,ncid,'PS',PS_in              ;read surface pressure array, PS[lon,lat]
  ncdf_varget,ncid,'P0',P0_in              ;read reference surface pressure
  ncdf_varget,ncid,'T',T_in              ;read reference surface pressure
;  ncdf_varget,ncid,'TS',TS_in              ;read reference surface pressure

  if (do_dry eq 1) then begin
    ncdf_varget,ncid,'Q',Q   & Q(*,*,*) = 0.0
    ncdf_varget,ncid,'CLDLIQ',CLDLIQ  & CLDLIQ(*,*,*) = 0.0
    ncdf_varget,ncid,'CLDICE',CLDICE  & CLDICE(*,*,*) = 0.0
  endif

;  if (do_iso eq 1) then begin   
;  endif
  ncdf_varget,ncid,'hyai',hyai
  ncdf_varget,ncid,'hybi',hybi
  ncdf_varget,ncid,'hyam',hyam
  ncdf_varget,ncid,'hybm',hybm

  ncdf_varget,ncid,'lon',lon
  ncdf_varget,ncid,'lat',lat
  ncdf_varget,ncid,'lev',lev
  ncdf_close,ncid

  nlon=n_elements(lon) 
  nlat=n_elements(lat) 
  nlev=n_elements(lev) 
  nilev = nlev + 1
  PS=fltarr(nlon,nlat)

;  area_weighted_avg_gen, lon, lat, PS_in, PS_avg
;  print, PS_avg

  for x=0,nlon-1 do begin
    for y=0,nlat-1 do begin
      for z=0, nlev-1 do begin
      endfor
;      PS(x,y) = PS_in(x,y)*psbar*1.0e5/P0_in
       PS(x,y) = psbar*1.0e5
;      PS(x,y) = PS_in(x,y)*psbar*1.0e5/PS_avg
    endfor
  endfor
  P0 = psbar*1.0e5


  ;; testing ---
  lev_P=fltarr(nlon,nlat,nlev)
  ilev_P=fltarr(nlon,nlat,nilev)
  hybrid2pressure,nlon,nlat,nlev,PS,P0,hyam,hybm,hyai,hybi,lev_P,ilev_P
  lev_P(*,*,*) = lev_P(*,*,*)/100.0   ; & lev_Pmulti(index, *,*,*) = lev_P(*,*,*)
  ilev_P(*,*,*) = ilev_P(*,*,*)/100.0  ;& ilev_Pmulti(index, *,*,*) = ilev_P(*,*,*)
  ;; testing ---


if (do_write eq 1) then begin
  ;--------- update ic file ----------
  ncid = ncdf_open(file_out, /WRITE)
  print, "updating atmosphere initial conditions file"
  print, file_out

  ;straight scaling
  ncdf_varput, ncid, 'P0',psbar*1.0e5
  ncdf_varput, ncid, 'PS',PS
  if (do_dry eq 1) then begin
    ncdf_varput, ncid, 'Q',Q
    ncdf_varput, ncid, 'CLDICE',CLDICE
    ncdf_varput, ncid, 'CLDLIQ',CLDLIQ
  endif

  ncdf_close, ncid
endif

end
