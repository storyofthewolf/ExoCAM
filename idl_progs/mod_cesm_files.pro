pro mod_cesm_files
;-------------------------------
; IDL script used to modify initial conditions files
; This code is a mess, use at your own peril, support is limited.
; -- e t w --
;-------------------------------

make_bndtopo = 0 ; topograph files  
make_popfrc = 0 ; modify existing popfrc file
make_ncdata = 0
make_tropopause_climo_file = 0
make_ozone = 0
make_aerosoldep = 0
make_prescribed_aerofil = 0
make_lnd_domain = 0
make_ocn_domain = 0
make_ocn_master = 1   ; make new popfrc files, because one doesn't exist on file

make_ne5np4_popfile = 0
make_ne5np4_icfile = 0

make_ne16np4_popfile = 0
make_ne16np4_icfile = 0
make_n16np4_landdomian = 0
make_n16np4_surdat = 0



;; spectral element ne16np4_icfile
;file_ne16np4_icfile = '/lustre/janus_scratch/cesm/inputdata/atm/cam/inic/homme/cami_0000-01-01_ne16np4_L26_c120525.nc'
;file_ne16np4_icfile_out = '/projects/btoon/wolfet/exofiles/ocn/aquaplanet/cami_0000-01-01_ne16np4_L26_c120525_aquaplanet.nc'
;file_ne16np4_icfile_out = '/projects/btoon/wolfet/exofiles/ocn/aquaplanet/cami_0000-01-01_ne16np4_L26_c120525_aquaplanet_ice.nc'
;file_bndtopo = "/projects/btoon/wolfet/exofiles/ocn/aquaplanet/USGS-gtopo30_1.9x2.5_smooth500-50_ne16np4_c050602.nc"
;file_bndtopo_out = '/projects/btoon/wolfet/exofiles/ocn/aquaplanet/USGS-gtopo30_1.9x2.5_smooth500-50_ne16np4_c050602_aquaplanet.nc'
;file_dom_domain = '/gpfs/summit/datasets/CESM/inputdata/share/domains/domain.ocn.ne16np4_gx3v7.120406.nc'
;file_dom_domain_out = '/projects/btoon/wolfet/exofiles/ocn/aquaplanet/domain.ocn.ne16np4_gx3v7.120406_aquaplanet.nc'


; spectral element ne5np4_icfile
;file_ne5np4_icfile = '/projects/btoon/wolfet/exofiles/ocn/ne5np4/cami_0000-01-01_ne5np4_L30_ape_c101201.nc'
;file_ne5np4_icfile_out = '/projects/btoon/wolfet/exofiles/ocn/ne5np4/cami_0000-01-01_ne5np4_L30_ape_c101201_aquaplanet.nc'
;file_bndtopo = "/projects/btoon/wolfet/exofiles/ocn/ne5np4/USGS-gtopo30_ne5np4_c140808.nc"
;file_bndtopo_out = "/projects/btoon/wolfet/exofiles/ocn/ne5np4/USGS-gtopo30_ne5np4_c140808_aquaplanet.nc"
;file_dom_domain = "/projects/btoon/wolfet/exofiles/ocn/ne5np4/domain.ocn.ne5np4_gx3v7.140810.nc"
;file_dom_domain_out = "/projects/btoon/wolfet/exofiles/ocn/ne5np4/domain.ocn.ne5np4_gx3v7.140810_aquaplanet.nc"




file_tropopause_climo = '/lustre/janus_scratch/cesm/inputdata/atm/cam/chem/trop_mozart/ub/clim_p_trop.nc'
file_tropopause_climo_out = '/projects/btoon/wolfet/exofiles/ocn/aquaplanet/clim_p_trop_aquaplanet.nc'
file_ozone = '/lustre/janus_scratch/cesm/inputdata/atm/cam/ozone/ozone_1.9x2.5_L26_2000clim_c091112.nc'
file_ozone_out = '/projects/btoon/wolfet/exofiles/ocn/aquaplanet/ozone_1.9x2.5_L26_2000clim_c091112_aquaplanet_ZEROD.nc'
file_prescribed_aero = '/lustre/janus_scratch/cesm/inputdata/atm/cam/chem/trop_mozart_aero/aero/aero_1.9x2.5_L26_2000clim_c091112.nc'
file_prescribed_aero_out = '/projects/btoon/wolfet/exofiles/ocn/aquaplanet/aero_1.9x2.5_L26_2000clim_c091112_ZEROD.nc'
file_aerosoldep = '/lustre/janus_scratch/cesm/inputdata/atm/cam/chem/trop_mozart_aero/aero/aerosoldep_monthly_1849-2006_1.9x2.5_c090803.nc'
file_aerosoldep_out = '/projects/btoon/wolfet/exofiles/ocn/aquaplanet/aerosoldep_monthly_1849-2006_1.9x2.5_c090803_ZEROD.nc'


;files lists 
 ;--- 1.9x2.5 ---
;file_bndtopo = '/lustr0e/janus_scratch/cesm/inputdata/atm/cam/topo/USGS-gtopo30_1.9x2.5_remap_c050602.nc'
;file_dom_domain = '/lustre/janus_scratch/cesm/inputdata/atm/cam/ocnfrac/domain.camocn.1.9x2.5_gx1v6_090403.nc'
;file_domfile = '/lustre/janus_scratch/cesm/inputdata/ocn/docn7/SOM/pop_frc.1x1d.090130.nc'
;file_ncdata = '/lustre/janus_scratch/cesm/inputdata/atm/cam/inic/fv/cami_0000-01-01_1.9x2.5_L26_c070408.nc'
;file_tropopause_climo = '/lustre/janus_scratch/cesm/inputdata/atm/cam/chem/trop_mozart/ub/clim_p_trop.nc'
;file_ozone = '/lustre/janus_scratch/cesm/inputdata/atm/cam/ozone/ozone_1.9x2.5_L26_2000clim_c090803.nc'
;file_prescribaero = '/lustre/janus_scratch/cesm/inputdata/atm/cam/chem/trop_mozart_aero/aero/aero_1.9x2.5_L26_2000clim_c090803.nc'
;file_aerosoldep = '/lustre/janus_scratch/cesm/inputdata/atm/cam/chem/trop_mozart_aero/aero/aerosoldep_monthly_1849-2006_1.9x2.5_c090803.nc'

;--- 4x5 ---  from my simulations
;file_bndtopo = '/lustre/janus_scratch/cesm/inputdata/atm/cam/topo/USGS-gtopo30_4x5_remap_c050520.nc'
;file_dom_domain = '/lustre/janus_scratch/cesm/inputdata/share/domains/domain.ocn.gx3v7.120323.nc'
;file_domfile = '/lustre/janus_scratch/cesm/inputdata/ocn/docn7/SOM/pop_frc.gx3v7.110128.nc'
;;;;file_ncdata = '/lustre/janus_scratch/cesm/inputdata/atm/cam/inic/fv/cami_0001-01-01_4x5_L26_c060608.nc'
;;;file_ncdata ='/projects/wolfet/EXO_RESTART/control_L26_0053-01-01-00000/control.cam.i.0053-01-01-00000.nc'
;;;file_ncdata = '/projects/wolfet/EXO_RESTART/control_L45_0040-01-01-00000/control_L45.cam.i.0040-01-01-00000.nc'

;file_tropopause_climo = '/lustre/janus_scratch/cesm/inputdata/atm/cam/chem/trop_mozart/ub/clim_p_trop.nc'
;file_ozone = '/lustre/janus_scratch/cesm/inputdata/atm/cam/ozone/ozone_1.9x2.5_L26_2000clim_c091112.nc'
;file_prescribaero = '/lustre/janus_scratch/cesm/inputdata/atm/cam/chem/trop_mozart_aero/aero/aero_1.9x2.5_L26_2000clim_c091112.nc'
;file_aerosoldep = '/lustre/janus_scratch/cesm/inputdata/atm/cam/chem/trop_mozart_aero/aero/aerosoldep_monthly_1849-2006_1.9x2.5_c090803.nc'
;file_lnd_domain = '/gpfs/summit/datasets/CESM/inputdata/share/domains/domain.lnd.fv4x5_gx3v7.091218.nc'
;file_lnd_domain = '/gpfs/summit/datasets/CESM/inputdata/share/domains/domain.lnd.fv4x5_gx3v7.091218.nc'
;file_ocn_domain = '/gpfs/summit/datasets/CESM/inputdata/share/domains/domain.ocn.0.47x0.63_gx1v6_090408.nc'
;file_ocn_master = '/lustre/janus_scratch/cesm/inputdata/ocn/docn7/SOM/pop_frc.1x1d.090130.nc'


;---  topography files in/out ---
;file_bndtopo = "/gpfs/summit/datasets/CESM/inputdata/atm/cam/topo/USGS_gtopo30_0.47x0.63_remap_c061106.nc"
;file_bndtopo_out = "/projects/wolfet/models/ExoCAM/cesm1.2.1/initial_files/cam4_aqua_fv/topo_0.47x0.63_aquaplanet.nc"

;--- ocean domain files in/out ---
;file_ocn_domain = '/gpfs/summit/datasets/CESM/inputdata/share/domains/domain.ocn.0.47x0.63_gx1v6_090408.nc'
;file_ocn_domain_out =  "/projects/wolfet/models/ExoCAM/cesm1.2.1/initial_files/cam4_aqua_fv/domain.ocn.0.47x0.63_aquaplanet.nc"

;--- ocn master files in/out ---
;file_ocn_master = "/gpfs/summit/datasets/CESM/inputdata/share/domains/domain.ocn.0.47x0.63_gx1v6_090408.nc"
;file_ocn_master_out = "/projects/wolfet/models/ExoCAM/cesm1.2.1/initial_files/cam4_aqua_fv/pop_frc.0.47x0.63d_aquaplanet_0OHT.nc"


;--- pop.frc files in/out ---
file_pop_frc = "/gpfsm/dnb53/etwolf/models/ExoCAM/cesm1.2.1/initial_files/cam4_aqua_fv/pop_frc.4x5d.090130_aquaplanet_300Kiso.nc"
file_pop_frc_out = "/gpfsm/dnb53/etwolf/models/ExoCAM/cesm1.2.1/initial_files/cam4_aqua_fv/pop_frc.4x5d.090130_aquaplanet_200Kiso.nc"


;--- ncdata files in/oout --
;file_ncdata = "/gpfsm/dnb53/etwolf/models/ExoCAM/cesm1.2.1/initial_files/cam4_aqua_fv/ic_P1bar_L40_300Kiso_ic.nc"
file_ncdata = "/gpfsm/dnb53/etwolf/models/ExoCAM/cesm1.2.1/initial_files/cam4_aqua_fv/ic_P1bar_L40_ic.nc"
file_ncdata_out = "/gpfsm/dnb53/etwolf/models/ExoCAM/cesm1.2.1/initial_files/cam4_aqua_fv/ic_P1bar_L40_test_ic.nc"

;-- new file names --
;file_domfile_out = '/projects/btoon/wolfet/exofiles/ocn/aquaplanet/pop_frc.gx3v7.110128.nc_0OHT.nc'
;;file_ncdata_out = '/projects/btoon/wolfet/exofiles/ocn/aquaplanet/control_L45.cam.i.0040-01-01-00000_aquaplanet.nc'
;;file_ncdata_out = '/projects/btoon/wolfet/exofiles/ocn/aquaplanet/control_L26_aquaplanet_0053-01-01-00000.nc'
;file_ncdata_out = '/projects/btoon/wolfet/exofiles/ocn/aquaplanet/aquaplanet_288K_4x5_L45.cam.i.nc'
;file_tropopause_climo_out = '/projects/btoon/wolfet/exofiles/ocn/aquaplanet/clim_p_trop_aquaplanet.nc'
;file_ozone_out = '/projects/btoon/wolfet/exofiles/ocn/aquaplanet/ozone_1.9x2.5_L26_2000clim_c091112_aquaplanet.nc'
;file_prescribaero_out = '/projects/btoon/wolfet/exofiles/ocn/aquaplanet/aero_1.9x2.5_L26_2000clim_c091112_aquaplanet.nc'
;file_aerosoldep_out = '/projects/btoon/wolfet/exofiles/ocn/aquaplanet/aerosoldep_monthly_1849-2006_1.9x2.5_c090803_aquaplanet.nc'
;file_lnd_domain_out = '/projects/btoon/wolfet/exofiles/ocn/aquaplanet/domain.lnd.fv4x5_gx3v7.091218_aquaplanet.nc'
;file_ocn_master_out = '/projects/btoon/wolfet/exofiles/ocn/aquaplanet/pop_frc.4x5d.090130_aquaplanet_Earth.nc'
;file_ocn_master_out = '/projects/btoon/wolfet/exofiles/ocn/aquaplanet/pop_frc.1x1d.090130_aquaplanet.nc'
;file_ocn_master_out = '/projects/btoon/wolfet/exofiles/ocn/aquaplanet/pop_frc.4x5d_warm_aqua_300K_1barN2.nc'

;========== BND_TOPO =====================
if (make_bndtopo eq 1) then begin

  oprstr="cp -r " + file_bndtopo + " " + file_bndtopo_out
  spawn, oprstr
  
  ; load metrics
  ncid=ncdf_open(file_bndtopo, /nowrite)
  ncdf_varget,ncid,'PHIS',PHIS
  ncdf_varget,ncid,'SGH',SGH
  ncdf_varget,ncid,'SGH30',SGH30
  ncdf_varget,ncid,'LANDFRAC',LANDFRAC
  ncdf_varget,ncid,'LANDM_COSLAT',LANDM_COSLAT
  ncdf_varget,ncid,'lon',lon
  ncdf_varget,ncid,'lat',lat
  ncdf_close,ncid

  nlon=n_elements(lon)
  nlat=n_elements(lat)

  ;set everything to zero FV
  PHIS(*,*) = 0.0
  SGH(*,*) = 0.0
  SGH30(*,*) = 0.0
  LANDFRAC(*,*) = 0.0
  LANDM_COSLAT(*,*) = 0.0

; cubed sphered
;  PHIS(*) = 0.0
;  SGH(*) = 0.0
;  SGH30(*) = 0.0
;  LANDFRAC(*) = 0.0
;  LANDM_COSLAT(*) = 0.0

  ;---- update bndtopo file ----
  ncid = ncdf_open(file_bndtopo_out, /WRITE)
  print, "updating bndtopo"
  print, file_bndtopo_out
  ncdf_varput, ncid, 'PHIS',PHIS
  ncdf_varput, ncid, 'SGH',SGH
  ncdf_varput, ncid, 'SGH30',SGH30
  ncdf_varput, ncid, 'LANDFRAC',LANDFRAC
  ncdf_varput, ncid, 'LANDM_COSLAT',LANDM_COSLAT
  ncdf_close, ncid

endif

if (make_popfrc eq 1) then begin
  oprstr="cp -r " + file_pop_frc + " " + file_pop_frc_out
  spawn, oprstr

  ncid=ncdf_open(file_pop_frc, /nowrite)
  ncdf_varget,ncid,'mask',mask
  ncdf_varget,ncid,'T',T     ; temperature
  ncdf_varget,ncid,'hblt',hblt  ; boundary later depth
  ncdf_varget,ncid,'qdp',qdp    ; ocean heat flux convergence
  ncdf_varget,ncid,'dhdx',dhdx  ; ocean surfae slope zonal
  ncdf_varget,ncid,'dhdy',dhdy  ; ocean surface slope meridional
  ncdf_varget,ncid,'U',U  ; u ocean current
  ncdf_varget,ncid,'V',V  ; v ocean current
  ncdf_varget,ncid,'S',S  ; salinity ppt
  ncdf_close,ncid

  ;set everything to zero
;  mask(*,*) = 1.0   
  T(*,*,*) = -73.16  ; degrees celcius
  hblt(*,*,*) = 50.0 ; meter mixed layer depth
  qdp(*,*,*) = 0.0
  dhdx(*,*,*) = 0.0
  dhdy(*,*,*) = 0.0
  U(*,*,*) = 0.0
  V(*,*,*) = 0.0
  S(*,*,*) = 0.0  

  ;---- update dom file ----
  ncid = ncdf_open(file_pop_frc_out, /WRITE)
  print, "updating pop_frc_out"
  print, file_pop_frc_out
;  ncdf_varput, ncid, 'mask',mask
  ncdf_varput, ncid, 'T',T
  ncdf_varput, ncid, 'hblt',hblt
  ncdf_varput, ncid, 'qdp',qdp
  ncdf_varput, ncid, 'dhdx',dhdx
  ncdf_varput, ncid, 'dhdy',dhdy
  ncdf_varput, ncid, 'U',U
  ncdf_varput, ncid, 'V',V
  ncdf_varput, ncid, 'S',S
  ncdf_close, ncid

endif

if (make_ncdata eq 1) then begin
  oprstr="cp -r " + file_ncdata + " " + file_ncdata_out
  spawn, oprstr


  ncid=ncdf_open(file_ncdata, /nowrite)
  ncdf_varget,ncid,'lon',lon
  ncdf_varget,ncid,'lat',lat
  ncdf_varget,ncid,'slat',slat
  ncdf_varget,ncid,'lev',lev
  ncdf_varget,ncid,'P0',P0
  ncdf_varget,ncid,'CLDICE',CLDICE ;[lev,lat,lon]
  ncdf_varget,ncid,'CLDLIQ',CLDLIQ ;[lev,lat,lon]
  ncdf_varget,ncid,'ICEFRAC',ICEFRAC ;[lat,lon]
  ncdf_varget,ncid,'PS',PS ;[lat,lon]
  ncdf_varget,ncid,'Q',Q ;[lev,lat,lon]
  ncdf_varget,ncid,'SICTHK',SICTHK ;[lat,lon]
  ncdf_varget,ncid,'SNOWHICE',SNOWHICE ;[lat,lon
  ncdf_varget,ncid,'T',T ;[lev,lat, lon]
;  ncdf_varget,ncid,'TS',TS ;[lat, lon]
  ncdf_varget,ncid,'TS1',TS1 ;[lat,lon]
  ncdf_varget,ncid,'TS2',TS2 ;[lat,lon]
  ncdf_varget,ncid,'TS3',TS3 ;[lat,lon]
  ncdf_varget,ncid,'TS4',TS4 ;[lat,lon]
  ncdf_varget,ncid,'TSICE',TSICE ;[lat,lon]
  ncdf_varget,ncid,'US',US ;[lev, slat,lon]
  ncdf_varget,ncid,'VS',VS ;[lev,lat,lon]
  ncdf_close,ncid

  nlat=n_elements(lat)
  nlon=n_elements(lon)
  nlev=n_elements(lev)
  nslat=n_elements(slat)

  PS_OUT = fltarr(nlon,nlat)
  TSICE_OUT = fltarr(nlon,nlat)
  SICTHK_OUT = fltarr(nlon,nlat)
  SNOWHICE_OUT = fltarr(nlon,nlat)
  ICEFRAC_OUT = fltarr(nlon,nlat)
  TS_OUT = fltarr(nlon,nlat)
  TS1_OUT = fltarr(nlon,nlat)
  TS2_OUT = fltarr(nlon,nlat)
  TS3_OUT = fltarr(nlon,nlat)
  TS4_OUT = fltarr(nlon,nlat)
  CLDICE_OUT = fltarr(nlon,nlat,nlev)
  CLDLIQ_OUT = fltarr(nlon,nlat,nlev)
  Q_OUT = fltarr(nlon,nlat,nlev)
  T_OUT = fltarr(nlon,nlat,nlev)
  VS_OUT = fltarr(nlon,nlat,nlev)
  US_OUT = fltarr(nlon,nslat,nlev)

  ;set variables
  P0_OUT=1.e5
  PS_OUT(*,*) =1.e5
  TSICE_OUT(*,*)=0.0   ;273.15
  ICEFRAC_OUT(*,*)=0.0
  SICTHK_OUT(*,*)=0.0
  SNOWHICE_OUT(*,*)=0.0
  T_OUT(*,*,*) = 200.0
  TS1_OUT(*,*) = 200.0
  TS2_OUT(*,*) = 200.0
  TS3_OUT(*,*) = 200.0
  TS4_OUT(*,*) = 200.0
  CLDICE_OUT(*,*,*) = 0.0
  CLDLIQ_OUT(*,*,*) = 0.0
  VS_OUT(*,*,*) = 0.0
  US_OUT(*,*,*) = 0.0


  iv=nlat
  for i=0,nlat-1 do begin
    iv=iv-1
    ;TS_OUT(*,i) = mean((TS(*,i)+TS(*,iv))/2.)
    TS1_OUT(*,i) = mean((TS1(*,i)+TS1(*,iv))/2.)
    TS2_OUT(*,i) = mean((TS2(*,i)+TS2(*,iv))/2.)
    TS3_OUT(*,i) = mean((TS3(*,i)+TS3(*,iv))/2.)
    TS4_OUT(*,i) = mean((TS4(*,i)+TS4(*,iv))/2.)
    for k=0,nlev-1 do begin
    ;  CLDICE_OUT(*,i,k) = mean((CLDICE(*,i,k)+CLDICE(*,iv,*))/2.)
    ;  CLDLIQ_OUT(*,i,k) = mean((CLDLIQ(*,i,k)+CLDLIQ(*,iv,k))/2.)
      Q_OUT(*,i,k) = mean((Q(*,i,k)+Q(*,iv,k))/2.)
      T_OUT(*,i,k) = mean((T(*,i,k)+T(*,iv,k))/2.)
;      VS_OUT(*,i,k) = mean((VS(*,i,k)+VS(*,iv,k))/2.)
    endfor
   endfor

;  iv=nslat
;  for i=0,nslat-1 do begin
;    iv=iv-1
;    for k=0,nlev-1 do begin
;  ;    US_OUT(*,i,k) = mean((US(*,i,k)+US(*,iv,k))/2.)
;    endfor
;  endfor



  ;---- update ncdata file ----
  ncid = ncdf_open(file_ncdata_out, /WRITE)
  print, "updatinf ncdata file"
  print, file_ncdata_out
  ncdf_varput, ncid, 'P0',P0_OUT
  ncdf_varput, ncid, 'PS',PS_OUT
 ; ncdf_varput, ncid, 'TSICE',TSICE_OUT
 ; ncdf_varput, ncid, 'SICTHK',SICTHK_OUT
 ; ncdf_varput, ncid, 'TS',TS_OUT
  ncdf_varput, ncid, 'TS1',TS1_OUT
  ncdf_varput, ncid, 'TS2',TS2_OUT
  ncdf_varput, ncid, 'TS3',TS3_OUT
  ncdf_varput, ncid, 'TS4',TS4_OUT
  ncdf_varput, ncid, 'CLDICE',CLDICE_OUT
  ncdf_varput, ncid, 'CLDLIQ',CLDLIQ_OUT
  ncdf_varput, ncid, 'Q',Q_OUT
  ncdf_varput, ncid, 'T',T_OUT
  ncdf_varput, ncid, 'VS',VS_OUT
  ncdf_varput, ncid, 'US',US_OUT
  ncdf_close, ncid




endif

if (make_lnd_domain eq 1) then begin

  oprstr="cp -r " + file_lnd_domain + " " + file_lnd_domain_out
  spawn, oprstr
  
  ; load metrics
  ncid=ncdf_open(file_lnd_domain, /nowrite)
  ncdf_varget,ncid,'mask',mask
  ncdf_varget,ncid,'frac',frac
  ncdf_close,ncid

  ;set everything to zero
  mask(*,*) = 0.0
  frac(*,*) = 0.0

  ;---- update lnd domain file ----
  ncid = ncdf_open(file_lnd_domain_out, /WRITE)
  print, "updating lnd_domain"
  print, file_lnd_domain_out
  ncdf_varput, ncid, 'mask',mask
  ncdf_varput, ncid, 'frac',frac
  ncdf_close, ncid
endif


if (make_ocn_domain eq 1) then begin

  oprstr="cp -r " + file_ocn_domain + " " + file_ocn_domain_out
  spawn, oprstr
  
  ; load metrics
  ncid=ncdf_open(file_ocn_domain, /nowrite)
  ncdf_varget,ncid,'mask',mask
  ncdf_varget,ncid,'frac',frac
  ncdf_close,ncid


  mask(*,*) = 1
  frac(*,*) = 1.0

  ;---- update ocn domain file ----
  ncid = ncdf_open(file_ocn_domain_out, /WRITE)
  print, "updating ocn_domain"
  print, file_ocn_domain_out
  ncdf_varput, ncid, 'mask',mask
  ncdf_varput, ncid, 'frac',frac
  ncdf_close, ncid
endif


; This option is different since the ocean 1x1d file must be coarsed
; down to 4x5d

if (make_ocn_master eq 1) then begin

  ;optiona read in .cam.i. file for exact temperature match
  copyT_file = "/projects/btoon/wolfet/exofiles/atm/warm_aqua_300K_1barN2.cam.i.nc"
  ;copyT_file = "/gpfs/summit/datasets/CESM/inputdata/atm/cam/inic/fv/cami_0000-09-01_0.47x0.63_L26_c061106.nc"
  ncid=ncdf_open(copyT_file, /nowrite)
  ncdf_varget,ncid,'TS',T_match_in
  ncdf_close,ncid
help, T_match_in

  ; get some data for domain
  ncid=ncdf_open(file_ocn_master, /nowrite)
  ncdf_varget,ncid,'xc',xc_in
  ncdf_varget,ncid,'yc',yc_in
  ncdf_varget,ncid,'area',area_in
  ncdf_close,ncid

  ;dimensions
  ;kludge
  ;nlon = 576
  ;nlat = 384
  ntime = 12

  ;data
  area_out = fltarr(nlon, nlat)
  mask_out = intarr(nlon, nlat)
  yc_out = fltarr(nlat)
  xc_out = fltarr(nlon)
  time_out = fltarr(ntime)
  S_out = fltarr(nlon, nlat,ntime)
  T_out = fltarr(nlon, nlat,ntime)
  U_out = fltarr(nlon, nlat,ntime)
  V_out = fltarr(nlon, nlat,ntime)
  dhdx_out = fltarr(nlon, nlat,ntime)
  dhdy_out = fltarr(nlon, nlat,ntime)
  hblt_out = fltarr(nlon, nlat,ntime)
  qdp_out = fltarr(nlon, nlat,ntime)


  time_out = [14, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]
  area_out = area_in
  yc_out(*) = yc_in(0,*)
  xc_out(*)= xc_in(*,0)

; uncomment for T_match
  for x=0,nlon-1 do begin
    for y=0,nlat-1 do begin 
     ; T_out(x,y,*) = T(x,y,nlev-1)-273.15
     T_out(x,y,*) = T_match_in(x,y)-273.15
    endfor
 endfor

  mask_out(*,*) = 1
  ;T_out(*,*,*) = 15.0  ; degrees celcius
  hblt_out(*,*,*) = 50.0 ; meter mixed layer depth
  qdp_out(*,*,*) = 0.0
  dhdx_out(*,*,*) = 0.0
  dhdy_out(*,*,*) = 0.0
  U_out(*,*,*) = 0.0
  V_out(*,*,*) = 0.0
  S_out(*,*,*) = 0.0  

  outname = file_ocn_master_out

   print, "update ", outname
  id = NCDF_CREATE(outname, /CLOBBER)
  dim1 = NCDF_DIMDEF(id,'time',ntime)
  dim2 = NCDF_DIMDEF(id,'lat',nlat)
  dim3 = NCDF_DIMDEF(id,'lon',nlon)

  varid1 = NCDF_VARDEF(id, 'area',[dim3,dim2], /float)
  varid2 = NCDF_VARDEF(id, 'mask', [dim3,dim2])
  varid3 = NCDF_VARDEF(id, 'xc', [dim3], /float)
  varid4 = NCDF_VARDEF(id, 'yc', [dim2], /float)
  varid5 = NCDF_VARDEF(id, 'time', [dim1], /float)
  varid6 = NCDF_VARDEF(id, 'S', [dim3,dim2,dim1], /float)
  varid7 = NCDF_VARDEF(id, 'T', [dim3,dim2,dim1], /float)
  varid8 = NCDF_VARDEF(id, 'U', [dim3,dim2,dim1], /float)
  varid9 = NCDF_VARDEF(id, 'V', [dim3,dim2,dim1], /float)
  varid10 = NCDF_VARDEF(id, 'dhdx', [dim3,dim2,dim1], /float)
  varid11 = NCDF_VARDEF(id, 'dhdy', [dim3,dim2,dim1], /float)
  varid12 = NCDF_VARDEF(id, 'hblt', [dim3,dim2,dim1], /float)
  varid13 = NCDF_VARDEF(id, 'qdp', [dim3,dim2,dim1], /float)


  NCDF_ATTPUT, id, varid1, "long_name", "area of grid cell"
  NCDF_ATTPUT, id, varid1, "units", "area"
  NCDF_ATTPUT, id, varid2, "long_name", "domain maskr"
  NCDF_ATTPUT, id, varid2, "units", "unitless"
;  NCDF_ATTPUT, id, varid2, "missing_value",-2147483647, /long
;  NCDF_ATTPUT, id, varid2, "_FillValue",-2147483647, /long
  NCDF_ATTPUT, id, varid3, "long_name", "degrees east"
  NCDF_ATTPUT, id, varid3, "units", "longitude of grid cell center"
  NCDF_ATTPUT, id, varid4, "long_name", "degrees north"
  NCDF_ATTPUT, id, varid4, "units", "latitude of grid cell center"
  NCDF_ATTPUT, id, varid5, "long_name", "observationb time"
  NCDF_ATTPUT, id, varid5, "calendar", "noleap"
  NCDF_ATTPUT, id, varid5, "units", "days since 0001-01-01 00:00:00"
  NCDF_ATTPUT, id, varid6, "long_name", "salinity"
  NCDF_ATTPUT, id, varid6, "units", "ppt"
;  NCDF_ATTPUT, id, varid6, "_FillValue",-999.0, /float
  NCDF_ATTPUT, id, varid7, "long_name", "temperature"
  NCDF_ATTPUT, id, varid7, "units", "degC"
;  NCDF_ATTPUT, id, varid7, "_FillValue",-999.0, /float
  NCDF_ATTPUT, id, varid8, "long_name", "u ocean current"
  NCDF_ATTPUT, id, varid8, "units", "m/s"
;  NCDF_ATTPUT, id, varid8, "_FillValue",-999.0, /float
  NCDF_ATTPUT, id, varid9, "long_name", "v ocean current"
  NCDF_ATTPUT, id, varid9, "units", "m/s"
;  NCDF_ATTPUT, id, varid9, "_FillValue",-999.0, /float
  NCDF_ATTPUT, id, varid10, "long_name", "ocean surface slope: zonal"
  NCDF_ATTPUT, id, varid10, "units", "m/m"
;  NCDF_ATTPUT, id, varid10, "_FillValue",-999.0, /float
  NCDF_ATTPUT, id, varid11, "long_name", "ocean surface slope: meridional"
  NCDF_ATTPUT, id, varid11, "units", "m/m"
;  NCDF_ATTPUT, id, varid11, "_FillValue",-999.0, /float
  NCDF_ATTPUT, id, varid12, "long_name", "boundary layer depth"
  NCDF_ATTPUT, id, varid12, "units", "m"
;  NCDF_ATTPUT, id, varid12, "_FillValue",-999.0, /float
  NCDF_ATTPUT, id, varid13, "long_name", "ocean heat flux convergence"
  NCDF_ATTPUT, id, varid13, "units", "W/m^2"
;  NCDF_ATTPUT, id, varid13, "_FillValue",-999.0, /float



  NCDF_CONTROL, id, /ENDEF

  NCDF_VARPUT, id, varid1, area_out
  NCDF_VARPUT, id, varid2, mask_out
  NCDF_VARPUT, id, varid3, xc_out
  NCDF_VARPUT, id, varid4, yc_out
  NCDF_VARPUT, id, varid5, time_out
  NCDF_VARPUT, id, varid6, S_out
  NCDF_VARPUT, id, varid7, T_out
  NCDF_VARPUT, id, varid8, U_out
  NCDF_VARPUT, id, varid9, V_out
  NCDF_VARPUT, id, varid10, dhdx_out
  NCDF_VARPUT, id, varid11, dhdy_out
  NCDF_VARPUT, id, varid12, hblt_out
  NCDF_VARPUT, id, varid13, qdp_out

  NCDF_CLOSE, id



endif




if (make_tropopause_climo_file eq 1) then begin
  oprstr="cp -r " + file_tropopause_climo + " " + file_tropopause_climo_out
  spawn, oprstr

  ; updating tropopause climo
  ncid = ncdf_open(file_tropopause_climo_out, /WRITE)
  ncdf_varget, ncid, 'month',month
  ncdf_varget, ncid, 'lat',lat
  ncdf_varget, ncid, 'lon',lon
  ncdf_varget, ncid, 'trop_p',trop_p
  ncdf_close, ncid

  trop_p_out=fltarr(n_elements(lon), n_elements(lat), n_elements(month))


  for d=0,n_elements(month)-1 do begin
  iv=n_elements(lat)
  for i=0,n_elements(lat)-1 do begin
     iv=iv-1
     trop_p_out(*,i,d) =  mean((trop_p(*,i,d)+trop_p(*,iv,d))/2.0)
  endfor
  endfor

  ;---- update ozone domain file ----
  ncid = ncdf_open(file_tropopause_climo_out, /WRITE)
  print, "updating tropoapuse_climo"
  print, file_tropopause_climo_out
  ncdf_varput, ncid, 'trop_p',trop_p_out
  ncdf_close, ncid


endif

if (make_ozone eq 1) then begin
  oprstr="cp -r " + file_ozone + " " + file_ozone_out
  spawn, oprstr

  ncid=ncdf_open(file_ozone, /nowrite)
  ncdf_varget,ncid,'O3',O3
  ncdf_close,ncid

  O3(*,*,*,*) = 1.0e-20
  ;---- update ozone domain file ----
  ncid = ncdf_open(file_ozone_out, /WRITE)
  print, "updating ozone"
  print, file_ozone_out
  ncdf_varput, ncid, 'O3',O3
  ncdf_close, ncid

endif

if (make_aerosoldep eq 1) then begin
  oprstr="cp -r " + file_aerosoldep + " " + file_aerosoldep_out
  spawn, oprstr


  ncid=ncdf_open(file_aerosoldep, /nowrite)
  ncdf_varget,ncid,'BCDEPWET',BCDEPWET & BCDEPWET(*,*,*) = 0.0
  ncdf_varget,ncid,'BCPHODRY',BCPHODRY & BCPHODRY(*,*,*) = 0.0
  ncdf_varget,ncid,'BCPHIDRY',BCPHIDRY & BCPHIDRY(*,*,*) = 0.0
  ncdf_varget,ncid,'OCDEPWET',OCDEPWET & OCDEPWET(*,*,*) = 0.0
  ncdf_varget,ncid,'OCPHODRY',OCPHODRY & OCPHODRY(*,*,*) = 0.0
  ncdf_varget,ncid,'OCPHIDRY',OCPHIDRY & OCPHIDRY(*,*,*) = 0.0
  ncdf_varget,ncid,'DSTX01DD',DSTX01DD & DSTX01DD(*,*,*) = 0.0
  ncdf_varget,ncid,'DSTX02DD',DSTX02DD & DSTX02DD(*,*,*) = 0.0
  ncdf_varget,ncid,'DSTX03DD',DSTX03DD & DSTX03DD(*,*,*) = 0.0
  ncdf_varget,ncid,'DSTX04DD',DSTX04DD & DSTX04DD(*,*,*) = 0.0
  ncdf_varget,ncid,'DSTX01WD',DSTX01WD & DSTX01WD(*,*,*) = 0.0
  ncdf_varget,ncid,'DSTX02WD',DSTX02WD & DSTX02WD(*,*,*) = 0.0
  ncdf_varget,ncid,'DSTX03WD',DSTX03WD & DSTX03WD(*,*,*) = 0.0
  ncdf_varget,ncid,'DSTX04WD',DSTX04WD & DSTX04WD(*,*,*) = 0.0
  ncdf_close,ncid



  ;---- update aerosoldep file ----
  ncid = ncdf_open(file_aerosoldep_out, /WRITE)
  print, "updating aerosoldep"
  print, file_aerosoldep_out
  ncdf_varput,ncid,'BCDEPWET',BCDEPWET 
  ncdf_varput,ncid,'BCPHODRY',BCPHODRY 
  ncdf_varput,ncid,'BCPHIDRY',BCPHIDRY 
  ncdf_varput,ncid,'OCDEPWET',OCDEPWET 
  ncdf_varput,ncid,'OCPHODRY',OCPHODRY 
  ncdf_varput,ncid,'OCPHIDRY',OCPHIDRY 
  ncdf_varput,ncid,'DSTX01DD',DSTX01DD 
  ncdf_varput,ncid,'DSTX02DD',DSTX02DD 
  ncdf_varput,ncid,'DSTX03DD',DSTX03DD 
  ncdf_varput,ncid,'DSTX04DD',DSTX04DD 
  ncdf_varput,ncid,'DSTX01WD',DSTX01WD 
  ncdf_varput,ncid,'DSTX02WD',DSTX02WD 
  ncdf_varput,ncid,'DSTX03WD',DSTX03WD 
  ncdf_varput,ncid,'DSTX04WD',DSTX04WD 
  ncdf_close, ncid
endif

if (make_prescribed_aerofil eq 1) then begin
  oprstr="cp -r " + file_prescribed_aero + " " + file_prescribed_aero_out
  spawn, oprstr

  ncid=ncdf_open(file_prescribed_aero, /nowrite)
  ncdf_varget,ncid,'CB1',CB1 & CB1(*,*,*,*) = 0.0
  ncdf_varget,ncid,'CB2',CB2 & CB2(*,*,*,*) = 0.0
  ncdf_varget,ncid,'DST01',DST01 & DST01(*,*,*,*) = 0.0
  ncdf_varget,ncid,'DST02',DST02 & DST02(*,*,*,*) = 0.0
  ncdf_varget,ncid,'DST03',DST03 & DST03(*,*,*,*) = 0.0
  ncdf_varget,ncid,'DST04',DST04 & DST04(*,*,*,*) = 0.0
  ncdf_varget,ncid,'OC1',OC1 & OC1(*,*,*,*) = 0.0
  ncdf_varget,ncid,'OC2',OC2 & OC2(*,*,*,*) = 0.0
  ncdf_varget,ncid,'SO4',SO4 & SO4(*,*,*,*) = 0.0
  ncdf_varget,ncid,'SSLT01',SSLT01 & SSLT01(*,*,*,*) = 0.0
  ncdf_varget,ncid,'SSLT02',SSLT02 & SSLT02(*,*,*,*) = 0.0
  ncdf_varget,ncid,'SSLT03',SSLT03 & SSLT03(*,*,*,*) = 0.0
  ncdf_varget,ncid,'SSLT04',SSLT04 & SSLT04(*,*,*,*) = 0.0
  ncdf_close,ncid

  ;---- update prescribed file ----
  ncid = ncdf_open(file_prescribed_aero_out, /WRITE)
  print, "updating prescribed_aero"
  print, file_prescribed_aero_out
  ncdf_varput,ncid,'CB1',CB1 
  ncdf_varput,ncid,'CB2',CB2
  ncdf_varput,ncid,'DST01',DST01
  ncdf_varput,ncid,'DST02',DST02
  ncdf_varput,ncid,'DST03',DST03
  ncdf_varput,ncid,'DST04',DST04
  ncdf_varput,ncid,'OC1',OC1
  ncdf_varput,ncid,'OC2',OC2
  ncdf_varput,ncid,'SO4',SO4
  ncdf_varput,ncid,'SSLT01',SSLT01
  ncdf_varput,ncid,'SSLT02',SSLT02
  ncdf_varput,ncid,'SSLT03',SSLT03
  ncdf_varput,ncid,'SSLT04',SSLT04
  ncdf_close, ncid
endif


if (make_ne16np4_popfile eq 1) then begin
;  oprstr="cp -r " + file_domfile + " " + file_domfile_out
;  spawn, oprstr
  ; get some data from domain file
  file = '/lustre/janus_scratch/cesm/inputdata/share/domains/domain.ocn.ne16np4_gx3v7.120406.nc'
  ncid=ncdf_open(file, /nowrite)
  ncdf_varget,ncid,'xc',xc_in
  ncdf_varget,ncid,'yc',yc_in
  ncdf_varget,ncid,'xv',xv_in
  ncdf_varget,ncid,'yv',yv_in
  ncdf_varget,ncid,'area',area_in
  ncdf_varget,ncid,'mask',mask_in
  ncdf_close,ncid

  ;dimensions
  n = 13826
  ni = 13826
  nj = 1
  nv = 4
  ntime=1

  ;data
;  area_out = fltarr(nlon, nlat)
;  mask_out = intarr(nlon, nlat)
;  yc_out = fltarr(nlat)
;  xc_out = fltarr(nlon)
;  time_out = fltarr(ntime)
  S_out = fltarr(ni, ntime)  & S_out(*,*,*) = 0.0
  T_out = fltarr(ni, ntime)  & T_out(*,*,*) = 10.0
  U_out = fltarr(ni, ntime)  & U_out(*,*,*) = 0.0
  V_out = fltarr(ni, ntime)  & V_out(*,*,*) = 0.0
  dhdx_out = fltarr(ni, ntime) & dhdx_out(*,*,*) = 0.0
  dhdy_out = fltarr(ni, ntime) & dhdy_out(*,*,*) = 0.0
  hblt_out = fltarr(ni, ntime) & hblt_out(*,*,*) = 50.0
  qdp_out = fltarr(ni, ntime)  & qdp_out(*,*,*) = 0.0


;  time_out = [14, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]
  time_out = 0.
  xc_out= xc_in
  yc_out = yc_in
  xv_out= xv_in
  yv_out = yv_in
  mask_out = mask_in
  mask_out(*) = 1.0
  area_out = area_in




  outname = "ocean_ne16np4.nc"

   print, "update ", outname
  id = NCDF_CREATE(outname, /CLOBBER)
  dim1 = NCDF_DIMDEF(id,'n',n)
  dim2 = NCDF_DIMDEF(id,'ni',ni)
;  dim3 = NCDF_DIMDEF(id,'nj',nj)
  dim4 = NCDF_DIMDEF(id,'nv',nv)
  dim5 = NCDF_DIMDEF(id,'time',ntime)

  varid1 = NCDF_VARDEF(id, 'area',[dim2], /float)
  varid2 = NCDF_VARDEF(id, 'mask', [dim2])
  varid3 = NCDF_VARDEF(id, 'xc', [dim2], /float)
  varid4 = NCDF_VARDEF(id, 'yc', [dim2], /float)
  varid5 = NCDF_VARDEF(id, 'xv', [dim4,dim2], /float)
  varid6 = NCDF_VARDEF(id, 'yv', [dim4,dim2], /float)
  varid7 = NCDF_VARDEF(id, 'time', [dim5], /float)
  varid8 = NCDF_VARDEF(id, 'S', [dim2,dim5], /float)
  varid9 = NCDF_VARDEF(id, 'T', [dim2,dim5], /float)
  varid10 = NCDF_VARDEF(id, 'U', [dim2,dim5], /float)
  varid11 = NCDF_VARDEF(id, 'V', [dim2,dim5], /float)
  varid12 = NCDF_VARDEF(id, 'dhdx', [dim2,dim5], /float)
  varid13 = NCDF_VARDEF(id, 'dhdy', [dim2,dim5], /float)
  varid14 = NCDF_VARDEF(id, 'hblt', [dim2,dim5], /float)
  varid15 = NCDF_VARDEF(id, 'qdp', [dim2,dim5], /float)


  NCDF_ATTPUT, id, varid1, "long_name", "area of grid cell"
  NCDF_ATTPUT, id, varid1, "units", "area"
  NCDF_ATTPUT, id, varid2, "long_name", "domain maskr"
  NCDF_ATTPUT, id, varid2, "units", "unitless"
;  NCDF_ATTPUT, id, varid2, "missing_value",-2147483647, /long
;  NCDF_ATTPUT, id, varid2, "_FillValue",-2147483647, /long
  NCDF_ATTPUT, id, varid3, "long_name", "degrees east"
  NCDF_ATTPUT, id, varid3, "units", "longitude of grid cell center"
  NCDF_ATTPUT, id, varid4, "long_name", "degrees north"
  NCDF_ATTPUT, id, varid4, "units", "latitude of grid cell center"

  NCDF_ATTPUT, id, varid5, "long_name", "degrees east"
  NCDF_ATTPUT, id, varid5, "units", "longitude of grid cell vertices"
  NCDF_ATTPUT, id, varid6, "long_name", "degrees north"
  NCDF_ATTPUT, id, varid6, "units", "latitude of grid cell vertices"



  NCDF_ATTPUT, id, varid7, "long_name", "observationb time"
  NCDF_ATTPUT, id, varid7, "calendar", "noleap"
  NCDF_ATTPUT, id, varid7, "units", "days since 0001-01-01 00:00:00"
  NCDF_ATTPUT, id, varid8, "long_name", "salinity"
  NCDF_ATTPUT, id, varid8, "units", "ppt"
;  NCDF_ATTPUT, id, varid8, "_FillValue",-999.0, /float
  NCDF_ATTPUT, id, varid9, "long_name", "temperature"
  NCDF_ATTPUT, id, varid9, "units", "degC"
;  NCDF_ATTPUT, id, varid9, "_FillValue",-999.0, /float
  NCDF_ATTPUT, id, varid10, "long_name", "u ocean current"
  NCDF_ATTPUT, id, varid10, "units", "m/s"
;  NCDF_ATTPUT, id, varid10, "_FillValue",-999.0, /float
  NCDF_ATTPUT, id, varid11, "long_name", "v ocean current"
  NCDF_ATTPUT, id, varid11, "units", "m/s"
;  NCDF_ATTPUT, id, varid11, "_FillValue",-999.0, /float
  NCDF_ATTPUT, id, varid12, "long_name", "ocean surface slope: zonal"
  NCDF_ATTPUT, id, varid12, "units", "m/m"
;  NCDF_ATTPUT, id, varid12, "_FillValue",-999.0, /float
  NCDF_ATTPUT, id, varid13, "long_name", "ocean surface slope: meridional"
  NCDF_ATTPUT, id, varid13, "units", "m/m"
;  NCDF_ATTPUT, id, varid13, "_FillValue",-999.0, /float
  NCDF_ATTPUT, id, varid14, "long_name", "boundary layer depth"
  NCDF_ATTPUT, id, varid14, "units", "m"
;  NCDF_ATTPUT, id, varid14, "_FillValue",-999.0, /float
  NCDF_ATTPUT, id, varid15, "long_name", "ocean heat flux convergence"
  NCDF_ATTPUT, id, varid15, "units", "W/m^2"
;  NCDF_ATTPUT, id, varid15, "_FillValue",-999.0, /float



  NCDF_CONTROL, id, /ENDEF

  NCDF_VARPUT, id, varid1, area_out
  NCDF_VARPUT, id, varid2, mask_out
  NCDF_VARPUT, id, varid3, xc_out
  NCDF_VARPUT, id, varid4, yc_out
  NCDF_VARPUT, id, varid5, xv_out
  NCDF_VARPUT, id, varid6, yv_out
  NCDF_VARPUT, id, varid7, time_out
  NCDF_VARPUT, id, varid8, S_out
  NCDF_VARPUT, id, varid9, T_out
  NCDF_VARPUT, id, varid10, U_out
  NCDF_VARPUT, id, varid11, V_out
  NCDF_VARPUT, id, varid12, dhdx_out
  NCDF_VARPUT, id, varid13, dhdy_out
  NCDF_VARPUT, id, varid14, hblt_out
  NCDF_VARPUT, id, varid15, qdp_out

  NCDF_CLOSE, id


endif



if (make_ne5np4_popfile eq 1) then begin
;  oprstr="cp -r " + file_domfile + " " + file_domfile_out
;  spawn, oprstr
  ; get some data from domain file
  file = "/projects/btoon/wolfet/exofiles/ocn/ne5np4/domain.ocn.ne5np4_gx3v7.140810.nc"
  ncid=ncdf_open(file, /nowrite)
  ncdf_varget,ncid,'xc',xc_in
  ncdf_varget,ncid,'yc',yc_in
  ncdf_varget,ncid,'xv',xv_in
  ncdf_varget,ncid,'yv',yv_in
  ncdf_varget,ncid,'area',area_in
  ncdf_varget,ncid,'mask',mask_in
  ncdf_close,ncid

  ;dimensions
  n = 1352
  ni = 1352
  nj = 1
  nv = 4
  ntime=1

  ;data
;  area_out = fltarr(nlon, nlat)
;  mask_out = intarr(nlon, nlat)
;  yc_out = fltarr(nlat)
;  xc_out = fltarr(nlon)
;  time_out = fltarr(ntime)
  S_out = fltarr(ni, ntime)  & S_out(*,*,*) = 0.0
  T_out = fltarr(ni, ntime)  & T_out(*,*,*) = 10.0
  U_out = fltarr(ni, ntime)  & U_out(*,*,*) = 0.0
  V_out = fltarr(ni, ntime)  & V_out(*,*,*) = 0.0
  dhdx_out = fltarr(ni, ntime) & dhdx_out(*,*,*) = 0.0
  dhdy_out = fltarr(ni, ntime) & dhdy_out(*,*,*) = 0.0
  hblt_out = fltarr(ni, ntime) & hblt_out(*,*,*) = 50.0
  qdp_out = fltarr(ni, ntime)  & qdp_out(*,*,*) = 0.0


;  time_out = [14, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]
  time_out = 0.
  xc_out= xc_in
  yc_out = yc_in
  xv_out= xv_in
  yv_out = yv_in
  mask_out = mask_in
  mask_out(*) = 1.0
  area_out = area_in

  outname = "ocean_ne5np4.nc"

   print, "update ", outname
  id = NCDF_CREATE(outname, /CLOBBER)
  dim1 = NCDF_DIMDEF(id,'n',n)
  dim2 = NCDF_DIMDEF(id,'ni',ni)
;  dim3 = NCDF_DIMDEF(id,'nj',nj)
  dim4 = NCDF_DIMDEF(id,'nv',nv)
  dim5 = NCDF_DIMDEF(id,'time',ntime)

  varid1 = NCDF_VARDEF(id, 'area',[dim2], /float)
  varid2 = NCDF_VARDEF(id, 'mask', [dim2])
  varid3 = NCDF_VARDEF(id, 'xc', [dim2], /float)
  varid4 = NCDF_VARDEF(id, 'yc', [dim2], /float)
  varid5 = NCDF_VARDEF(id, 'xv', [dim4,dim2], /float)
  varid6 = NCDF_VARDEF(id, 'yv', [dim4,dim2], /float)
  varid7 = NCDF_VARDEF(id, 'time', [dim5], /float)
  varid8 = NCDF_VARDEF(id, 'S', [dim2,dim5], /float)
  varid9 = NCDF_VARDEF(id, 'T', [dim2,dim5], /float)
  varid10 = NCDF_VARDEF(id, 'U', [dim2,dim5], /float)
  varid11 = NCDF_VARDEF(id, 'V', [dim2,dim5], /float)
  varid12 = NCDF_VARDEF(id, 'dhdx', [dim2,dim5], /float)
  varid13 = NCDF_VARDEF(id, 'dhdy', [dim2,dim5], /float)
  varid14 = NCDF_VARDEF(id, 'hblt', [dim2,dim5], /float)
  varid15 = NCDF_VARDEF(id, 'qdp', [dim2,dim5], /float)


  NCDF_ATTPUT, id, varid1, "long_name", "area of grid cell"
  NCDF_ATTPUT, id, varid1, "units", "area"
  NCDF_ATTPUT, id, varid2, "long_name", "domain maskr"
  NCDF_ATTPUT, id, varid2, "units", "unitless"
;  NCDF_ATTPUT, id, varid2, "missing_value",-2147483647, /long
;  NCDF_ATTPUT, id, varid2, "_FillValue",-2147483647, /long
  NCDF_ATTPUT, id, varid3, "long_name", "degrees east"
  NCDF_ATTPUT, id, varid3, "units", "longitude of grid cell center"
  NCDF_ATTPUT, id, varid4, "long_name", "degrees north"
  NCDF_ATTPUT, id, varid4, "units", "latitude of grid cell center"

  NCDF_ATTPUT, id, varid5, "long_name", "degrees east"
  NCDF_ATTPUT, id, varid5, "units", "longitude of grid cell vertices"
  NCDF_ATTPUT, id, varid6, "long_name", "degrees north"
  NCDF_ATTPUT, id, varid6, "units", "latitude of grid cell vertices"



  NCDF_ATTPUT, id, varid7, "long_name", "observation time"
  NCDF_ATTPUT, id, varid7, "calendar", "noleap"
  NCDF_ATTPUT, id, varid7, "units", "days since 0001-01-01 00:00:00"
  NCDF_ATTPUT, id, varid8, "long_name", "salinity"
  NCDF_ATTPUT, id, varid8, "units", "ppt"
;  NCDF_ATTPUT, id, varid8, "_FillValue",-999.0, /float
  NCDF_ATTPUT, id, varid9, "long_name", "temperature"
  NCDF_ATTPUT, id, varid9, "units", "degC"
;  NCDF_ATTPUT, id, varid9, "_FillValue",-999.0, /float
  NCDF_ATTPUT, id, varid10, "long_name", "u ocean current"
  NCDF_ATTPUT, id, varid10, "units", "m/s"
;  NCDF_ATTPUT, id, varid10, "_FillValue",-999.0, /float
  NCDF_ATTPUT, id, varid11, "long_name", "v ocean current"
  NCDF_ATTPUT, id, varid11, "units", "m/s"
;  NCDF_ATTPUT, id, varid11, "_FillValue",-999.0, /float
  NCDF_ATTPUT, id, varid12, "long_name", "ocean surface slope: zonal"
  NCDF_ATTPUT, id, varid12, "units", "m/m"
;  NCDF_ATTPUT, id, varid12, "_FillValue",-999.0, /float
  NCDF_ATTPUT, id, varid13, "long_name", "ocean surface slope: meridional"
  NCDF_ATTPUT, id, varid13, "units", "m/m"
;  NCDF_ATTPUT, id, varid13, "_FillValue",-999.0, /float
  NCDF_ATTPUT, id, varid14, "long_name", "boundary layer depth"
  NCDF_ATTPUT, id, varid14, "units", "m"
;  NCDF_ATTPUT, id, varid14, "_FillValue",-999.0, /float
  NCDF_ATTPUT, id, varid15, "long_name", "ocean heat flux convergence"
  NCDF_ATTPUT, id, varid15, "units", "W/m^2"
;  NCDF_ATTPUT, id, varid15, "_FillValue",-999.0, /float



  NCDF_CONTROL, id, /ENDEF

  NCDF_VARPUT, id, varid1, area_out
  NCDF_VARPUT, id, varid2, mask_out
  NCDF_VARPUT, id, varid3, xc_out
  NCDF_VARPUT, id, varid4, yc_out
  NCDF_VARPUT, id, varid5, xv_out
  NCDF_VARPUT, id, varid6, yv_out
  NCDF_VARPUT, id, varid7, time_out
  NCDF_VARPUT, id, varid8, S_out
  NCDF_VARPUT, id, varid9, T_out
  NCDF_VARPUT, id, varid10, U_out
  NCDF_VARPUT, id, varid11, V_out
  NCDF_VARPUT, id, varid12, dhdx_out
  NCDF_VARPUT, id, varid13, dhdy_out
  NCDF_VARPUT, id, varid14, hblt_out
  NCDF_VARPUT, id, varid15, qdp_out

  NCDF_CLOSE, id


endif


if (make_ne16np4_icfile eq 1) then begin
 ; oprstr="cp -r " + file_ne16np4_icfile + " " + file_ne16np4_icfile_out
 ; spawn, oprstr


  ncid=ncdf_open(file_ne16np4_icfile, /nowrite)
  ncdf_varget,ncid,'lon',lon
  ncdf_varget,ncid,'lat',lat
  ncdf_varget,ncid,'lev',lev
  ncdf_varget,ncid,'P0',P0
  ncdf_varget,ncid,'CLDICE',CLDICE ;[lev,ncol]
  ncdf_varget,ncid,'CLDLIQ',CLDLIQ 
  ncdf_varget,ncid,'ICEFRAC',ICEFRAC ;[ncol]
  ncdf_varget,ncid,'PS',PS ;[ncol]
  ncdf_varget,ncid,'Q',Q 
  ncdf_varget,ncid,'SICTHK',SICTHK 
  ncdf_varget,ncid,'SNOWHICE',SNOWHICE 
  ncdf_varget,ncid,'T',T 
  ncdf_varget,ncid,'TS1',TS1 
  ncdf_varget,ncid,'TS2',TS2 
  ncdf_varget,ncid,'TS3',TS3 
  ncdf_varget,ncid,'TS4',TS4 
  ncdf_varget,ncid,'TSICE',TSICE 
  ncdf_varget,ncid,'U',U 
  ncdf_varget,ncid,'V',V
  ncdf_close,ncid

  nlat=n_elements(lat)
  nlon=n_elements(lon)
  nlev=n_elements(lev)

;fix this
;  PS_OUT = fltarr(nlon,nlat)
;  TSICE_OUT = fltarr(nlon,nlat)
;  SICTHK_OUT = fltarr(nlon,nlat)
;  SNOWHICE_OUT = fltarr(nlon,nlat)
;  ICEFRAC_OUT = fltarr(nlon,nlat)
;  TS_OUT = fltarr(nlon,nlat)
;  TS1_OUT = fltarr(nlon,nlat)
;  TS2_OUT = fltarr(nlon,nlat)
;  TS3_OUT = fltarr(nlon,nlat)
;  TS4_OUT = fltarr(nlon,nlat)
;  CLDICE_OUT = fltarr(nlon,nlat,nlev)
;  CLDLIQ_OUT = fltarr(nlon,nlat,nlev)
;  Q_OUT = fltarr(nlon,nlat,nlev)
;  T_OUT = fltarr(nlon,nlat,nlev)
;  VS_OUT = fltarr(nlon,nlat,nlev)
;  US_OUT = fltarr(nlon,nslat,nlev)

  ;set variables
  P0=1.e5
  PS(*) =1.e5
;  TSICE(*)=0.0   ;273.15
;  ICEFRAC(*)=0.0
;  SICTHK(*)=0.0
;  SNOWHICE(*)=0.0


;  iv=nlat
;  for i=0,nlat-1 do begin
;    iv=iv-1
;    TS_OUT(*,i) = mean((TS(*,i)+TS(*,iv))/2.)
;    TS1_OUT(*,i) = mean((TS1(*,i)+TS1(*,iv))/2.)
;    TS2_OUT(*,i) = mean((TS2(*,i)+TS2(*,iv))/2.)
;    TS3_OUT(*,i) = mean((TS3(*,i)+TS3(*,iv))/2.)
;    TS4_OUT(*,i) = mean((TS4(*,i)+TS4(*,iv))/2.)
;    for k=0,nlev-1 do begin
;      CLDICE_OUT(*,i,k) = mean((CLDICE(*,i,k)+CLDICE(*,iv,*))/2.)
;      CLDLIQ_OUT(*,i,k) = mean((CLDLIQ(*,i,k)+CLDLIQ(*,iv,k))/2.)
;      Q_OUT(*,i,k) = mean((Q(*,i,k)+Q(*,iv,k))/2.)
;      T_OUT(*,i,k) = mean((T(*,i,k)+T(*,iv,k))/2.)
;      VS_OUT(*,i,k) = mean((VS(*,i,k)+VS(*,iv,k))/2.)
;    endfor
;   endfor




  ;---- update ncdata file ----
  ncid = ncdf_open(file_ne16np4_icfile_out, /WRITE)
  print, "updatinf ncdata file"
  ncdf_varput, ncid, 'P0',P0
  ncdf_varput, ncid, 'PS',PS
;  ncdf_varput, ncid, 'TSICE',TSICE
;  ncdf_varput, ncid, 'SICTHK',SICTHK
;  ncdf_varput, ncid, 'ICEFRAC',ICEFRAC
;  ncdf_varput, ncid, 'SNOWHICE',SNOWHICE
 ; ncdf_varput, ncid, 'TS',TS_OUT
 ; ncdf_varput, ncid, 'TS1',TS1_OUT
 ; ncdf_varput, ncid, 'TS2',TS2_OUT
 ; ncdf_varput, ncid, 'TS3',TS3_OUT
 ; ncdf_varput, ncid, 'TS4',TS4_OUT
 ; ncdf_varput, ncid, 'CLDICE',CLDICE_OUT
 ; ncdf_varput, ncid, 'CLDLIQ',CLDLIQ_OUT
 ; ncdf_varput, ncid, 'Q',Q_OUT
 ; ncdf_varput, ncid, 'T',T_OUT
 ; ncdf_varput, ncid, 'VS',VS_OUT
 ; ncdf_varput, ncid, 'US',US_OUT
  ncdf_close, ncid




endif





if (make_ne5np4_icfile eq 1) then begin
 ; oprstr="cp -r " + file_ne16np4_icfile + " " + file_ne16np4_icfile_out
 ; spawn, oprstr


  ncid=ncdf_open(file_ne5np4_icfile, /nowrite)
  ncdf_varget,ncid,'lon',lon
  ncdf_varget,ncid,'lat',lat
  ncdf_varget,ncid,'lev',lev
  ncdf_varget,ncid,'P0',P0
  ncdf_varget,ncid,'CLDICE',CLDICE ;[lev,ncol]
  ncdf_varget,ncid,'CLDLIQ',CLDLIQ 
  ncdf_varget,ncid,'ICEFRAC',ICEFRAC ;[ncol]
  ncdf_varget,ncid,'PS',PS ;[ncol]
  ncdf_varget,ncid,'Q',Q 
  ncdf_varget,ncid,'SICTHK',SICTHK 
  ncdf_varget,ncid,'SNOWHICE',SNOWHICE 
  ncdf_varget,ncid,'T',T 
  ncdf_varget,ncid,'TS1',TS1 
  ncdf_varget,ncid,'TS2',TS2 
  ncdf_varget,ncid,'TS3',TS3 
  ncdf_varget,ncid,'TS4',TS4 
  ncdf_varget,ncid,'TSICE',TSICE 
  ncdf_varget,ncid,'U',U 
  ncdf_varget,ncid,'V',V
  ncdf_close,ncid

  nlat=n_elements(lat)
  nlon=n_elements(lon)
  nlev=n_elements(lev)

;fix this
;  PS_OUT = fltarr(nlon,nlat)
;  TSICE_OUT = fltarr(nlon,nlat)
;  SICTHK_OUT = fltarr(nlon,nlat)
;  SNOWHICE_OUT = fltarr(nlon,nlat)
;  ICEFRAC_OUT = fltarr(nlon,nlat)
;  TS_OUT = fltarr(nlon,nlat)
;  TS1_OUT = fltarr(nlon,nlat)
;  TS2_OUT = fltarr(nlon,nlat)
;  TS3_OUT = fltarr(nlon,nlat)
;  TS4_OUT = fltarr(nlon,nlat)
;  CLDICE_OUT = fltarr(nlon,nlat,nlev)
;  CLDLIQ_OUT = fltarr(nlon,nlat,nlev)
;  Q_OUT = fltarr(nlon,nlat,nlev)
;  T_OUT = fltarr(nlon,nlat,nlev)
;  VS_OUT = fltarr(nlon,nlat,nlev)
;  US_OUT = fltarr(nlon,nslat,nlev)

  ;set variables
  P0=1.e5
  PS(*) =1.e5
  TSICE(*)=0.0   ;273.15
  ICEFRAC(*)=0.0
  SICTHK(*)=0.0
  SNOWHICE(*)=0.0


;  iv=nlat
;  for i=0,nlat-1 do begin
;    iv=iv-1
;    TS_OUT(*,i) = mean((TS(*,i)+TS(*,iv))/2.)
;    TS1_OUT(*,i) = mean((TS1(*,i)+TS1(*,iv))/2.)
;    TS2_OUT(*,i) = mean((TS2(*,i)+TS2(*,iv))/2.)
;    TS3_OUT(*,i) = mean((TS3(*,i)+TS3(*,iv))/2.)
;    TS4_OUT(*,i) = mean((TS4(*,i)+TS4(*,iv))/2.)
;    for k=0,nlev-1 do begin
;      CLDICE_OUT(*,i,k) = mean((CLDICE(*,i,k)+CLDICE(*,iv,*))/2.)
;      CLDLIQ_OUT(*,i,k) = mean((CLDLIQ(*,i,k)+CLDLIQ(*,iv,k))/2.)
;      Q_OUT(*,i,k) = mean((Q(*,i,k)+Q(*,iv,k))/2.)
;      T_OUT(*,i,k) = mean((T(*,i,k)+T(*,iv,k))/2.)
;      VS_OUT(*,i,k) = mean((VS(*,i,k)+VS(*,iv,k))/2.)
;    endfor
;   endfor




  ;---- update ncdata file ----
  ncid = ncdf_open(file_ne5np4_icfile_out, /WRITE)
  print, "updatinf ncdata file"
  ncdf_varput, ncid, 'P0',P0
  ncdf_varput, ncid, 'PS',PS
  ncdf_varput, ncid, 'TSICE',TSICE
  ncdf_varput, ncid, 'SICTHK',SICTHK
  ncdf_varput, ncid, 'ICEFRAC',ICEFRAC
  ncdf_varput, ncid, 'SNOWHICE',SNOWHICE
 ; ncdf_varput, ncid, 'TS',TS_OUT
 ; ncdf_varput, ncid, 'TS1',TS1_OUT
 ; ncdf_varput, ncid, 'TS2',TS2_OUT
 ; ncdf_varput, ncid, 'TS3',TS3_OUT
 ; ncdf_varput, ncid, 'TS4',TS4_OUT
 ; ncdf_varput, ncid, 'CLDICE',CLDICE_OUT
 ; ncdf_varput, ncid, 'CLDLIQ',CLDLIQ_OUT
 ; ncdf_varput, ncid, 'Q',Q_OUT
 ; ncdf_varput, ncid, 'T',T_OUT
 ; ncdf_varput, ncid, 'VS',VS_OUT
 ; ncdf_varput, ncid, 'US',US_OUT
  ncdf_close, ncid




endif

end
