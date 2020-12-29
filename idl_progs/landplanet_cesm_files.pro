pro landplanet_cesm_files
;-------------------------------


make_bndtopo = 0 ;1
make_lnd_domain = 0 
make_ncdata = 1
make_finidat = 0


;--- 4x5 ---  from my simulations
file_bndtopo_in = '/Users/wolfe/Desktop/NASA_Planetary_Atmospheres/Mstar_Kopparapu2016/aqua_4x5/USGS-gtopo_aquaplanet_4x5.nc'
file_bndtopo_out = '/Users/wolfe/Desktop/landplanet/topo_landplanet.nc'

;========== BND_TOPO =====================
if (make_bndtopo eq 1) then begin

  oprstr="cp -r " + file_bndtopo_in + " " + file_bndtopo_out
  spawn, oprstr
  
  ; load metrics
  ncid=ncdf_open(file_bndtopo_in, /nowrite)
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
  LANDFRAC(*,*) = 1.0
  LANDM_COSLAT(*,*) = 0.0


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


file_lnd_domain_in =  '/Users/wolfe/Desktop/landplanet/domain.lnd.fv4x5_gx3v7.091218.nc'
file_lnd_domain_out = '/Users/wolfe/Desktop/landplanet/domain.lnd.fv4x5_lanplanet.nc'

if (make_lnd_domain eq 1) then begin

  oprstr="cp -r " + file_lnd_domain_in + " " + file_lnd_domain_out
  spawn, oprstr
  
  ; load metrics
  ncid=ncdf_open(file_lnd_domain_in, /nowrite)
  ncdf_varget,ncid,'mask',mask
  ncdf_varget,ncid,'frac',frac
  ncdf_close,ncid

  ;set everything to one
  mask(*,*) = 1.0
  frac(*,*) = 1.0

  ;---- update lnd domain file ----
  ncid = ncdf_open(file_lnd_domain_out, /WRITE)
  print, "updating lnd_domain"
  print, file_lnd_domain_out
  ncdf_varput, ncid, 'mask',mask
  ncdf_varput, ncid, 'frac',frac
  ncdf_close, ncid
endif


;--- 4x5 ---  from my simulations
file_ncdata_in = '/Users/wolfe/Desktop/landplanet/initial_files/control_L40.cam.i.0009-01-01-00000.nc'
;file_ncdata_out = '/Users/wolfe/Desktop/landplanet/initial_files/control_L40_dry.cam.i.nc'
file_ncdata_out = '/Users/wolfe/Desktop/landplanet/initial_files/control_L40_dry.cam.i.nc'
;========== ncdata =====================
if (make_ncdata eq 1) then begin

  oprstr="cp -r " + file_ncdata_in + " " + file_ncdata_out
  spawn, oprstr
  
  ; load metrics
  ncid=ncdf_open(file_ncdata_in, /nowrite)
  ncdf_varget,ncid,'lon',lon
  ncdf_varget,ncid,'lat',lat
  ncdf_varget,ncid,'CLDICE',CLDICE
  ncdf_varget,ncid,'CLDLIQ',CLDLIQ
  ncdf_varget,ncid,'Q',Q
  ncdf_varget,ncid,'T',T
  ncdf_varget,ncid,'PS',PS
  nlon=n_elements(lon)
  nlat=n_elements(lat)

  ;T(*,*) = 400.
  PS(*,*) = 100000.
  ;set everything to zero FV
  CLDICE(*,*,*) = 0.0
  CLDLIQ(*,*,*) = 0.0
  Q(*,*,*) = 0.0

  ;---- update bndtopo file ----
  ncid = ncdf_open(file_ncdata_out, /WRITE)
  print, "updating ncdata"
  print, file_ncdata_out
  ncdf_varput, ncid, 'CLDICE',CLDICE
  ncdf_varput, ncid, 'CLDLIQ',CLDLIQ
  ncdf_varput, ncid, 'Q',Q
  ncdf_varput, ncid, 'T',T
  ncdf_varput, ncid, 'PS',PS
  ncdf_close, ncid

endif



file_finidat_in = '/Users/wolfe/Desktop/landplanet/initial_files/landplanet_q0_wt10.clm2.r.nc'
file_finidat_out = '/Users/wolfe/Desktop/landplanet/initial_files/landplanet_small.clm2.r.nc'
;========== finidat =====================
if (make_finidat eq 1) then begin

  oprstr="cp -r " + file_finidat_in + " " + file_finidat_out
  spawn, oprstr
  
  ; load metrics
  ncid=ncdf_open(file_finidat_in, /nowrite)
  ncdf_varget,ncid,'SNOWDP',SNOWDP &   SNOWDP(*) = 0.0   ; snow depth
  ncdf_varget,ncid,'SNLSNO',SNLSNO &   SNLSNO(*) = 0.0   ; number of snow layers
  ncdf_varget,ncid,'WA',WA  & WA(*) = 0.0                ; water in confined aquifer
  ncdf_varget,ncid,'WT',WT  & WT(*) = 0.0                ; total water storage
  ncdf_varget,ncid,'ZWT',ZWT &  ZWT(*) = 0.0             ; water table depth
  ncdf_varget,ncid,'frac_sno',frac_sno & frac_sno(*) = 0.0 ; fraaction of ground covered by snow 
  ncdf_varget,ncid,'DZSNO',DZSNO & DZSNO(*,*) = 0.0      ; snow layer thickness
  ncdf_varget,ncid,'ZSNO',ZSNO  & ZSNO(*,*) = 0.0        ; snow layer depth
  ncdf_varget,ncid,'ZISNO',ZISNO  & ZISNO(*,*) = 0.0     ; snow interface depth
  ncdf_varget,ncid,'H2OSNO',H2OSNO  & H2OSNO(*) = 0.0    ; snow water
  ncdf_varget,ncid,'H2OSOI_LIQ',H2OSOI_LIQ &  H2OSOI_LIQ(*,*) = 1.0e-20   ; soil liquid water  [kg m-2]
  ncdf_varget,ncid,'H2OSOI_ICE',H2OSOI_ICE &  H2OSOI_ICE(*,*) = 1.0e-20   ; ice lquid water
  ; there are a number of these bullshit variables with some non-zero 
  ; value for where snow exists
  ;ncdf_varget,ncid,'snw_rds',snw_rds & snw_rds(*,*) = 0.0
  ;ncdf_varget,ncid,'mss_bcpho', mss_bcpho & mss_bcpho(*,*) = 0.0
  ;ncdf_varget,ncid,'mss_bcphi', mss_bcphi & mss_bcphi(*,*) = 0.0
  ;ncdf_varget,ncid,'mss_ocpho', mss_ocpho & mss_ocpho(*,*) = 0.0
  ;ncdf_varget,ncid,'mss_ocphi', mss_ocphi & mss_ocphi(*,*) = 0.0

  nlon=n_elements(lon)
  nlat=n_elements(lat)

  ;---- update bndtopo file ----
  ncid = ncdf_open(file_finidat_out, /WRITE)
  print, "updating finidat"
  print, file_finidat_out
  ncdf_varput,ncid,'SNOWDP',SNOWDP   ; snow depth
  ncdf_varput,ncid,'SNLSNO',SNLSNO   ; number of snow layers
  ncdf_varput,ncid,'WA',WA           ; water in confined aquifer
  ncdf_varput,ncid,'WT',WT           ; total water storage
  ncdf_varput,ncid,'ZWT',ZWT         ; water table depth
  ncdf_varput,ncid,'frac_sno',frac_sno  ; fraaction of ground covered by snow 
  ncdf_varput,ncid,'DZSNO',DZSNO       ; snow layer thickness
  ncdf_varput,ncid,'ZSNO',ZSNO          ; snow layer depth
  ncdf_varput,ncid,'ZISNO',ZISNO       ; snow interface depth
  ncdf_varput,ncid,'H2OSNO',H2OSNO      ; snow water
  ncdf_varput,ncid,'H2OSOI_LIQ',H2OSOI_LIQ   ; soil liquid water
  ncdf_varput,ncid,'H2OSOI_ICE',H2OSOI_ICE   ; ice lquid water
  ncdf_close, ncid

endif




end
