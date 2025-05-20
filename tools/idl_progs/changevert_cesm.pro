pro changevert_cesm
;====================================================================
; Author: Eric T. Wolf
; Date:  some time long ago
; needs fixing
;====================================================================
; Change number of vertical levels in a CAM initial condition file.
; This is achieved starting with a 66 level WACCM initial condition
; file and cutting it down to size.  Climate data is imported from 
; separate file.
;
;
; 
;====================================================================

do_write_file = 1
do_print_diagnostic = 1
do_dry = 1
do_snowice = 0  ; requires read of cam.h0 file in addition to .i. file

;================================
;=======  name of new file ======
;================================
fname_out = '/discover/nobackup/etwolf/models/ExoCAM/cesm1.2.1/initial_files/cam_land_fv/ic_0.1bar_L40_zmean_dry_ic.nc'
;fname_out = '/discover/nobackup/etwolf/models/ExoCAM/cesm1.2.1/initial_files/cam_aqua_fv/gj251c_land_h2_L55.nc'
diagnostic_out = fname_out


;L40 --- 3 orders of magnitude pressure
;L46 --- 4 orders of magnitude pressure
;L51 --- 5 orders of magnitude pressure

nlev_out = 40
nilev_out = nlev_out + 1

;read in appropriate grid coordinates from WACCM 66 level grid 
lev_fname_new = '/discover/nobackup/etwolf/models/ExoCAM/cesm1.2.1/initial_files/other/oxygen_CE.cam2.avg.nc'

ncid=ncdf_open(lev_fname_new, /nowrite)
ncdf_varget,ncid,'lev',lev_new
ncdf_varget,ncid,'ilev',ilev_new
ncdf_varget,ncid,'lat',lat_new
ncdf_varget,ncid,'lon',lon_new
ncdf_varget,ncid,'hyai',hyai_new
ncdf_varget,ncid,'hybi',hybi_new
ncdf_varget,ncid,'hyam',hyam_new
ncdf_varget,ncid,'hybm',hybm_new
ncdf_varget,ncid,'P0',P0_new 
ncdf_varget,ncid,'PS',PS_new 
ncdf_close,ncid

nlev_new = n_elements(lev_new)
nilev_new = n_elements(ilev_new)
nlat_new = n_elements(lat_new)
nlon_new = n_elements(lon_new)



;===========================================
;=======  name file with climate data ======
;===========================================
;clim_fname_in ="/discover/nobackup/etwolf/cesm_scratch/archive/ExoCAM_thai_hab2_L51_n68equiv/rest/0197-01-01-00000/ExoCAM_thai_hab2_L51_n68equiv.cam.i.0197-01-01-00000.nc" 
;clim_fname_in = "/discover/nobackup/etwolf/cesm_scratch/rundir/samosa5/run/samosa5.cam.i.0004-01-01-00000.nc"
clim_i_fname_in = "/discover/nobackup/etwolf/models/ExoCAM/cesm1.2.1/initial_files/cam_aqua_fv/ic_1bar_L51_zmean_ic.nc" 
;clim_fname_in = "/discover/nobackup/etwolf/models/ExoCAM/cesm1.2.1/initial_files/cam_aqua_fv/wolf1069_land_planet.i.nc"
;clim_fname_in = "/discover/nobackup/etwolf/models/ExoCAM/cesm1.2.1/initial_files/mars/atm/present_day_mars.cam.i.L40.nc"
;clim_i_fname_in = "/gpfsm/dnb33/etwolf/cesm_scratch/archive/gj251c_hycean/rest/0191-08-01-00000/gj251c_hycean.cam.i.0191-01-01-00000.nc"
;clim_h_fname_in = "/gpfsm/dnb33/etwolf/cesm_scratch/archive/gj251c_hycean/rest/0191-08-01-00000/gj251c_hycean.cam.h0.0191-07.nc"
;clim_i_fname_in = "/gpfsm/dnb33/etwolf/cesm_scratch/archive/gj251c_land_h2/rest/0548-12-01-00000/gj251c_land_h2.cam.i.0548-01-01-00000.nc"
;clim_h_fname_in = "/gpfsm/dnb33/etwolf/cesm_scratch/archive/gj251c_land_h2/rest/0548-12-01-00000/gj251c_land_h2.cam.h0.0548-11.nc"




; read clim_i file (cam.i. initial condition file)
ncid=ncdf_open(clim_i_fname_in, /nowrite)
ncdf_varget,ncid,'lev',lev_clim
ncdf_varget,ncid,'ilev',ilev_clim
ncdf_varget,ncid,'lat',lat
ncdf_varget,ncid,'lon',lon
ncdf_varget,ncid,'hyai',hyai_clim
ncdf_varget,ncid,'hybi',hybi_clim
ncdf_varget,ncid,'hyam',hyam_clim
ncdf_varget,ncid,'hybm',hybm_clim

nlev_clim = n_elements(lev_clim)
nilev_clim = n_elements(ilev_clim)
nlat = n_elements(lat)
nlon = n_elements(lon)
;kludge, I am keeping with a 4x5 grid
nslat=nlat-1
nslon=nlon

ncdf_varget,ncid,'CLDICE',CLDICE_clim
ncdf_varget,ncid,'CLDLIQ',CLDLIQ_clim
ncdf_varget,ncid,'Q',Q_clim
ncdf_varget,ncid,'T',T_clim
ncdf_varget,ncid,'US',US_clim
ncdf_varget,ncid,'VS',VS_clim
ncdf_varget,ncid,'PS',PS_clim
ncdf_varget,ncid,'P0',P0_clim

;stuff that doesn't need to be changed
ncdf_varget,ncid,'slat',slat
ncdf_varget,ncid,'slon',slon
ncdf_varget,ncid,'w_stag',w_stag
ncdf_varget,ncid,'time',time
ncdf_varget,ncid,'time_bnds',time_bnds
ncdf_varget,ncid,'date_written',date_written
ncdf_varget,ncid,'time_written',time_written
ncdf_varget,ncid,'ntrm',ntrm
ncdf_varget,ncid,'ntrn',ntrn
ncdf_varget,ncid,'ntrk',ntrk
ncdf_varget,ncid,'ndbase',ndbase
ncdf_varget,ncid,'nsbase',nsbase
ncdf_varget,ncid,'nbdate',nbdate
ncdf_varget,ncid,'nbsec',nbsec
ncdf_varget,ncid,'mdt',mdt
ncdf_varget,ncid,'gw',gw
ncdf_varget,ncid,'ndcur',ndcur
ncdf_varget,ncid,'nscur',nscur
ncdf_varget,ncid,'date',date
ncdf_varget,ncid,'datesec',datesec
ncdf_varget,ncid,'nsteph',nsteph
; these are all zeros, why?
ncdf_varget,ncid,'TS1',TS1
ncdf_varget,ncid,'TS2',TS2
ncdf_varget,ncid,'TS3',TS3
ncdf_varget,ncid,'TS4',TS4
ncdf_varget,ncid,'TSICE',TSICE
ncdf_varget,ncid,'SICTHK',SICTHK_clim
ncdf_close,ncid

; read clim_h (cam.h0) history file
if (do_snowice eq 1) then begin
  ncid=ncdf_open(clim_h_fname_in, /nowrite)
  ncdf_varget,ncid,'TS',TS_clim
  ncdf_varget,ncid,'SNOWHICE',SNOWHICE_clim
  ncdf_varget,ncid,'ICEFRAC',ICEFRAC_clim
  ncdf_close,ncid
endif 

  
; create pressure grids from hybrid sigma coordinates
lev_P_new=fltarr(nlon, nlat, nlev_new)    ;[Pa] pressure coordinate matrix, layer midpoints    
ilev_P_new=fltarr(nlon, nlat, nilev_new)   ;[Pa] pressure coordinate matrix, layer interfaces  
hybrid2pressure,nlon,nlat,nlev_new,PS_new,P0_new,hyam_new,hybm_new,hyai_new,hybi_new,lev_P_new,ilev_P_new

lev_P_clim=fltarr(nlon, nlat, nlev_clim)    ;[Pa] pressure coordinate matrix, layer midpoints     
ilev_P_clim=fltarr(nlon, nlat, nilev_clim)   ;[Pa] pressure coordinate matrix, layer interfaces 
hybrid2pressure,nlon,nlat,nlev_clim,PS_clim,P0_clim,hyam_clim,hybm_clim,hyai_clim,hybi_clim,lev_P_clim,ilev_P_clim

;out variables
CLDICE_out=fltarr(nlon,nlat,nlev_out)
CLDLIQ_out=fltarr(nlon,nlat,nlev_out)
CLOUD_out=fltarr(nlon,nlat,nlev_out)
Q_out=fltarr(nlon,nlat,nlev_out)
T_out=fltarr(nlon,nlat,nlev_out)
US_out=fltarr(nlon,nslat,nlev_out)
VS_out=fltarr(nslon,nlat,nlev_out)
PS_out=fltarr(nlon,nlat)
SICTHK_out=fltarr(nlon,nlat)
ICEFRAC_out=fltarr(nlon,nlat)
SNOWHICE_out=fltarr(nlon,nlat)


lev_out=fltarr(nlev_out)
ilev_out=fltarr(nilev_out)
hyai_out=fltarr(nilev_out)
hybi_out=fltarr(nilev_out)
hyam_out=fltarr(nlev_out)
hybm_out=fltarr(nlev_out)



n=66-nlev_out
for x=0,nlon-1 do begin
  for y=0,nlat-1 do begin
    T_out(x,y,*) = interpol(T_clim(x,y,*), lev_P_clim(x,y,*), lev_P_new(x,y,n:65))  
    CLDICE_out(x,y,*) = interpol(CLDICE_clim(x,y,*), lev_P_clim(x,y,*), lev_P_new(x,y,n:65))  
    CLDLIQ_out(x,y,*) = interpol(CLDLIQ_clim(x,y,*), lev_P_clim(x,y,*), lev_P_new(x,y,n:65))  
    Q_out(x,y,*)      = interpol(Q_clim(x,y,*), lev_P_clim(x,y,*), lev_P_new(x,y,n:65))  
 ;   VS_out(x,y,*)     = interpol(VS_clim(x,y,*), lev_P_clim(x,y,*), lev_P_new(x,y,n:65))  
    VS_out(x,y,*) = 0.0

    if (do_dry eq 1) then begin
      CLDICE_out(x,y,*) = 0.0
      CLDLIQ_out(x,y,*) = 0.0
      Q_out(x,y,*) = 0.0
    endif

    PS_out(x,y)       = PS_clim(x,y)
    SICTHK_out(x,y)   = 5.0  ;SICTHK_clim(x,y)
    if (do_snowice eq 1) then begin
      ICEFRAC_out(x,y)  = ICEFRAC_clim(x,y)
      SNOWHICE_out(x,y) = SNOWHICE_clim(x,y)
    endif else begin
      ICEFRAC_out(x,y)  = 0.0
      SNOWHICE_out(x,y) = 0.0
    endelse  
  endfor
endfor
P0_out = P0_clim

for x=0,nlon-1 do begin
  for y=0,nlat-2 do begin
  ;  US_out(x,y,*) = interpol(US_old(x,y,*), lev_P_old(x,y,*), lev_P_new(x,y,n:65))  
    US_out(x,y,*) = 0.0
  endfor
endfor


for i=nilev_out-1,0,-1 do begin
  j=nilev_new-nilev_out+i
  hyai_out(i)        = hyai_new(j)
  hybi_out(i)        = hybi_new(j)
  ilev_out(i)        = ilev_new(j)
endfor

for i=0,nlev_out-1 do begin
  j=nlev_new-nlev_out+i
  hyam_out(i)        = hyam_new(j)
  hybm_out(i)        = hybm_new(j)
  lev_out(i)         = lev_new(j)
;  CLDICE_out(*,*,i)  = CLDICE_in(*,*,j)
;  CLDLIQ_out(*,*,i)  = CLDLIQ_in(*,*,j)
;  Q_out(*,*,i)       = Q_in(*,*,j)
;  T_out(*,*,i)       = T_in(*,*,j)
;  US_out(*,*,i)      = US_in(*,*,j)
;  VS_out(*,*,i)      = VS_in(*,*,j)
endfor

;write

if (do_write_file eq 1) then begin
 spawn, "rm -r -f fname_out"
 print, "creating file ...." 
 print, fname_out
  id = ncdf_create(fname_out,/clobber)
  dim1 = NCDF_DIMDEF(id, 'lat', nlat)
  dim2 = NCDF_DIMDEF(id, 'lon', nlon)  
  dim3 = NCDF_DIMDEF(id, 'slat', nslat)  
  dim4 = NCDF_DIMDEF(id, 'slon', nlon)  
  dim5 = NCDF_DIMDEF(id, 'lev', nlev_out)  
  dim6 = NCDF_DIMDEF(id, 'ilev', nilev_out)  
  dim10 = NCDF_DIMDEF(id, 'time', 1)
  dim11 = NCDF_DIMDEF(id, 'nbnd', 2)
  dim12 = NCDF_DIMDEF(id, 'chars', 8)
    
  varid1 = NCDF_VARDEF(id,'P0',/double)
  varid2 = NCDF_VARDEF(id,'lat',dim1,/double)
  varid3 = NCDF_VARDEF(id,'lon',dim2,/double)

  varid4 = NCDF_VARDEF(id,'slat',dim3,/double)
  varid5 = NCDF_VARDEF(id,'slon',dim4,/double)
  varid6 = NCDF_VARDEF(id,'w_stag',dim3,/double)
  varid7 = NCDF_VARDEF(id,'lev',dim5,/double)
  varid8 = NCDF_VARDEF(id,'ilev',dim6,/double)
  varid12 = NCDF_VARDEF(id,'time',dim10,/double)
  varid13 = NCDF_VARDEF(id,'time_bnds',[dim10,dim1],/double)
  varid14 = NCDF_VARDEF(id,'date_written',[dim10,dim12],/char)
  varid15 = NCDF_VARDEF(id,'time_written',[dim10,dim12],/char)
  varid16 = NCDF_VARDEF(id,'ntrm')
  varid17 = NCDF_VARDEF(id,'ntrn')
  varid18 = NCDF_VARDEF(id,'ntrk')
  varid19 = NCDF_VARDEF(id,'ndbase')
  varid20 = NCDF_VARDEF(id,'nsbase')
  varid21 = NCDF_VARDEF(id,'nbdate')
  varid22 = NCDF_VARDEF(id,'nbsec')
  varid23 = NCDF_VARDEF(id,'mdt')
  varid24 = NCDF_VARDEF(id,'hyai',dim6,/double)
  varid25 = NCDF_VARDEF(id,'hybi',dim6,/double)
  varid26 = NCDF_VARDEF(id,'hyam',dim5,/double)
  varid27 = NCDF_VARDEF(id,'hybm',dim5,/double)
  varid28 = NCDF_VARDEF(id,'gw',dim1,/double)
  varid29 = NCDF_VARDEF(id,'ndcur',dim10)
  varid30 = NCDF_VARDEF(id,'nscur',dim10)
  varid31 = NCDF_VARDEF(id,'date',dim10)
  varid33 = NCDF_VARDEF(id,'datesec',dim10)
  varid34 = NCDF_VARDEF(id,'nsteph',dim10)
  varid57 = NCDF_VARDEF(id,'CLDICE',[dim2,dim1,dim5,dim10],/double)
  varid58 = NCDF_VARDEF(id,'CLDLIQ',[dim2,dim1,dim5,dim10],/double)
  varid76 = NCDF_VARDEF(id,'ICEFRAC',[dim2,dim1,dim10],/double)
  varid99 = NCDF_VARDEF(id,'PS',[dim2,dim1,dim10],/double)
  varid100 = NCDF_VARDEF(id,'Q',[dim2,dim1,dim5,dim10],/double)
  varid103 = NCDF_VARDEF(id,'SICTHK',[dim2,dim1,dim10],/double)
  varid104 = NCDF_VARDEF(id,'SNOWHICE',[dim2,dim1,dim10],/double)
  varid105 = NCDF_VARDEF(id,'T',[dim2,dim1,dim5,dim10],/double)
  varid110 = NCDF_VARDEF(id,'TS1',[dim2,dim1,dim10],/double)
  varid111 = NCDF_VARDEF(id,'TS2',[dim2,dim1,dim10],/double)
  varid112 = NCDF_VARDEF(id,'TS3',[dim2,dim1,dim10],/double)
  varid113 = NCDF_VARDEF(id,'TS4',[dim2,dim1,dim10],/double)
  varid114 = NCDF_VARDEF(id,'TSICE',[dim2,dim1,dim10],/double)
  varid117 = NCDF_VARDEF(id,'US',[dim2,dim3,dim5,dim10],/double)
  varid119 = NCDF_VARDEF(id,'VS',[dim4,dim1,dim5,dim10],/double)

  ;attributes
  NCDF_ATTPUT, id, varid1, "long_name", "reference pressure"
  NCDF_ATTPUT, id, varid1, "units", "Pa"
  NCDF_ATTPUT, id, varid2, "long_name", "latitude"
  NCDF_ATTPUT, id, varid2, "units", "degrees_north"
  NCDF_ATTPUT, id, varid3, "long_name", "longitude"
  NCDF_ATTPUT, id, varid3, "units", "degrees_east"
  NCDF_ATTPUT, id, varid4, "long_name", "staggered latitude"
  NCDF_ATTPUT, id, varid4, "units", "degrees_north"
  NCDF_ATTPUT, id, varid5, "long_name", "staggered longitude"
  NCDF_ATTPUT, id, varid5, "units", "degrees_east"
  NCDF_ATTPUT, id, varid6, "long_name", "staggered latitude weights"
  NCDF_ATTPUT, id, varid7, "long_name", "hybrid level at midpoints (1000*(A+B))"
  NCDF_ATTPUT, id, varid7, "units", "level"
  NCDF_ATTPUT, id, varid7, "positive", "down"
  NCDF_ATTPUT, id, varid7, "standard_name", "atmosphere_hybrid_sigma_pressure_coordinate"
  NCDF_ATTPUT, id, varid7, "formula", "a: hyam b: hybm p0: P0 ps: PS"
  NCDF_ATTPUT, id, varid8, "long_name", "hybrid level at interfaces (1000*(A+B))"
  NCDF_ATTPUT, id, varid8, "units", "level"
  NCDF_ATTPUT, id, varid8, "positive", "down"
  NCDF_ATTPUT, id, varid8, "standard_name", "atmosphere_hybrid_sigma_pressure_coordinate"
  NCDF_ATTPUT, id, varid8, "formula", "a: hyai b: hybi p0: P0 ps: PS"
  NCDF_ATTPUT, id, varid12, "long_name", "time"
  NCDF_ATTPUT, id, varid12, "units", "days since 1995-01-01 00:00:00"
  NCDF_ATTPUT, id, varid12, "calendar", "noleap"
  NCDF_ATTPUT, id, varid12, "bounds", "time_bnds"
  NCDF_ATTPUT, id, varid13, "long_name", "time interval endpoints"
  NCDF_ATTPUT, id, varid16, "long_name", "spectral truncation parameter M"
  NCDF_ATTPUT, id, varid17, "long_name", "spectral truncation parameter N"
  NCDF_ATTPUT, id, varid18, "long_name", "spectral truncation parameter K"
  NCDF_ATTPUT, id, varid19, "long_name", "base day"
  NCDF_ATTPUT, id, varid20, "long_name", "seconds of base day"
  NCDF_ATTPUT, id, varid21, "long_name", "base date (YYYYMMDD)"
  NCDF_ATTPUT, id, varid22, "long_name", "seconds of base date"
  NCDF_ATTPUT, id, varid23, "long_name", "timestep"
  NCDF_ATTPUT, id, varid23, "units", "s"
  NCDF_ATTPUT, id, varid24, "long_name", "hybrid A coefficient at layer interfaces"
  NCDF_ATTPUT, id, varid25, "long_name", "hybrid B coefficient at layer interfaces"
  NCDF_ATTPUT, id, varid26, "long_name", "hybrid A coefficient at layer midpoints"
  NCDF_ATTPUT, id, varid27, "long_name", "hybrid B coefficient at layer midpoints"
  NCDF_ATTPUT, id, varid28, "long_name", "gauss weights"
  NCDF_ATTPUT, id, varid29, "long_name", "current day (from base day)"
  NCDF_ATTPUT, id, varid30, "long_name", "current seconds of current day"
  NCDF_ATTPUT, id, varid31, "long_name", "current date (YYYYMMDD)"
  NCDF_ATTPUT, id, varid33, "long_name", "current seconds of current date"
  NCDF_ATTPUT, id, varid34, "long_name", "current timestep"
  NCDF_ATTPUT, id, varid57, "units", "kg/kg"
  NCDF_ATTPUT, id, varid57, "long_name", "Grid box averaged ice condensate amount"
  NCDF_ATTPUT, id, varid58, "units", "kg/kg"
  NCDF_ATTPUT, id, varid58, "long_name", "Grid box averaged liquid condensate amount"
  NCDF_ATTPUT, id, varid76, "units", "fraction"
  NCDF_ATTPUT, id, varid76, "long_name", "Fraction of sfc area covered by sea-ice"
  NCDF_ATTPUT, id, varid99, "units", "Pa"
  NCDF_ATTPUT, id, varid99, "long_name", "surface pressure"
  NCDF_ATTPUT, id, varid100, "units", "kg/kg"
  NCDF_ATTPUT, id, varid100, "long_name", "Specific humidity"
  NCDF_ATTPUT, id, varid103, "units", "m"
  NCDF_ATTPUT, id, varid103, "long_name", "Sea ice thickness"
  NCDF_ATTPUT, id, varid104, "units", "m"
  NCDF_ATTPUT, id, varid104, "long_name", "Water equivalent snow depth"
  NCDF_ATTPUT, id, varid105, "units", "K"
  NCDF_ATTPUT, id, varid105, "long_name", "Temperature"
  NCDF_ATTPUT, id, varid110, "units", "K"
  NCDF_ATTPUT, id, varid110, "long_name", "TS1 subsoil temperature"
  NCDF_ATTPUT, id, varid111, "units", "K"
  NCDF_ATTPUT, id, varid111, "long_name", "TS2 subsoil temperature"
  NCDF_ATTPUT, id, varid112, "units", "K"
  NCDF_ATTPUT, id, varid112, "long_name", "TS3 subsoil temperature"
  NCDF_ATTPUT, id, varid113, "units", "K"
  NCDF_ATTPUT, id, varid113, "long_name", "TS4 subsoil temperature"
  NCDF_ATTPUT, id, varid114, "units", "K"
  NCDF_ATTPUT, id, varid114, "long_name", "Ice temperature"
  NCDF_ATTPUT, id, varid117, "units", "m/s"
  NCDF_ATTPUT, id, varid117, "long_name", "Zonal wind staggered"
  NCDF_ATTPUT, id, varid119, "units", "m/s"
  NCDF_ATTPUT, id, varid119, "long_name", "Meridional wind staggered"

  ncdf_control, id, /ENDEF

  
  NCDF_VARPUT, id, varid1, P0_out
  NCDF_VARPUT, id, varid2, lat
  NCDF_VARPUT, id, varid3, lon
  NCDF_VARPUT, id, varid4, slat
  NCDF_VARPUT, id, varid5, slon
  NCDF_VARPUT, id, varid6, w_stag
  NCDF_VARPUT, id, varid7, lev_out
  NCDF_VARPUT, id, varid8, ilev_out


  NCDF_VARPUT, id, varid12, time
  NCDF_VARPUT, id, varid13, time_bnds
  NCDF_VARPUT, id, varid14, date_written
  NCDF_VARPUT, id, varid15, time_written
  NCDF_VARPUT, id, varid16, ntrm
  NCDF_VARPUT, id, varid17, ntrn
  NCDF_VARPUT, id, varid18, ntrk
  NCDF_VARPUT, id, varid19, ndbase
  NCDF_VARPUT, id, varid20, nsbase
  NCDF_VARPUT, id, varid21, nbdate
  NCDF_VARPUT, id, varid22, nbsec
  NCDF_VARPUT, id, varid23, mdt
  NCDF_VARPUT, id, varid24, hyai_out
  NCDF_VARPUT, id, varid25, hybi_out
  NCDF_VARPUT, id, varid26, hyam_out
  NCDF_VARPUT, id, varid27, hybm_out
  NCDF_VARPUT, id, varid28, gw
  NCDF_VARPUT, id, varid29, ndcur
  NCDF_VARPUT, id, varid30, nscur
  NCDF_VARPUT, id, varid31, date
  NCDF_VARPUT, id, varid33, datesec
  NCDF_VARPUT, id, varid34, nsteph


  NCDF_VARPUT, id, varid57, CLDICE_out
  NCDF_VARPUT, id, varid58, CLDLIQ_out
  NCDF_VARPUT, id, varid76, ICEFRAC_out
  NCDF_VARPUT, id, varid99, PS_out
  NCDF_VARPUT, id, varid100, Q_out
  NCDF_VARPUT, id, varid103, SICTHK_out
  NCDF_VARPUT, id, varid104, SNOWHICE_out
  NCDF_VARPUT, id, varid105, T_out
  NCDF_VARPUT, id, varid110, TS1
  NCDF_VARPUT, id, varid111, TS2
  NCDF_VARPUT, id, varid112, TS3
  NCDF_VARPUT, id, varid113, TS4
  NCDF_VARPUT, id, varid114, TSICE
  
  NCDF_VARPUT, id, varid117, US_out
  NCDF_VARPUT, id, varid119, VS_out

  ncdf_close, id

endif

if (do_print_diagnostic eq 1) then begin

  ncid=ncdf_open(diagnostic_out, /nowrite)
  ncdf_varget,ncid,'lev',lev
  ncdf_varget,ncid,'ilev',ilev
  ncdf_varget,ncid,'lat',lat
  ncdf_varget,ncid,'lon',lon
  ncdf_varget,ncid,'hyai',hyai
  ncdf_varget,ncid,'hybi',hybi
  ncdf_varget,ncid,'hyam',hyam
  ncdf_varget,ncid,'hybm',hybm
  ncdf_varget,ncid,'PS',PS
  ncdf_varget,ncid,'P0',P0
  ncdf_varget,ncid,'T',T
  ncdf_varget,ncid,'Q',Q

  nlon=n_elements(lon)
  nlat=n_elements(lat)
  nlev=n_elements(lev)
  nilev=n_elements(ilev)

  lev_P=fltarr(nlon,nlat,nlev)
  ilev_P=fltarr(nlon,nlat,nilev)
  hybrid2pressure,nlon,nlat,nlev,PS,P0,hyam,hybm,hyai,hybi,lev_P,ilev_P
  lev_P(*,*,*) = lev_P(*,*,*)/100.0

  PLEV = fltarr(nlev)
  TLEV = fltarr(nlev)
  QLEV = fltarr(nlev)
  print, "PavgZ, TavgZ, QavgZ"
  for z=0, nlev-1 do begin
     area_weighted_avg_gen, lon, lat, lev_P(*,*,z), PavgZ        &   PLEV(z) = PavgZ
     area_weighted_avg_gen, lon, lat, T(*,*,z), TavgZ            &   TLEV(z) = TavgZ
     area_weighted_avg_gen, lon, lat, Q(*,*,z), QavgZ            &   QLEV(z) = QavgZ
     print, z, PavgZ, TavgZ, QavgZ
  endfor

endif



end
