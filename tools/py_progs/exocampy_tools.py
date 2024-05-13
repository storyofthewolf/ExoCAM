import numpy as np
import math


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# // calculate pressure coordinates from hybrid-sigma levels
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def hybrid2pressure(nlon,nlat,nlev,PS,P0,hyam,hybm,hyai,hybi):
  """
  AUTHOR WOLF E.T.
  7/11/2008
  Translated to Python by R. Deitrick, October 2022
  -------------------------------------------------------------
  ARGUMENTS:
  nlon                => number of elements, longitude array
  nlat                => number of elements, latitude array
  nlev                => number of elements, hybrid vertical level midpoint coordinates
  PS[lon,lat]         => surface pressure array at each grid point
  P0                  => reference pressre
  hyam                => hybrid A coefficient at layer midpoints
  hybm                => hybrid B coefficient at layer midpoints
  hyai                => hybrid A coefficient at layer interfaces
  hybi                => hybrid B coefficient at layer interfaces
  lev_P[lon,lat,lev]  => returns mid layer pressure coordinate matrix, must be
                         defined with proper dimensions before calling this procedure
  ilev_P[lon,lat,lev+1] => returns interface pressure coordinate matrix, must be
                         defined with proper dimensions before calling this procedure
  ---------------------------------------------------------------
  """

  lev_P = hyam[:,None,None]*P0 + hybm[:,None,None]*PS[None,:,:]
  ilev_P = hyai[:,None,None]*P0 + hybi[:,None,None]*PS[None,:,:]

  return (lev_P, ilev_P)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# // calculate height coordinates
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def hybrid2height(nlon,nlat,nlev,PS,P0,hyam,hybm,hyai,hybi,T,G,R):
  """
  ;AUTHOR: WOLF, E.T.
  ;7/14/2008
  ; NOTES:
  ; Needs incorporation of surface geopotential heigh from topography data
  ; Currently surface is set to zero.
  ;---------------------------------------------------------------------
  ;ARGUMENTS:
  ;nlon                => number of elements,  longitude array
  ;nlat                => numbrt of elements, latitude array
  ;nilev               => number of elements, hybrid vertical levels at interfaces
  ;PS[lon,lat]         => surface pressure array at each grid point
  ;P0                  => reference pressure
  ;hyai                => hybrid A coefficient at layer interfaces
  ;hybi                => hybrid B coefficient at layer interfaces
  ;hyam                => hybrid A coefficient at midlayer
  ;hybm                => hybrid B coefficient at midlayer
  ;T                   => temperature field (x,y,z)
  ;G                   => gravity of the planet ;[m s-2] acceleration due to gravity
  ;R                   => gas constant [J kg-1 K-1]   gas constant for dry air
  ;lev_Z[lon,lat,lev]  => returns mid layer height coordinate matrix in METERS,
  ;                       must be defined with proper dimensions before being passed
  ;ilev_Z[lon,lat,nilev]  => returns interface height coordinate matrix in METERS,
  ;                       must be defined with proper dimensions before being passed
  ;----------------------------------------------------------------------                      
  """

  # geopotential height of surace
  # currently set to zero
  # will need to make this an input argument to account
  # for topography
  PHIS   = np.zeros((nlon, nlat), dtype=float)
  z_surf = np.zeros((nlon, nlat), dtype=float)
  
  # first calculate pressures
  lev_P = hyam[:,None,None]*P0 + hybm[:,None,None]*PS[None,:,:]
  ilev_P = hyai[:,None,None]*P0 + hybi[:,None,None]*PS[None,:,:]

  # define dimensions from Temperature array
  nlev = T.shape[0]
  nlat = T.shape[1]
  nlon = T.shape[2]
  nilev = nlev+1


  ilev_Z = np.zeros((nilev, nlat, nlon), dtype=float) 
  lev_Z  = np.zeros((nlev, nlat, nlon), dtype=float) 
  
  z_surf[:,:]=PHIS[:,:]/G

  for z in range(nlev,0,-1): 
    for y in range(nlat):
      for x in range(nlon):
        zi=z-1
        ptop = ilev_P[zi,y,x]
        pbot = ilev_P[zi+1,y,x]
        pmid = lev_P[zi,y,x] 
        delta_Z = R*T[zi,y,x]/G*math.log(pbot/ptop)
        ilev_Z[zi,y,x]=ilev_Z[zi+1,y,x]+delta_Z
        Z_SCALE=-R*T[zi,y,x]/G*math.log(pmid/pbot)
        lev_Z[zi,y,x]=ilev_Z[zi+1,y,x]+Z_SCALE

  return (lev_Z, ilev_Z)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# // calculate area weighted averaged of a longitude-latitude field
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def area_weighted_avg(lon, lat, var):
  """
  AUTHOR: WOLF E.T.
  Cinco de Mayo, 2023
  Thanks ChatGPT!
  -------------------------------------------------------------
  PURPOSE: Calculate area weighted average of geophysical quantities 
           from any 2D ExoCAM output data.  Can originate from atmosphere,
           ice, or land models.  Code is able to accept variable arrays
           in either [lon,lat], or [lat,lon] formats

  NOTES: Error handling: program expects fill_value of -999.0.  
         Multiple pole counting for cartesian grids removed.
  -------------------------------------------------------------
  Arguments
  lon            =>  longitude array from CAM 
  lat            =>  lattitude array from CAM
  var            =>  atmospheric variable we wish to calculate the 
                     area weighted averaged of
                     Python input args are dynamically alocated 
                     at runtime. 

  weighted_avg   =>  returns the area weighted average
  -------------------------------------------------------------
  """

  # Define the fill value
  fill_value = -999.0

  # Determine the number of rows and columns in variable array
  n0 = var.shape[0]
  n1 = var.shape[1]

  # Determine the number of elements in lat and lon arrays
  nlat = len(lat)  
  nlon = len(lon)

  # Determine which dimension is latitude, which is longitude
  if (nlon == n0):
    nl0=nlon
    nl1=nlat
  else:
    nl0=nlat
    nl1=nlon


  #print(n0,n1,nlon,nlat,nl0,nl1)

  # Set summing quantities to zero
  weighted_avg = 0.0
  sumArea = 0.0
  missing_area = 0.0

  # Define some arrays
  area = np.zeros((n0, n1), dtype=float)
  slon = np.zeros(nlon + 1, dtype=float)
  slat = np.zeros(nlat + 1, dtype=float)

  # Define staggered grid
  slon[:nlon] = lon[:]
  slon[nlon] = 360.0
  slat[0] = -90.0
  slat[nlat] = 90.0

  # Create staggered latitudes
  for ya in range(1, nlat):
    slat[ya] = (lat[ya - 1] + lat[ya]) / 2.0

  for x in range(nl0):
    for y in range(nl1):
      if (nlon == n0): 
        p=x
        q=y
      else:
        p=y
        q=x
      area[x,y] =  np.pi/180 * (slon[p+1] - slon[p]) * (np.sin(slat[q+1] * np.pi/180) - np.sin(slat[q] * np.pi/180))

  for x in range(nl0):
    for y in range(nl1):
        if var[x,y] != fill_value:
            weighted_avg += area[x,y] * var[x,y]
            sumArea += area[x,y]
        else:
            missing_area += area[x,y]

  weighted_avg /= sumArea

  return weighted_avg


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# // calculate global mean profiles //
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def calc_gmean_profiles(lon, lat, var):
    nlev = var.shape[0]
    var_gmean = np.zeros((nlev), dtype=float)

    for z in range(nlev):
        temp_in      = var[z,:,:] ; temp_in = np.squeeze(temp_in)
        temp_out     =  exo.area_weighted_avg(lon, lat, temp_in)
        var_gmean[z] = temp_out

    return var_gmean
