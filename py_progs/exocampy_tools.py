import numpy as np


def hybrid2pressure(nlon,nlat,nlev,PS,P0,hyam,hybm,hyai,hybi):
  """
  AUTHOR WOLF E.T.
  7/11/2008
  Translated to Python by R. Deitrick, October 2022
  -------------------------------------------------------------
  PURPOSE:  Convert WACCM hybrid level coordinates to pressure coordinates

  NOTES:  for implementation can change lon, lat, lev passing
  parameters to nlon, nlat, nlev to improve speed.  This is meant to be
  used as a function within larger data analysis packages

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
  ilev_P = hyam[:,None,None]*P0 + hybm[:,None,None]*PS[None,:,:]

  return (lev_P, ilev_P)
