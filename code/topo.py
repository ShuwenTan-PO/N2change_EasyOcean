# -*- coding: utf-8 -*-
"""
Created on Wed Jul 1 12:09:43 2020

@author: shuwentan-po

Convert topography data from ETOPO1 (.tiff) to netcdf

"""

import numpy as np
import xarray as xr



def read_topo(datafilepath, filename, method=0):
    """read geotiff file using GDAL, and build DataArray
    method = 0: negative longitude (visulize Atlantic Ocean)
    method = 1: all positive longitude (visulize Pacific Ocean)"""
    from osgeo import gdal
    im = gdal.Open(datafilepath + filename)
    topo = im.ReadAsArray()

    nrows, ncols = topo.shape
    x0, dx, dxdy, y0, dydx, dy = im.GetGeoTransform()
    x1 = x0 + dx * ncols
    y1 = y0 + dy * nrows

    lon = np.arange(x0, x1, dx)
    lat = np.arange(y0, y1, dy)

    if method == 0:
        # create DataArray (see more: https://xarray.pydata.org/en/stable/user-guide/data-structures.html)
        Topo = xr.DataArray(topo.T, dims=['longitude', 'latitude'],
                                 coords={'longitude': lon,
                                         'latitude': lat},
                                 attrs={'long_name': 'Depth (m)'})
    elif method == 1:
        lon[lon<0] += 360
        lon_idx = np.argsort(lon)
        Topo = xr.DataArray(topo[:,lon_idx].T, dims=['longitude', 'latitude'],
                                 coords={'longitude': lon[lon_idx],
                                         'latitude': lat},
                                 attrs={'long_name': 'Depth (m)'})
    else:
        raise ValueError('method must be integer from 0 or 1')
        
    return Topo