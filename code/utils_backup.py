# -*- coding: utf-8 -*-
"""
Created on Wed Jul 1 12:09:43 2020

@author: shuwentan-po

Codes that I find useful in futher codings.

To keep things simple this should only import modules from the python standard
library or numpy and scipy.

Inspired by Jesse Cusack

"""

import numpy as np
import scipy.signal as sig
import scipy.io as io
import scipy.stats as stats
from datetime import datetime, timedelta


class Bunch(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        


class data_bin:
    def get_blocks(var, z, z_span, z_center):
        """
        return var & z within z_center +- z_span, and indexes
        both var and z are 1d
        """
        var_blocks = [] 
        z_blocks = [] 
        var_blocks_idx = [] 
        for i in range(len(z_center)):
            mask_block = (z>=z_center[i]-z_span) & (z<=z_center[i]+z_span)
            var_blocks.append(var[mask_block])
            z_blocks.append(z[mask_block])
            var_blocks_idx.append(np.where(mask_block==1)[0])

        return var_blocks, z_blocks, var_blocks_idx



    def mean_blocks(var_blocks):
        """
        compute mean value within blocks
        """
        var_mean = np.zeros(len(var_blocks),) + np.nan
        for i in range(len(var_blocks)):
            var_mean[i] = var_blocks[i].mean()    

        return var_mean



class time_related:
    
    
    
    def sec_in_nyr(n):
        """ return how many seconds in n years (365 days in a year convension) """
        return 3600*24*365*n

    
    def matlab_to_datetime64(t):
        """https://stackoverflow.com/questions/13965740/converting-matlabs-datenum-format-to-python"""
        origin = np.datetime64('0000-01-01', 'D') - np.timedelta64(1, 'D')
        return t * np.timedelta64(1, 'D') + origin



def find_nearest(array, value):
    """ Find the nearest point of
        code stolen from https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx



def sort_profile(var, descend = True):
    """descend = True/False : whether the profile reduces """
    var_sort = var.copy()
    if descend: # np.nanmean(var[0:10])> np.nanmean(var[-10::]):
        idx_sort = np.argsort(var_sort[~np.isnan(var_sort)])[::-1]
        mask_sort = ~np.isnan(var_sort)
        var_sort[~np.isnan(var_sort)] = np.sort(var_sort[~np.isnan(var_sort)])[::-1]
    else:
        idx_sort = np.argsort(var_sort[~np.isnan(var_sort)])
        mask_sort = ~np.isnan(var_sort)
        var_sort[~np.isnan(var_sort)] = np.sort(var_sort[~np.isnan(var_sort)])

    return var_sort, idx_sort, mask_sort



    
    

















def track_dis_lon(trk_lon, trk_lat):
    import gsw
    trk_dis = np.zeros(len(trk_lon), )
    for i in np.arange(1, len(trk_lon)):
        trk_dis[i] = gsw.distance([trk_lon[i], trk_lon[i-1]], [trk_lat[i-1], trk_lat[i-1]], p=0, axis=0)

    return trk_dis



def track_dis_lat(trk_lon, trk_lat):
    import gsw
    trk_dis = np.zeros(len(trk_lon), )
    for i in np.arange(1, len(trk_lon)):
        trk_dis[i] = gsw.distance([trk_lon[i-1], trk_lon[i-1]], [trk_lat[i], trk_lat[i-1]], p=0, axis=0)

    return trk_dis




def lon_per_nm(lon_ref, lat_ref):
    return 1/(track_dis([lon_ref-.5, lon_ref+.5], [lat_ref, lat_ref])[1]/m_per_nm())



def lon_per_m(lon_ref, lat_ref):
    return 1/(track_dis([lon_ref-.5, lon_ref+.5], [lat_ref, lat_ref])[1])




def depth_trk(ds, lon_trk, lat_trk):
    depth = np.zeros(lon_trk.shape) + np.nan
    for i in range(len(depth)):
        if (~np.isnan(lat_trk[i])) & (~np.isnan(lon_trk[i])):
            depth[i] = ds.z.interp(lat=lat_trk[i], lon=lon_trk[i])    
        
    return depth




