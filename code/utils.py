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



    
    
    
def get_nans_blocks_length(self, method = True):
    """
    https://stackoverflow.com/questions/15062205/finding-start-and-stops-of-consecutive-values-block-in-python-numpy-pandas
    method = True: Returns 1D length of np.nan s block in sequence depth wise (last axis).
    method = False: Returns 1D length of ~np.nan s block in sequence depth wise (last axis).
    """
    if method:
        nan_mask = np.isnan(self)
    else:
        nan_mask = ~np.isnan(self)
    start_nans_mask = np.concatenate((np.resize(nan_mask[...,0],self.shape[:-1]+(1,)),
                                    np.logical_and(np.logical_not(nan_mask[...,:-1]), nan_mask[...,1:])
                                    ), axis=self.ndim-1)
    stop_nans_mask = np.concatenate((np.logical_and(nan_mask[...,:-1], np.logical_not(nan_mask[...,1:])),
                                np.resize(nan_mask[...,-1], self.shape[:-1]+(1,))
                                ), axis=self.ndim-1)

    start_idxs = np.where(start_nans_mask)
    stop_idxs = np.where(stop_nans_mask)
    return stop_idxs[-1] - start_idxs[-1] + 1



def butter(cutoff, fs, btype="low", order=4):
    """Return Butterworth filter coefficients. See scipy.signal.butter for a
    more thorough documentation.

    Parameters
    ----------
    cutoff : array
        Cutoff frequency, e.g. roughly speaking, the frequency at which the
        filter acts. Units should be same as for fs paramter.
    fs : float
        Sampling frequency of signal. Units should be same as for cutoff
        parameter.
    btype : {‘lowpass’, ‘highpass’, ‘bandpass’, ‘bandstop’}, optional
        Default is 'low'.
    order : optional, int
        Default is 4. The order of the Butterworth filter.

    Returns
    -------
    b : numpy array
        Filter b coefficients.
    a : numpy array
        Filter a coefficients.

    """
    import scipy.signal as sig
    cutoff = np.asarray(cutoff)
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = sig.butter(order, normal_cutoff, btype=btype, analog=False)
    return b, a   



def butter_filter(x, cutoff, fs, btype="low", order=4, **kwargs):
    """Apply Butterworth filter to data using scipy.signal.filtfilt.

    Parameters
    ----------
    x : array
        The data to be filtered. Should be evenly sampled.
    cutoff : array
        Cutoff frequency, e.g. roughly speaking, the frequency at which the
        filter acts. Units should be same as for fs paramter.
    fs : float
        Sampling frequency of signal. Units should be same as for cutoff
        parameter.
    btype : optional, string
        Default is 'low'. Filter type can be 'low', 'high' or 'band'.
    order : optional, int
        Default is 4. The order of the Butterworth filter.

    Returns
    -------
    y : numpy array
        The filtered data.

    """
    import scipy.signal as sig
    b, a = util.butter(cutoff, fs, btype=btype, order=order)
    y = sig.filtfilt(b, a, x, **kwargs)
    return y


def extract_num_from_str(txt):
    l = ''
    for s in txt:
        if s.isdigit() or s == '.' or s == '-':
            try:
                l += s
            except ValueError:
                pass
    return float(l)



def nan_index_in_var(var):
    """ find leading and trailing indexes for nans in a 1-d array (var) """
    if len(var.shape)>1:
        raise ValueError('var must be a 1-d array')
    else:
        mask = np.isnan(var)
        leading_nans = mask.argmin()
        trailing_nans = mask[::-1].argmin()

        leading_nans_idx = np.arange(leading_nans)
        trailing_nans_idx = np.arange(var.size - trailing_nans, var.size)

        mask[leading_nans_idx] = 0
        mask[trailing_nans_idx] = 0
        middle_nans_idx = np.where(mask==1)[0]

    return leading_nans_idx, trailing_nans_idx, middle_nans_idx






def track_dis(trk_lon, trk_lat):
    import gsw
    trk_dis = np.zeros(len(trk_lon), )
    for i in np.arange(1, len(trk_lon)):
        trk_dis[i] = gsw.distance([trk_lon[i], trk_lon[i-1]], [trk_lat[i], trk_lat[i-1]], p=0, axis=0)

    return trk_dis



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



def nan_index_in_var(var):
    """ find leading and trailing indexes for nans in a 1-d array (var) """
    if len(var.shape)>1:
        raise ValueError('var must be a 1-d array')
    else:
        mask = np.isnan(var)
        leading_nans = mask.argmin()
        trailing_nans = mask[::-1].argmin()

        leading_nans_idx = np.arange(leading_nans)
        trailing_nans_idx = np.arange(var.size - trailing_nans, var.size)

        mask[leading_nans_idx] = 0
        mask[trailing_nans_idx] = 0
        middle_nans_idx = np.where(mask==1)[0]

    return leading_nans_idx, trailing_nans_idx, middle_nans_idx



def m_per_nm():
    return 1851.85



def fine_trk(lon_trk, lat_trk, n):
    lon_trk_new = []
    lat_trk_new = []
    for i in range(len(lon_trk)-1):
        lon_trk_ = np.linspace(lon_trk[i], lon_trk[i+1], n)[0:-1]
        lat_trk_ = np.linspace(lat_trk[i], lat_trk[i+1], n)[0:-1]
        lon_trk_new = np.append(lon_trk_new, lon_trk_)
        lat_trk_new = np.append(lat_trk_new, lat_trk_)
    return lon_trk_new, lat_trk_new



def depth_trk(ds, lon_trk, lat_trk):
    depth = np.zeros(lon_trk.shape) + np.nan
    for i in range(len(depth)):
        if (~np.isnan(lat_trk[i])) & (~np.isnan(lon_trk[i])):
            depth[i] = ds.z.interp(lat=lat_trk[i], lon=lon_trk[i])    
        
    return depth



def lon_per_nm(lon_ref, lat_ref):
    return 1/(track_dis([lon_ref-.5, lon_ref+.5], [lat_ref, lat_ref])[1]/m_per_nm())



def lon_per_m(lon_ref, lat_ref):
    return 1/(track_dis([lon_ref-.5, lon_ref+.5], [lat_ref, lat_ref])[1])



