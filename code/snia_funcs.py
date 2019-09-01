#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 16:22:05 2019

@author: yuhanyao
"""
import os
import sncosmo
import numpy as np


def add_ZTFfilters():
    cwd = os.getcwd()
    directory = cwd[:-4] + '/data/'
    bandsP48 = {'p48i': 'P48_I.dat',
                'p48r': 'P48_R.dat',
                'p48g': 'P48_g.dat'}

    fileDirectory = directory+'filters/P48/'
    for bandName, fileName in bandsP48.items():
        filePath = os.path.join(fileDirectory, fileName)
        if not os.path.exists(filePath):
            raise IOError("No such file: %s" % filePath)
        b = np.loadtxt(filePath)
        band = sncosmo.Bandpass(b[:, 0], b[:, 1], name=bandName)
        sncosmo.registry.register(band, force=True)
        
        
def mylinear_fit(x, y, yerr, npar = 2):
    '''
    Ref: 
        1. Numerical Recipes, 3rd Edition, p745, 781 - 782
        2. http://web.ipac.caltech.edu/staff/fmasci/ztf/ztf_pipelines_deliverables.pdf, p38
    '''
    assert len(x) == len(y)
    assert len(y) == len(yerr)
    
    Sx = np.sum(x)
    Sy = np.sum(y)
    Sxy = np.sum(x * y)
    Sxx = np.sum(x**2)
    N = len(x)
    
    Sx_sigma = np.sum(x * yerr**2)
    Sxx_sigma = np.sum(x**2 * yerr**2)
    S_sigma = np.sum(yerr**2)
    
    if npar==1:
        Fpsf = Sxy / Sxx
        e_Fpsf = np.sqrt(Sxx_sigma) / Sxx
        a = 0
    elif npar==2:
        Fpsf = (N * Sxy - Sx * Sy) / (N * Sxx - Sx**2)
        a = (Sxx * Sy - Sx * Sxy) / (N * Sxx - Sx**2)
        e_Fpsf = np.sqrt(N**2*Sxx_sigma - 2*N*Sx*Sx_sigma + Sx**2*S_sigma) / (N * Sxx - Sx**2)
    # x_mean = np.mean(x)
    # y_mean = np.mean(y)
    # pearson_r = np.sum( (x - x_mean) * (y - y_mean) ) / np.sqrt(np.sum( (x - x_mean)**2 )) / np.sqrt(np.sum( (y - y_mean)**2 ))
    return Fpsf, e_Fpsf, a
