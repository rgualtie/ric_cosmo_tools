import os, sys
import scipy
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pickle as pkl
from glob import glob
from importlib import reload
sys.path.append('/home/rgualtie/spt3g/RHEL_7_x86_64_py3-v2')
import spt3g
from spt3g import core, mapmaker, util, mapspectra, maps
from spt3g.mapspectra import map_analysis
from copy import copy
import pandas as pd
import healpy as hp
from ispice import *
import camb, pysm

def average_N_spectra(spectra,N_spectra,N_ells):
     avgSpectra = np.zeros(N_ells)
     rmsSpectra = np.zeros(N_ells)
     # calcuate the average spectrum
     i = 0
     while (i < N_spectra):
         avgSpectra = avgSpectra + spectra[i,:]
         i = i + 1
     avgSpectra = avgSpectra/(1. * N_spectra)
     #calculate the rms of the spectrum
     i =0
     while (i < N_spectra):
         rmsSpectra = rmsSpectra +  (spectra[i,:] - avgSpectra)**2
         i = i + 1
     rmsSpectra = np.sqrt(rmsSpectra/(1.*N_spectra))
     #rmsSpectra = np.std(spectra, axis=0)
     return(avgSpectra, rmsSpectra)

def bin_spectrum(cls, lmin=8, lmax=None, binwidth=25, return_error=False):
     cls = np.atleast_2d(cls)
     if lmax is None:
         lmax = cls.shape[-1] - 1
     ell = np.arange(lmax + 1)
     bins = np.arange(lmin, lmax + 1, binwidth)
     ellb = stats.binned_statistic(ell, ell, statistic=np.mean, bins=bins)[0]
     clsb = np.array([stats.binned_statistic(
             ell, C, statistic=np.mean, bins=bins)[0] for C in cls]).squeeze()
     if return_error:
         clse = np.array([stats.binned_statistic(
             ell, C, statistic= np.std, bins=bins)[0] for C in cls]).squeeze()
         return ellb, clsb, clse
     return ellb, clsb

Nsims = 128
freqs = ['90','150','220']
path = '/home/rgualtie/mock_signal_TnoP/total/'
lmin = 8
lmax = 500
ell = np.arange(lmax+1)
lfac = ell*(ell+1)/2./np.pi
bin_width = 25
binEEspecs = np.zeros([len(freqs), Nsims, ell.max()//bin_width-1])

for i,fr in enumerate(freqs):
    for j in range(Nsims):
        m = hp.read_map(path+'total_'+fr+'ghz_map_3g_'+str(j).zfill(4)+'.fits', field=None)
	cls = hp.anafast(map1=m, map2=m, lmax=lmax, pol=True)
	ellb, clsb = bin_spectrum(cls, lmin=lmin, lmax=lmax, binwidth=bin_width)
	binEEspecs[i,j,:] = clsb[1]

