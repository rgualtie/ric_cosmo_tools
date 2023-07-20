import numpy as np
import os, argparse
from glob import glob

parser = argparse.ArgumentParser()
parser.add_argument('--band', type=str, default='90')
parser.add_argument('--TPmap', action = 'store_true')
parser.add_argument('--store', action = 'store_true')
pargs = parser.parse_args()

TTfiles = list(np.sort(glob('mixmat/'+pargs.band+'GHz/spectra_TT_*.npz')))
EEfiles = list(np.sort(glob('mixmat/'+pargs.band+'GHz/spectra_EE_*.npz')))
BBfiles = list(np.sort(glob('mixmat/'+pargs.band+'GHz/spectra_BB_*.npz')))

nbins = 21
TTblock = np.zeros(6*nbins)
for i,f in enumerate(TTfiles):
    if i==0:
        TTblock = np.load(f)['std'].flatten()
    else:
        TTblock = np.vstack((TTblock, np.load(f)['std'].flatten()))

EEblock = np.zeros(6*nbins)
for i,f in enumerate(EEfiles):
    if i==0:
        EEblock = np.load(f)['std'].flatten()
    else:
        EEblock = np.vstack((EEblock, np.load(f)['std'].flatten()))

BBblock = np.zeros(6*nbins)
for i,f in enumerate(BBfiles):
    if i==0:
        BBblock = np.load(f)['std'].flatten()
    else:
        BBblock = np.vstack((BBblock, np.load(f)['std'].flatten()))

matrix = np.hstack((np.hstack((TTblock.T, EEblock.T)), BBblock.T))

if pargs.TPmap:
    matrix[nbins:, :nbins] = 0.
    matrix[:nbins, nbins:] = 0.

# TE, EB, TB
crossmatrix = np.zeros((6*nbins, 3*nbins))

#Those are nbins x nbins objects
TEblock = np.sqrt(np.abs(TTblock.T[0*nbins: 1*nbins, :]) * np.abs(EEblock.T[1*nbins:2*nbins, :]))
EBblock = np.sqrt(np.abs(EEblock.T[1*nbins: 2*nbins, :]) * np.abs(BBblock.T[2*nbins:3*nbins, :]))
TBblock = np.sqrt(np.abs(TTblock.T[0*nbins: 1*nbins, :]) * np.abs(BBblock.T[2*nbins:3*nbins, :]))

crossmatrix[3*nbins:4*nbins, 0*nbins:1*nbins] = TEblock
crossmatrix[4*nbins:5*nbins, 1*nbins:2*nbins] = EBblock
crossmatrix[5*nbins:6*nbins, 2*nbins:3*nbins] = TBblock

ubermat = np.hstack((matrix, crossmatrix))


if pargs.store:
    if pargs.TPmap:
        np.save('bincenter_mixmatrix_'+pargs.band+'GHz_TPmap_error', ubermat)
    else:
        np.save('bincenter_mixmatrix_'+pargs.band+'GHz_error', ubermat)


