import numpy as np
import pandas as pd
from glob import glob

import os, argparse
parser = argparse.ArgumentParser()
parser.add_argument('--band', type=str, default='150')
parser.add_argument('--store', action = 'store_true')
pargs = parser.parse_args()

# Signal only
files = sorted(glob('/sptgrid/user/rgualtie/EnoB/output_spectra/'+pargs.band+'GHz/spectra_EnoBlens_s1*_idx*_'+pargs.band+'GHz_nstubs256.npz'))
BB = []; TT = []; EE = []
TE = []; EB = []; TB = []
for i,f in enumerate(files):
    container = np.load(f)
    ellb = container['ellb']
    dlsb = container['dlsb']
    lfac = ellb*(ellb+1)/np.pi/4.
    TT.append(dlsb[0])
    EE.append(dlsb[1])
    BB.append(dlsb[2]) # Only BB is needed
    TE.append(dlsb[3])
    EB.append(dlsb[4])
    TB.append(dlsb[5])
TT = np.array(TT)/np.mean(TT, axis=0)
EE = np.array(EE)/np.mean(EE, axis=0)
BB = np.array(BB)/np.mean(BB, axis=0)
TE = np.array(TE)/np.mean(TE, axis=0)
EB = np.array(EB)/np.mean(EB, axis=0)
TB = np.array(TB)/np.mean(TB, axis=0)

dfTT = pd.DataFrame(TT)
dfEE = pd.DataFrame(EE)
dfBB = pd.DataFrame(BB)
dfTE = pd.DataFrame(TE)
dfEB = pd.DataFrame(EB)
dfTB = pd.DataFrame(TB)

cov_mat = np.array([dfTT.cov(), dfEE.cov(), dfBB.cov(), dfTE.cov(), dfEB.cov(), dfTB.cov()])

# Signflip noise covariance
nfiles = np.sort(glob('/sptgrid/user/rgualtie/spectra/'+pargs.band+'GHz_spectra_0*.npz'))
nBB = [];nTT = [];nEE = []
nTE = []; nEB = []; nTB = []
for i,f in enumerate(nfiles):
    container = np.load(f)
    ellb = container['ellb']
    dlsb = 1e6*container['dlsb']
    lfac = ellb*(ellb+1)/np.pi/4.
    nTT.append(dlsb[0])
    nEE.append(dlsb[1])
    nBB.append(dlsb[2]) # Only BB is needed
    nTE.append(dlsb[3])
    nEB.append(dlsb[4])
    nTB.append(dlsb[5])
nTT = np.array(nTT)#/np.mean(TT, axis=0)
nEE = np.array(nEE)#/np.mean(EE, axis=0)
nBB = np.array(nBB)#/np.mean(BB, axis=0)
nTE = np.array(nTE)#/np.mean(TE, axis=0)
nEB = np.array(nEB)#/np.mean(EB, axis=0)
nTB = np.array(nTB)#/np.mean(TB, axis=0)

dfnTT = pd.DataFrame(nTT)
dfnEE = pd.DataFrame(nEE)
dfnBB = pd.DataFrame(nBB)
dfnTE = pd.DataFrame(nTE)
dfnEB = pd.DataFrame(nEB)
dfnTB = pd.DataFrame(nTB)

ncov_mat = np.array([dfnTT.cov(), dfnEE.cov(), dfnBB.cov(), dfnTE.cov(), dfnEB.cov(), dfnTB.cov()])

# Total covariance assuming that cov(S+N)=cov(S)+cov(N)

tot_cov_mat = cov_mat + ncov_mat
np.save(pargs.band+'GHz_signal_covariance_matrix', cov_mat)
np.save(pargs.band+'GHz_noise_covariance_matrix', ncov_mat)
np.save(pargs.band+'GHz_total_covariance_matrix', tot_cov_mat)

plot = True
binmax = 21
ellb = np.load('/sptgrid/user/rgualtie/spectra/bundle_spectra.npz')['ellb']
if plot:
    titles=['TT','EE','BB']
    fig, ax = plt.subplots(4,3, figsize=(10,10))
    plt.suptitle('SPT-3G '+pargs.band+'GHz covariance matrices and $\sqrt{diagonal}$')
    for i in range(3):
        ax[0, i%3].set_title(titles[i])
        im0 = ax[0, i%3].imshow(np.log10(cov_mat[i][:binmax,:binmax]), origin='lower', extent=[ellb.min(), ellb.max(), ellb.min(),ellb.max()], cmap='Blues')
        ax[0, i%3].text(5, 5, 'Signal', bbox={'facecolor': 'white', 'pad': 10})
        plt.colorbar(im0, ax=ax[0, i%3])
        im1 = ax[1, i%3].imshow(np.log10(ncov_mat[i][:binmax,:binmax]), origin='lower', extent=[ellb.min(), ellb.max(), ellb.min(),ellb.max()], cmap='Blues')
        ax[1, i%3].text(5, 5, 'Noise', bbox={'facecolor': 'white', 'pad': 10})
        plt.colorbar(im1, ax=ax[1, i%3])
        im2 = ax[2, i%3].imshow(np.log10(tot_cov_mat[i][:binmax,:binmax]), origin='lower', extent=[ellb.min(), ellb.max(), ellb.min(),ellb.max()], cmap='Blues')
        ax[2, i%3].text(5, 5, 'Total', bbox={'facecolor': 'white', 'pad': 10})
        plt.colorbar(im2, ax=ax[2, i%3])
        ax[3, i%3].plot(ellb, np.sqrt(np.diag(cov_mat[i])), label='Signal')
        ax[3, i%3].plot(ellb, np.sqrt(np.diag(ncov_mat[i])), label='Noise')
        ax[3, i%3].plot(ellb, np.sqrt(np.diag(tot_cov_mat[i])), label='Total')
        ax[3, i%3].set_xlabel('Multipole $\ell$')
        ax[3, i%3].legend();ax[3, i%3].grid()
        ax[i, 0].set_ylabel('Multipole $\ell$')
    ax[3,0].set_ylabel('Power $\mu K^2$')

