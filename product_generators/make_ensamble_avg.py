# Average down the mix_matrix power spectra
import numpy as np
import os, argparse
import ric_tools as rt
import matplotlib.pylab as plt
from glob import glob

parser = argparse.ArgumentParser()
parser.add_argument('--infile', type=str, default='/sptgrid/user/rgualtie/mix_matrix/output_spectra/150GHz/spectra_TT_0008',
                   help = 'input npz file root, full path required')
parser.add_argument('--spec_type', type=str, default='TT')
parser.add_argument('--bb', type=str, default='008')
parser.add_argument('--band', type=str, default='150')
parser.add_argument('--nobs', type=str, default='256',
                   help = 'Number of ObsIDs to be used, active if --randomobs is used, default=256')
parser.add_argument('--nbins', type=int, default=21)
parser.add_argument('--store', action = 'store_true')
pargs = parser.parse_args()

bb = os.path.basename(pargs.infile)[-4:]

th_specs = np.genfromtxt('/home/rgualtie/spt3g/spt3g_software/simulations/python/data/camb/planck18_TTEEEE_lowl_lowE_lensing_highacc/planck2018_base_plikHM_TTTEEE_lowl_lowE_lensing_lensedCls.dat', unpack=True)
th_specs = np.vstack((th_specs, np.zeros_like(th_specs[0]), np.zeros_like(th_specs[0])))
lmin = 8
lmax = 6143
lbin = 25
th_specs = th_specs[1:]
th_ellb, th_spec_bin = rt.bin_spectrum(th_specs, lmin=lmin, lmax=None, binwidth=lbin, return_error=False)

bins_dict = {}
for i,j in enumerate(range(8,508+25,25)):
    bins_dict[str(j).zfill(3)] = i
spec_type_dict = {'TT':0, 'EE':1, 'BB':2}

avs = np.load('mixmat_input_averages.npy')

dlsb = {}
errb = {}
ellb = {}
skip = []
for i in range(200):
    try:
        temp = np.load(pargs.infile+'_s10'+str(i).zfill(2)+'_idx0'+str(i).zfill(2)+'_'+pargs.band+'GHz_nstubs'+pargs.nobs+'.npz')
    except FileNotFoundError:
        skip.append(i)
        continue
    ellb[i] = temp['ellb']
    dlsb[i] = temp['dlsb']
    errb[i] = temp['errb']

spec = np.zeros((len(ellb.keys())+1,len(dlsb[10]), len(dlsb[10][0]))) 
for i,k in enumerate(ellb.keys()):
    for j in range(6):
        spec[i,j] = spec[i][j]+dlsb[k][j]

std = np.std(spec, axis=0)
std = np.array(std[:,:pargs.nbins])
m = np.mean(spec, axis=0)
m = np.array(m[:,:pargs.nbins])

#bspec = np.load('/home/rgualtie/camb_models/binned_spectra.npz')

ratio = []
for i,s in enumerate(['TT','EE','BB','TE','EB','TB']):
    ratio.append(m[i]/th_spec_bin[spec_type_dict[pargs.spec_type]][bins_dict[pargs.bb]])
    #ratio.append(m[i]/avs[bins_dict[pargs.bb]+pargs.nbins*i])
ratio = np.array(ratio)
#ratio[np.abs(ratio)== np.inf] = 0
ratio = ratio.flatten()
column = m.flatten() 

if pargs.store:
    np.savez(os.path.basename(pargs.infile)+'_'+pargs.band+'GHz.npz', mean=m, std=std)
    np.savez(os.path.basename(pargs.infile)+'_'+pargs.band+'GHz_ratio.npz', ratio=ratio)
    np.savez(os.path.basename(pargs.infile)+'_'+pargs.band+'GHz_column.npz', column=column)

    os.system('gfal-copy '+os.path.basename(pargs.infile)+'_'+pargs.band+'GHz.npz gsiftp://osg-gridftp.grid.uchicago.edu/sptgrid/user/rgualtie/mix_matrix/ratios/'+pargs.band+'GHz/')
    os.system('gfal-copy '+os.path.basename(pargs.infile)+'_'+pargs.band+'GHz_ratio.npz gsiftp://osg-gridftp.grid.uchicago.edu/sptgrid/user/rgualtie/mix_matrix/ratios/'+pargs.band+'GHz/')
    os.system('gfal-copy '+os.path.basename(pargs.infile)+'_'+pargs.band+'GHz_column.npz gsiftp://osg-gridftp.grid.uchicago.edu/sptgrid/user/rgualtie/mix_matrix/ratios/'+pargs.band+'GHz/')
    os.system('rm '+os.path.basename(pargs.infile)+'_'+pargs.band+'GHz.npz')
    os.system('rm '+os.path.basename(pargs.infile)+'_'+pargs.band+'GHz_ratio.npz')
    os.system('rm '+os.path.basename(pargs.infile)+'_'+pargs.band+'GHz_column.npz')
