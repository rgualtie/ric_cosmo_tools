import numpy as np
from glob import glob
import os, argparse
import matplotlib.pylab as plt

parser = argparse.ArgumentParser()
parser.add_argument('--freq',type=str, default=90,
                   help = 'Frequency, default = 90')
parser.add_argument('--store', action = 'store_true',
                   help = 'Store the results on the grid')
parser.add_argument('--plot', action = 'store_true',
                   help = 'Create and store plots')
pargs = parser.parse_args()

files = sorted(glob('/sptgrid/user/rgualtie/EnoB/output_spectra/'+pargs.freq+'GHz/*'))
# file prototype: spectra_EnoBlens_s1084_idx084_90GHz.npz
print('%i files found to estimate the covariance matrix'%len(files))

BB = []
TT = []
EE = []
wonky = []
for i,f in enumerate(files):
    container = np.load(f)
    ellb = container['ellb']
    dlsb = container['dlsb']
    lfac = ellb*(ellb+1)/np.pi/4.
    #res = np.sum(np.abs(th_spec_bin[0][:59]-7.43e3*dlsb[0]/lfac))
    #if res < 60000:
    TT.append(dlsb[0])
    EE.append(dlsb[1])
    BB.append(dlsb[2]) # Only BB is needed 
    #else:
    #    wonky.append(f)
    #    continue
print('%i/%i files used'%(len(files)-len(wonky),len(files)))

TTtf = np.mean(TT/th_spec_bin[0][:len(ellb)], axis=0)
np.save(pargs.freq+'GHz_TTtf.npy', TTtf)
os.system('gfal-copy -f '+pargs.freq+'GHz_TTtf.npy gsiftp://osg-gridftp.grid.uchicago.edu/sptgrid/user/rgualtie/EnoB/')
os.system('rm '+pargs.freq+'GHz_TTtf.npy')

EEtf = np.mean(EE/th_spec_bin[1][:len(ellb)], axis=0)
np.save(pargs.freq+'GHz_EEtf.npy', EEtf)
os.system('gfal-copy -f '+pargs.freq+'GHz_EEtf.npy gsiftp://osg-gridftp.grid.uchicago.edu/sptgrid/user/rgualtie/EnoB/')
os.system('rm '+pargs.freq+'GHz_EEtf.npy')

BBtf = np.mean(BB/th_spec_bin[2][:len(ellb)], axis=0)
np.save(pargs.freq+'GHz_BBtf.npy', BBtf)
os.system('gfal-copy -f '+pargs.freq+'GHz_BBtf.npy gsiftp://osg-gridftp.grid.uchicago.edu/sptgrid/user/rgualtie/EnoB/')
os.system('rm '+pargs.freq+'GHz_BBtf.npy')
 
BB = np.array(BB)#/np.mean(BB, axis=0)
cov_mat = np.cov(BB.T)

store = True 
if pargs.store:
    store = True

if pargs.plot:
    plot = True
if plot:
    class MidpointNormalize(Normalize):
        def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
            self.midpoint = midpoint
            Normalize.__init__(self, vmin, vmax, clip)
        def __call__(self, value, clip=None):
            # I'm ignoring masked values and all kinds of edge cases to make a
            # simple example...
            x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
            return np.ma.masked_array(np.interp(value, x, y))
    fig, ax = plt.subplots()
    plt.suptitle('SPT-3G '+pargs.freq+'GHz $C_\ell^{BB}$ covariance matrix')
    norm = MidpointNormalize(midpoint=0)
    flip_cov = cov_mat
    im = ax.imshow(flip_cov[:20, :20], norm=norm, cmap=plt.cm.ocean, origin='lower', interpolation='nearest', extent=[ellb.min(), ellb[20], ellb.min(), ellb[20]])
    fig.colorbar(im, label='[$\mu K^2\cdot\mu K^2]$')
    plt.xlabel(r"Multipole $\ell$")
    plt.ylabel(r"Multipole $\ell'$")
    plt.figure()
    for i in range(BB.shape[0]):
        plt.plot(ellb, BB[i], 'k-', alpha=.1)
    plt.plot(np.arange(len(th_specs[2])), th_specs[2], 'g-', label='$\Lambda -CDM$')
    plt.errorbar(ellb, np.mean(BB, axis=0), yerr=np.std(BB, axis=0), fmt='o', label='mean$\pm \sigma$')
    plt.xlabel('Multipole $\ell$')
    plt.ylabel('$D_\ell\quad[\mu K^2]$')
    plt.xlim(0,800)
    plt.ylim(-.1,None)
    plt.suptitle('EnoB mock-observed ensamble')
    plt.grid();plt.legend()
    plt.figure()
    plt.plot(ellb, flip_cov.diagonal(0), label='Main diagonal')
    plt.plot(ellb, np.fliplr(flip_cov).diagonal(0), label='Secondary diagonal')
    plt.grid();plt.legend()
    if store:
        plt.savefig('lowEllBB_cov_mat.png')
        os.system('gfal-copy -f lowEllBB_cov_mat.png gsiftp://osg-gridftp.grid.uchicago.edu/sptgrid/user/rgualtie/EnoB/')
        os.system('rm lowEllBB_cov_mat.png')
        np.savetxt('lowEllBB_cov_mat.dat', cov_mat.T, delimiter=' ', newline='\n')
        os.system('gfal-copy -f lowEllBB_cov_mat.dat gsiftp://osg-gridftp.grid.uchicago.edu/sptgrid/user/rgualtie/EnoB/')
        os.system('rm lowEllBB_cov_mat.dat')
        plt.savefig('lowEllBB_ensamble.png')
        os.system('gfal-copy -f lowEllBB_ensamble.png gsiftp://osg-gridftp.grid.uchicago.edu/sptgrid/user/rgualtie/EnoB/')
        os.system('rm lowEllBB_ensamble.png')
    
