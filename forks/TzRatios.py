import numpy as np, sys, os, scipy as sc
import flatsky_routines as flatsky, misc_tools as tools, foregrounds as fg
data_folder = '/home/mck74/code/cmb_T_of_z/python/data/'
sys_path_folder = '/home/mck74/code/cmb_T_of_z/python/'
sys.path.append(sys_path_folder)

sys.path.append("/home/mck74/code/spt3g_software/build/")
from spt3g import util, core
from spt3g.beams import beam_analysis as beam_mod
import numpy as np
import numpy.random as rd
import matplotlib.pyplot as pl
from spt3g.util.genericutils import add_scratch_to_path
add_scratch_to_path("mck74")
import clusterfunctions as clu
import cmb_modules as cmb
import os as os
import importlib
import idlsave 
from spt3g import core, maps, sources
from spt3g.mapspectra import basicmaputils, apodmask, map_analysis
from spt3g.sources import source_fitting
import astropy.io.fits as fits
importlib.reload(clu)

pl.rcParams['figure.figsize'] = [8, 8]
cmap = pl.cm.RdYlBu

#params or supply a params file
dx = .25
boxsize_am = 840. #boxsize in arcmins
ny = nx = int(boxsize_am/dx)
mapparams = [ny, nx, dx, dx]
x1,x2 = -nx/2. * dx, nx/2. * dx
verbose = 0
#lmax = 10000
#el = np.arange(lmax)

#CMB power spectrum
Cls_file = '%s/camb/planck18_TTEEEE_lowl_lowE_lensing_highacc/planck2018_base_plikHM_TTTEEE_lowl_lowE_lensing_lensedCls.dat' %(data_folder)
#file with cluster profile + astro cls
tom_products=idlsave.read('/big_scratch/tcrawfor/full_clusterfinder_save_for_mk_gaussfit_noise_27feb22.sav',verbose=False)
profile=tom_products['profiles'][3]#l space cluster profile, coresize = 1 arcmin
profile_ls=np.arange(80001)

apod=tom_products['apm2use'][0]
lpf=None
isohpf=None
hpf=300
xfer_in=None
psize=dx

freq_arr = [90, 150, 220]
noise_val_dict = {90: 3.0, 150: 2.2, 220: 8.8} #uK-arcmin
beam_val_dict = {90: 1.7, 150: 1.2, 220: 1.} #arcmins
color_dict = {90: 'navy', 150: 'darkgreen', 220:'darkred'}
z_arr=np.linspace(0,1.5,8)#z values to use for grid
y_arr=np.linspace(10,20,8)#compton y values to use for grid (factor of 10**-6)
resrad=(dx*core.G3Units.arcmin)/core.G3Units.radians #psize in radians

#get ra, dec or map-pixel grid
ra = np.linspace(x1,x2, nx) #arcmins
dec = np.linspace(x1,x2, nx) #arcmins
ra_grid, dec_grid = np.meshgrid(ra,dec)

def dl_to_cl(el, cl_or_dl, inverse = 0):
    dl_fac = (el * (el+1)/2./np.pi)
    if inverse:
        return cl_or_dl*dl_fac
    else:
        return cl_or_dl/dl_fac

#read Cls now
el, dl_tt, dl_ee, dl_bb, dl_te  = np.loadtxt(Cls_file, unpack = 1)
dl_all = np.asarray( [dl_tt, dl_ee, dl_bb, dl_te] )
cl_all = dl_to_cl(el, dl_all)
cl_tt, cl_ee, cl_bb, cl_te = cl_all #Cls in uK
cl = cl_tt
cl_dic = {}
cl_dic['TT'] = cl_tt

cls_astro=np.array([tom_products['clcmb']*10**12]) #AstroCls in K
ls_astro=np.array([np.arange(len(tom_products['clcmb']))])

#loglog(el, cl_tt)
print(len(el))

#get beam and noise
bl_dict, nl_dict = {}, {}
for freq in freq_arr:
    beamval, noiseval = beam_val_dict[freq], noise_val_dict[freq]
    bl = tools.get_bl(beamval, el, make_2d = 1, mapparams = mapparams)
    nl = tools.get_nl(noiseval, el)
    nl_dict[freq] = nl
    bl_dict[freq] = bl
print(bl_dict.keys(), nl_dict.keys())

#transfer function info
#adding a l of 300 highpass filter
if xfer_in is not None:
    xfer=xfer_in
#if transfer function not given, and any filter cutoffs exist, use those to create transfer function
elif (np.array([lpf,hpf,isohpf])!=None).any():
    xfer=np.fft.fftshift(
        clu.make_filter_xfer(
            lpf=lpf,hpf=hpf,isohpf=isohpf,nx=nx,ny=ny,reso_arcmin=psize))
else:
    print('No Transfer Function Given! Assuming Unity')
    xfer=1
beam_tf={90:bl_dict[90]*xfer,
         150:bl_dict[150]*xfer,
         220:bl_dict[220]*xfer}

#create CMB and convolve with beam
sim_map_dict_1 = {}
for freq in freq_arr:
    #cmb_map = flatsky.make_gaussian_realisation(mapparams, ls_astro[0], cls_astro[0], bl = beam_tf[freq])
    cmb_map = flatsky.make_gaussian_realisation(mapparams, el, cl, bl = beam_tf[freq])    
    #noise map
    noise_map = flatsky.make_gaussian_realisation(mapparams, el, nl_dict[freq])
    sim_map = cmb_map + noise_map
    sim_map_dict_1[freq] = sim_map
    pl.clf()
    pl.figure(figsize=(10, 3.5))
    pl.subplot(131);pl.imshow(cmb_map, extent = [x1,x2,x1,x2], cmap = cmap); pl.colorbar(); pl.title(r'CMB: %s GHz' %(freq))
    pl.subplot(132);pl.imshow(noise_map, extent = [x1,x2,x1,x2], cmap = cmap); pl.colorbar(); pl.title(r'Noise: %s GHz' %(freq))
    pl.subplot(133);pl.imshow(sim_map, extent = [x1,x2,x1,x2], cmap = cmap); pl.colorbar(); pl.title(r'CMB + Noise: %s GHz' %(freq))
    pl.show()

#alternative method
profile_ft2=clu.gridtomap(clu.ell_grid(dx,ny,nx),profile,profile_ls)
real_profile2=np.real(np.fft.fftshift(np.fft.ifft2(profile_ft2,norm='ortho')))
real_profile2/=np.max(real_profile2)
cutout2=real_profile2[mid-200:mid+200,mid-200:mid+200]

#cluster locations
xlocs = np.linspace(700,nx-700, 8)
ylocs = np.linspace(700,nx-700, 8)
x_grid, y_grid = np.meshgrid(xlocs,ylocs)

#create two grids of clusters
#one with alpha = 0, one with alpha = .1
alpha=.1
cluster_map_dic_1={}#alpha=0
cluster_map_dic_2={}#alpha=.1
for freq in freq_arr:
    cluster_map_dic_1[freq]=np.zeros(np.shape(noise_map))
    cluster_map_dic_2[freq]=np.zeros(np.shape(noise_map))
for xi in range(0,8):
    for yi in range(0,8):
        xl=xlocs[xi]
        yl=ylocs[yi]
        z=z_arr[xi]
        cy=y_arr[yi]
        for freq in freq_arr:
            cluster_map_dic_1[freq][int(yl-200):
                                    int(yl+200),
                                    int(xl-200):
                                    int(xl+200)]+=(cutout2*clu.calc_fsz(freq)*2.73*cy)
                                    #int(xl+200)]+=(cutouts[freq]*clu.calc_fsz(freq)*2.73*cy)
            cluster_map_dic_2[freq][int(yl-200):
                                    int(yl+200),
                                    int(xl-200):
                                    int(xl+200)]+=(cutout2*clu.calc_fsz_alpha_z(freq,z,alpha)*2.73*cy)
                                    #int(xl+200)]+=(cutouts[freq]*clu.calc_fsz_alpha_z(freq,z,alpha)*2.73*cy)


for freq in freq_arr:
    cluster_map_dic_1[freq]=np.fft.ifft2(
        np.fft.fft2(cluster_map_dic_1[freq]) * (beam_tf[freq])).real
    cluster_map_dic_2[freq]=np.fft.ifft2(
        np.fft.fft2(cluster_map_dic_2[freq]) * (beam_tf[freq])).real

pl.rcParams['figure.figsize'] = [12, 12]
pl.imshow(cluster_map_dic_1[220])
pl.colorbar()
pl.title('90 ghz cluster grid')
pl.show()

#average SZ decrements (in uK)
avtz90=np.mean([(clu.calc_fsz(90)*2.73*cy) for cy in y_arr])
avtz150=np.mean([(clu.calc_fsz(150)*2.73*cy) for cy in y_arr])
avtz220=np.mean([(clu.calc_fsz(220)*2.73*cy) for cy in y_arr])
print(avtz90)
print(avtz150)
print(avtz220)

print([(clu.calc_fsz(90)*2.73*cy) for cy in y_arr])
print([(clu.calc_fsz(150)*2.73*cy) for cy in y_arr])
print([(clu.calc_fsz(220)*2.73*cy) for cy in y_arr])

#now add clusters to simulated cmb maps
skymaps1=np.array([(sim_map_dict_1[f]+cluster_map_dic_1[f])*apod for f in freq_arr])#alpha=0
skymaps2=np.array([(sim_map_dict_1[f]+cluster_map_dic_2[f])*apod for f in freq_arr])#alpha=.1

print(np.std(skymaps1[2]))

pl.imshow(skymaps1[1])
pl.colorbar()
pl.title('Simulated 150 GHz Map + clusters')
pl.show()

#now we can proceed with the analysis
#first, deconvolve the beams from the 150 and 220 ghz maps
#reconvolve the 150GHz map
skymaps1[1]= np.fft.ifft2( np.fft.fft2(skymaps1[1]) * (bl_dict[90]/bl_dict[150])).real#/factor150
skymaps1[2]=np.fft.ifft2( np.fft.fft2(skymaps1[2]) * (bl_dict[90]/bl_dict[220])).real#/factor220
skymaps2[1]= np.fft.ifft2( np.fft.fft2(skymaps2[1]) * (bl_dict[90]/bl_dict[150])).real#/factor150
skymaps2[2]=np.fft.ifft2( np.fft.fft2(skymaps2[2]) * (bl_dict[90]/bl_dict[220])).real#/factor220

ell_min = 0
ell_max = np.pi / (resrad)
delta_ell = 2 * np.pi / (np.max([nx,ny]) * resrad)
ell_bins = np.arange(ell_min, ell_max + delta_ell, delta_ell)
ell_plot=(np.array(ell_bins[1:]) + np.array(ell_bins[:-1]))/2

#create astro covariance matrix for filter
ncov1d=clu.create_ncov1d(ls_astro,
                         cls_astro,
                         ell_max,
                         1)

#set up instrument noise component for filter
ncl=nl_dict[90][0]
noise_cls=ncl*np.ones(80001)
noise_els=np.arange(80001)
whitenoise=clu.gridtomap(clu.ell_grid(dx,ny,nx),noise_cls,noise_els)
nmatinst=clu.create_N_d(
    [whitenoise],diag=True)

#combine instrument / astro noise to make full covariance matrix
nmatastro=clu.create_ncovbd(ncov1d,np.array([beam_tf[90]]),dx,1,ny,nx,ell_max)
ncovfull=nmatinst+nmatastro

#set up cluster profile for filter
profile_ft=clu.gridtomap(clu.ell_grid(dx,ny,nx),profile,profile_ls)
real_profile=np.fft.fftshift(np.fft.ifft2(profile_ft,norm='ortho'))/dx
norm_factor=1/np.max(np.real(real_profile))
profile_ft*=norm_factor    
ft_signal=clu.create_multi_band_s(
    [profile_ft *beam_tf[90]])

#create filter
psi,sigma=clu.psi_faster([1],ft_signal,ncovfull,ny,nx,1)

psiav=basicmaputils.av_ps(psi,.25*core.G3Units.arcmin,ell_bins,s=np.shape(skymaps1[0]),real=False)
pl.semilogx(ell_plot,psiav)
pl.grid()
pl.title('l space matched filter')
pl.show()

fmaps_1={}
fmaps_2={}

#apply the matched filter
for i in range(0,len(freq_arr)):
    fmap1=clu.multi_band_filter([skymaps1[i]],psi,dx,1)
    fmap2=clu.multi_band_filter([skymaps2[i]],psi,dx,1)
    fmaps_1[freq_arr[i]]=fmap1
    fmaps_2[freq_arr[i]]=fmap2

num1=150
den1=90
num2=220
den2=90

dtsz1={}
dtsz2={}
for freq in freq_arr:
    dtsz1[freq]=clu.ndh_sigma(fmaps_1[freq],apod,mask=None)
    dtsz2[freq]=clu.ndh_sigma(fmaps_2[freq],apod,mask=None)

def sigma_ratio(y,x,sy,sx):#calculate sigma on f=y/x, where y and x have std sy and sx
    return np.sqrt((y**2/x**4)*sx**2+(1/x**2)*sy**2)
def sigma_mean(sigs):#calculate sigma on the mean of an array with the given stds
    return(np.sqrt(np.sum(sigs**2))/len(sigs))

sigma_15090_1=sigma_ratio(fmaps_1[num1],fmaps_1[den1],dtsz1[num1],dtsz1[den1])#np.sqrt(((fmaps_1[150]**2/fmaps_1[90]**4)*dtsz1[90]**2)+(dtsz1[150]*(fmaps_1[90]**-2)))
sigma_15090_2=sigma_ratio(fmaps_2[num1],fmaps_2[den1],dtsz2[num1],dtsz2[den1])#np.sqrt(((fmaps_2[150]**2/fmaps_2[90]**4)*dtsz2[90]**2)+(dtsz2[150]*(fmaps_2[90]**-2)))
sigma_22090_1=sigma_ratio(fmaps_1[num2],fmaps_1[den2],dtsz1[num2],dtsz1[den2])#np.sqrt(((fmaps_1[220]**2/fmaps_1[90]**4)*dtsz1[90]**2)+(dtsz1[220]*(fmaps_1[90]**-2)))
sigma_22090_2=sigma_ratio(fmaps_2[num2],fmaps_2[den2],dtsz2[num2],dtsz2[den2])#np.sqrt(((fmaps_2[220]**2/fmaps_2[90]**4)*dtsz2[90]**2)+(dtsz2[220]*(fmaps_2[90]**-2)))

sigmas_1={150:sigma_15090_1,220:sigma_22090_1}
sigmas_2={150:sigma_15090_2,220:sigma_22090_2}

pl.rcParams['figure.figsize'] = [8, 6]
#at every cluster location, determine the ratio of 90/150 and 90/220
#average in bins of Z and plot
zbins={}
fig1, (ax1,bx1)=pl.subplots(ncols=2,nrows=1,sharey='row')
fig2, (ax2,bx2)=pl.subplots(ncols=2,nrows=1,sharey='row')

sigs=np.array([])
zbins1={}
zbins2={}
dr1={}
dr2={}
for xi in range(0,8):
    zbins1[z_arr[xi]]={150:np.ones(8),
                       220:np.ones(8)}
    zbins2[z_arr[xi]]={150:np.ones(8),
                       220:np.ones(8)}
    dr1[z_arr[xi]]={150:np.ones(8),
                    220:np.ones(8)}
    dr2[z_arr[xi]]={150:np.ones(8),
                    220:np.ones(8)}
    for yi in range(0,8):
        yl=int(ylocs[yi])
        xl=int(xlocs[xi])
        for freq in freq_arr[1:3]:
            if den1==den2:
                zbins1[z_arr[xi]][freq][yi]=fmaps_1[freq][yl,xl]/fmaps_1[den1][yl,xl]
                zbins2[z_arr[xi]][freq][yi]=fmaps_2[freq][yl,xl]/fmaps_2[den1][yl,xl]  
                dr1[z_arr[xi]][freq][yi]=sigmas_1[freq][yl,xl]
                dr2[z_arr[xi]][freq][yi]=sigmas_2[freq][yl,xl]
            else:
                zbins1[z_arr[xi]][freq][yi]=fmaps_1[num1][yl,xl]/fmaps_1[freq][yl,xl]
                zbins2[z_arr[xi]][freq][yi]=fmaps_2[num1][yl,xl]/fmaps_2[freq][yl,xl]  
                dr1[z_arr[xi]][freq][yi]=sigmas_1[freq][yl,xl]
                dr2[z_arr[xi]][freq][yi]=sigmas_2[freq][yl,xl]                
    for freq in freq_arr[1:3]:
        med1=np.median(zbins1[z_arr[xi]][freq])
        med2=np.median(zbins2[z_arr[xi]][freq])
    ax1.scatter(z_arr[xi],np.mean(zbins1[z_arr[xi]][150]))
    ax1.errorbar(z_arr[xi],
                 np.mean(zbins1[z_arr[xi]][150]),
                 sigma_mean(dr1[z_arr[xi]][150]))
    ax2.scatter(z_arr[xi],np.mean(zbins1[z_arr[xi]][220]))
    ax2.errorbar(z_arr[xi],
                 np.mean(zbins1[z_arr[xi]][220]),
                 sigma_mean(dr1[z_arr[xi]][220]))
    bx1.scatter(z_arr[xi],np.mean(zbins2[z_arr[xi]][150]))
    bx1.errorbar(z_arr[xi],
                 np.mean(zbins2[z_arr[xi]][150]),
                 sigma_mean(dr2[z_arr[xi]][150]))
    bx2.scatter(z_arr[xi],np.mean(zbins2[z_arr[xi]][220]))
    bx2.errorbar(z_arr[xi],
                 np.mean(zbins2[z_arr[xi]][220]),
                 sigma_mean(dr2[z_arr[xi]][220]))
ax1.set_title('90/150,alpha=0')
ax2.set_title('90/220, alpha=0')
bx1.set_title('90/150, alpha='+str(alpha))
bx2.set_title('90/220, alpha='+str(alpha))
# ax1.set_ylim([1.5,2.5])
# bx1.set_ylim([1.5,2.5])
#ax2.set_ylim([0,3])
ax1.plot([0,1.5],[(clu.calc_fsz(num1)/clu.calc_fsz(den1)),(clu.calc_fsz(num1)/clu.calc_fsz(den1))],label='alpha=0 prediction')
ax2.plot([0,1.5],[(clu.calc_fsz(num2)/clu.calc_fsz(den2)),(clu.calc_fsz(num2)/clu.calc_fsz(den2))],label='alpha=0.0 prediction')
bx1.plot(z_arr,clu.calc_fsz_alpha_z(num1,z_arr,.1)/clu.calc_fsz_alpha_z(den1,z_arr,.1),label='alpha=0.1 prediction')
bx2.plot(z_arr,(clu.calc_fsz_alpha_z(num2,z_arr,.1)/clu.calc_fsz_alpha_z(den2,z_arr,.1)),label='alpha=0.1 prediction')
ax1.grid()
ax1.legend(loc='best')
ax2.legend(loc='best')
bx1.legend(loc='best')
bx2.legend(loc='best')
ax2.grid()
bx1.grid()
bx2.grid()
fig1.suptitle('Average SZ Ratio vs Redshift')
pl.show()

for key in dr1.keys():
    print(dr1[key][220])

print(1-(.55**2))

arnaud=np.load(data_folder+'arnaud/arnaud_cluster_profiles.npy')
arnaud_dic = np.load(data_folder+'arnaud/arnaud_cluster_profiles.npy', allow_pickle=True, encoding='latin1').item()

#just looking at some values of compton y to see what's reasonable
for key in arnaud_dic.keys():
    print(np.max(arnaud_dic[key])*10**6)


