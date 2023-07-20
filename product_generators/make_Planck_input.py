import numpy as np
import healpy as hp
import os
import ric_tools as rt

path = '/sptgrid/user/rgualtie/PlanckDR3Maps/'
freqs = ['100','143','217','353']
fields= ['full','halfmission-1','halfmission-2']
mm = {}
for freq in freqs:
    mm[freq] = {}
    for field in fields:
        inmap = hp.read_map(path+'HFI_SkyMap_'+freq+'_2048_R3.01_'+field+'.fits', field=[0,1,2])
        print('Making map for %s at %s GHz'%(field,freq))
        mm[freq][field] = rt.rotate_map(inmap, coord=['G','C'])
        for i in range(3):
            mm[freq][field][i][np.abs(mm[freq][field][i])>1e8] = 0.
        hp.write_map('/big_scratch/rgualtie/P'+freq+'_'+field+'.fits', mm[freq][field], overwrite=True)
        os.system('gfal-copy /big_scratch/rgualtie/P'+freq+'_'+field+'.fits gsiftp://osg-gridftp.grid.uchicago.edu/sptgrid/user/rgualtie/PlanckMockInput/P'+freq+'_'+field+'.fits')
        os.system('rm /big_scratch/rgualtie/P'+freq+'_'+field+'.fits')
