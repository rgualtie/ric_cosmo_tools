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
        print('Loading map for %s at %s GHz'%(field,freq))
        mm[freq][field] = rt.rotate_map(inmap, coord=['G','C'])
        #for i in [1,2]:
        #    mm[freq][field][i] = 0.
        #hp.write_map('/big_scratch/rgualtie/P'+freq+'Tonly_'+field+'.fits', mm[freq][field], overwrite=True)
        #os.system('gfal-copy -p /big_scratch/rgualtie/P'+freq+'Tonly_'+field+'.fits gsiftp://osg-gridftp.grid.uchicago.edu/sptgrid/user/rgualtie/PlanckTonlyMaps/P'+freq+'Tonly_'+field+'.fits')
        #os.system('rm /big_scratch/rgualtie/P'+freq+'Tonly_'+field+'.fits')

for field in fields:
    for freq in ['100','143','217']:
        print('Building template for %s'%freq)
        template = mm['353'][field]-mm[freq][field]
        hp.write_map('/big_scratch/rgualtie/P353-'+freq+'_template_'+field+'.fits', template, overwrite=True)
        os.system('gfal-copy -p /big_scratch/rgualtie/P353-'+freq+'_template_'+field+'.fits gsiftp://osg-gridftp.grid.uchicago.edu/sptgrid/user/rgualtie/PlanckDR3Maps')
        os.system('rm /big_scratch/rgualtie/P353-'+freq+'_template_'+field+'.fits')
