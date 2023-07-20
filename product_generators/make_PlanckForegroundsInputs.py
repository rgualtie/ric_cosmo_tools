import numpy as np
import os
from spt3g import core, maps
import ric_tools as rt

P353 = [fr for fr in core.G3File('/sptgrid/simulation/planck_maps/2018/preprocessed/field/HFI_SkyMap_353_2048_R3.01_fullmission_winter.g3')][0]
maps.RemoveWeights(P353, zero_nans=True)
P100 = [fr for fr in core.G3File('/sptgrid/simulation/planck_maps/2018/preprocessed/field/HFI_SkyMap_100_2048_R3.01_fullmission_winter.g3')][0]
maps.RemoveWeights(P100, zero_nans=True)
DTemplate = rt.subtract_G3maps(P353, P100)

core.G3Writer('P353-P100_DustTemplate.g3')(DTemplate)
os.system('gfal-copy P353-P100_DustTemplate.g3 $URI/sptgrid/user/rgualtie/PlanckMockForegroundsInput/')
os.system('rm P353-P100_DustTemplate.g3')

