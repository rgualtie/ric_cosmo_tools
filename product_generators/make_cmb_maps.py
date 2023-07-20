import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import constants as cs # the constants module
from cmb_modules import * # the module of functions
N = cs.N
N_iterations = cs.N_iterations
c_min = cs.c_min
c_max = cs.c_max
X_width =cs.X_width
Y_width = cs.Y_width
beam_size_fwhp = cs.beam_size_fwhp
pix_size = cs.pix_size
Number_of_Sources  = cs.Number_of_Sources
Amplitude_of_Sources = cs.Amplitude_of_Sources
Number_of_Sources_EX = cs.Number_of_Sources_EX
Amplitude_of_Sources_EX = cs.Amplitude_of_Sources_EX
Number_of_SZ_Clusters  = cs.Number_of_SZ_Clusters
Mean_Amplitude_of_SZ_Clusters = cs.Mean_Amplitude_of_SZ_Clusters
SZ_beta = cs.SZ_beta
SZ_Theta_core = cs.SZ_Theta_core
white_noise_level = cs.white_noise_level
atmospheric_noise_level = cs.atmospheric_noise_level
one_over_f_noise_level = cs.one_over_f_noise_level

window = (cosine_window(N))

nbig = 10000
# read in the input CMB spectra
ell, DlTT,DlEE,DlBB, DlTE= np.loadtxt("CMB_fiducial_totalCls.dat", usecols=(0,1,2,3,4), unpack=True)
##
ell_big = np.arange(nbig)
DlTT_big = np.zeros(nbig)
DlTT_big[ell.astype(int)] = DlTT
DlEE_big = np.zeros(nbig)
DlEE_big[ell.astype(int)] = DlEE
DlBB_big = np.zeros(nbig)
DlBB_big[ell.astype(int)] = DlBB
DlTE_big = np.zeros(nbig)
DlTE_big[ell.astype(int)] = DlTE
ell = ell_big
DlTT = DlTT_big #+ 1e-3   ### the 1e-3 factor maps plotting easy
DlEE = DlEE_big #+ 1e-3
DlBB = DlBB_big #+ 1e-3
DlTE = DlTE_big

Temp_point_source_spectrum = DlTT[3000]*(ell/3000.)**2.
Pol_point_source_spectrum = DlEE[4500]*(ell/4500.)**2.

#plt.semilogy(Temp_point_source_spectrum, label='T_PS_spec')
#plt.semilogy(Pol_point_source_spectrum, label='E_PS_spec')
#plt.legend()

DlTT_PS = DlTT + Temp_point_source_spectrum   ### these are used for computing the transer functions
DlEE_PS = DlEE + Pol_point_source_spectrum
DlBB_PS = DlBB + Pol_point_source_spectrum
#plt.semilogy(ell,DlTT,'r')
#plt.semilogy(ell,DlEE,'g')
#plt.semilogy(ell,DlBB,'b')
#plt.xlim(0,2500)
#plt.title('TT (red), EE (green), BB (blue) spectra')
#plt.ylabel('$D^{XX}_{\ell}$ [$\mu$K$^2$]')
#plt.xlabel('$\ell$')
#plt.plot(ell,DlTE,'y')
#plt.ylabel('$D^{TE}_{\ell}$ [$\mu$K$^2$]')
#plt.xlabel('$\ell$')
#plt.title('TE spectrum')
#plt.xlim(0,2500)

def make_CMB_maps(N,pix_size,ell,DlTT,DlEE,DlTE,DlBB):
    "makes a realization of a simulated CMB sky map"
    # convert Dl to Cl
    ClTT = DlTT*2*np.pi/(ell*(ell+1.))
    ClEE = DlEE * 2 * np.pi / (ell*(ell+1.))
    ClTE = DlTE * 2 * np.pi / (ell*(ell+1.))
    ClBB = DlBB * 2 * np.pi / (ell*(ell+1.))
    ## set the \ell = 0 and \ell =1 modes to zero as these are unmeasurable and blow up with the above transform
    ClTT[0:2] = 0.
    ClEE[0:2] = 0.
    ClTE[0:2] = 0.
    ClBB[0:2] = 0.
    ## seperate the correlated and uncorrelated part of the EE spectrum
    correlated_part_of_E = ClTE**2./(ClTT+1e-12)
    uncorrelated_part_of_EE = ClEE - ClTE**2./(ClTT+1e-12)
    correlated_part_of_E[0:2] = 0.
    uncorrelated_part_of_EE[0:2] = 0.
    # make a 2d coordinate system
    ones = np.ones(N)
    inds  = (np.arange(N)+.5 - N/2.) /(N-1.)
    X = np.outer(ones,inds)
    Y = np.transpose(X)
    R = np.sqrt(X**2. + Y**2.)
    ang = np.arctan2(Y,X)   ## we now need this angle to handled the EB <--> QU rotation
    # now make a set of 2d CMB masks for the T, E, and B maps
    ell_scale_factor = 2. * np.pi / (pix_size/60. * np.pi/180.)
    ell2d = R * ell_scale_factor
    ClTT_expanded = np.zeros(int(ell2d.max())+1)
    ClTT_expanded[0:(ClTT.size)] = ClTT
    ClEE_uncor_expanded = np.zeros(int(ell2d.max())+1)
    ClEE_uncor_expanded[0:(uncorrelated_part_of_EE.size)] = uncorrelated_part_of_EE
    ClE_corr_expanded = np.zeros(int(ell2d.max())+1)
    ClE_corr_expanded[0:(correlated_part_of_E.size)] = correlated_part_of_E
    ClBB_expanded = np.zeros(int(ell2d.max())+1)
    ClBB_expanded[0:(ClBB.size)] = ClBB
    CLTT2d = ClTT_expanded[ell2d.astype(int)]
    ClEE_uncor_2d = ClEE_uncor_expanded[ell2d.astype(int)]
    ClE_corr2d = ClE_corr_expanded[ell2d.astype(int)]
    CLBB2d = ClBB_expanded[ell2d.astype(int)]
    # now make a set of gaussin random fields that will be turned into the CMB maps
    random_array_for_T = np.fft.fft2(np.random.normal(0,1,(N,N)))
    random_array_for_E = np.fft.fft2(np.random.normal(0,1,(N,N)))
    random_array_for_B = np.fft.fft2(np.random.normal(0,1,(N,N)))
    ## make the T, E, and B maps by multiplyign the masks against the random fields
    FT_2d = np.sqrt(CLTT2d) * random_array_for_T
    FE_2d = np.sqrt(ClEE_uncor_2d) * random_array_for_E + ClE_corr2d* random_array_for_T
    FB_2d = np.sqrt(CLBB2d) * random_array_for_B
    ## now conver E abd B to Q and U
    FQ_2d = FE_2d* np.cos(2.*ang) - FB_2d * np.sin(2. *ang)
    FU_2d = FE_2d* np.sin(2.*ang) + FB_2d * np.cos(2. *ang)
    ## convert from fourier space to real space
    CMB_T = np.fft.ifft2(np.fft.fftshift(FT_2d)) /(pix_size /60.* np.pi/180.)
    CMB_T = np.real(CMB_T)
    CMB_Q = np.fft.ifft2(np.fft.fftshift(FQ_2d)) /(pix_size /60.* np.pi/180.)
    CMB_Q = np.real(CMB_Q)
    CMB_U = np.fft.ifft2(np.fft.fftshift(FU_2d)) /(pix_size /60.* np.pi/180.)
    CMB_U = np.real(CMB_U)
    ## optional code for spitting out E and B maps
    CMB_E = np.fft.ifft2(np.fft.fftshift(FE_2d)) /(pix_size /60.* np.pi/180.)
    CMB_E = np.real(CMB_E)
    CMB_B = np.fft.ifft2(np.fft.fftshift(FB_2d)) /(pix_size /60.* np.pi/180.)
    CMB_B = np.real(CMB_B)
    ## return the maps
    return(CMB_T,CMB_Q,CMB_U,CMB_E,CMB_B)
  ###############################
## make a CMB map
CMB_T,CMB_Q,CMB_U,CMB_E,CMB_B = make_CMB_maps(N,pix_size,ell,DlTT,DlEE,DlTE,DlBB)
p = Plot_CMB_Map(CMB_T,c_min,c_max,X_width,Y_width)
p = Plot_CMB_Map(CMB_Q,c_min/20.,c_max/20.,X_width,Y_width)
p = Plot_CMB_Map(CMB_U,c_min/20.,c_max/20.,X_width,Y_width)
p = Plot_CMB_Map(CMB_E,c_min/20.,c_max/20.,X_width,Y_width)
p = Plot_CMB_Map(CMB_B,c_min/20.,c_max/20.,X_width,Y_width)
#Plot_CMB_Map(CMB_T,c_min,c_max,X_width,Y_width)
PSMap = Poisson_source_component(N,pix_size,Number_of_Sources,Amplitude_of_Sources)
PSMap += Exponential_source_component(N,pix_size,Number_of_Sources_EX,Amplitude_of_Sources_EX)
#Add polarization angle to the sources..??
## make an SZ map
SZMap,SZCat = SZ_source_component(N,pix_size,Number_of_SZ_Clusters,Mean_Amplitude_of_SZ_Clusters,SZ_beta,SZ_Theta_core,False)
nCMB_T = CMB_T+PSMap+SZMap
nCMB_Q = CMB_Q+PSMap*0.03#+SZMap
nCMB_U = CMB_U+PSMap*0.03#+SZMap
nCMB_E = CMB_E+PSMap*0.03
nCMB_B = CMB_B+PSMap*0.03
CMB_T_convolved = convolve_map_with_gaussian_beam(N,pix_size,beam_size_fwhp,nCMB_T)
CMB_Q_convolved = convolve_map_with_gaussian_beam(N,pix_size,beam_size_fwhp,nCMB_Q)
CMB_U_convolved = convolve_map_with_gaussian_beam(N,pix_size,beam_size_fwhp,nCMB_U)
CMB_E_convolved = convolve_map_with_gaussian_beam(N,pix_size,beam_size_fwhp,nCMB_E)
CMB_B_convolved = convolve_map_with_gaussian_beam(N,pix_size,beam_size_fwhp,nCMB_B)
Noise = make_noise_map(N,pix_size,white_noise_level,atmospheric_noise_level,one_over_f_noise_level)
UberTMap = CMB_T_convolved + Noise
UberQMap = CMB_Q_convolved + Noise*np.sqrt(2)
UberUMap = CMB_U_convolved + Noise*np.sqrt(2)
UberEMap = CMB_E_convolved + Noise*np.sqrt(2)
UberBMap = CMB_B_convolved + Noise*np.sqrt(2)
p = Plot_CMB_Map(UberBMap*window,c_min/10.,c_max/10.,X_width,Y_width)
