#This code is based on frq_scaling.m written by Colin bischoff and translated from Matlab to python. 
#Expanded by Caterina Umilta 
#More details on bandpasses can be found at https://cmb-s4.org/wiki/index.php/Bandpass_Convention_-_What_does_flat_mean

import numpy as np

h = 6.62606957e-34 # J*s
kB = 1.3806488e-23 # J/K
Tcmb = 2.72548 # K

def sync_scaling(bandpass, beta):

    #power law
    scale_fac = bandpass[:,0]**(2+beta)

    if bandpass.shape[0] > 1 :
        #Integrate over bandpass
        N = bandpass.shape[0]
        dnu=np.zeros((N))
        dnu[0] = bandpass[1,0] - bandpass[0,0];
        dnu[1:-1] = (bandpass[2:,0] - bandpass[1:-1,0]) / 2.0
        dnu[-1] = bandpass[-1,0] - bandpass[-2,0]

        scale_fac = np.sum(dnu * scale_fac * bandpass[:,1]) / np.sum(dnu * bandpass[:,1])

    return scale_fac


def dust_scaling(bandpass, beta, Td):

    #Grey-body
    scale_fac= bandpass[:,0]**(3+beta) / (np.exp(h*bandpass[:,0]*10**9 / (kB*Td)) - 1.)

    if bandpass.shape[0] > 1 :
        #Integrate over bandpass
        N = bandpass.shape[0]
        dnu=np.zeros((N))
        dnu[0] = bandpass[1,0] - bandpass[0,0];
        dnu[1:-1] = (bandpass[2:,0] - bandpass[1:-1,0]) / 2.0
        dnu[-1] = bandpass[-1,0] - bandpass[-2,0]

        scale_fac = np.sum(dnu * scale_fac * bandpass[:,1]) / np.sum(dnu * bandpass[:,1])

    return scale_fac


def convfac(bandpass):

    # Conversion factor for thermodynamic temperature.
    conv_fac = bandpass[:,0]**4 * np.exp(h*bandpass[:,0]*10**9 / (kB*Tcmb)) /(np.exp(h*bandpass[:,0]*10**9 / (kB*Tcmb)) - 1.)**2
    if bandpass.shape[0] > 1:
        #Integrate over bandpass
        N = bandpass.shape[0]
        dnu=np.zeros((N))
        dnu[0] = bandpass[1,0] - bandpass[0,0];
        dnu[1:-1] = (bandpass[2:,0] - bandpass[1:-1,0]) / 2.0
        dnu[-1] = bandpass[-1,0] - bandpass[-2,0]
        conv_fac = np.sum(dnu * conv_fac * bandpass[:,1]) / np.sum(dnu * bandpass[:,1])
    return conv_fac


def freq_scaling(bandpass, nu0, beta, greybody=True, Td=19.6):


    # Divide by scale factor at nu0.
    if greybody:
        scale_fac = dust_scaling(bandpass, beta, Td) / dust_scaling(nu0, beta, Td)

    else:
        scale_fac = sync_scaling(bandpass, beta) / sync_scaling(nu0, beta)
    
    conv_fac=convfac(bandpass)

    # Divide by thermodynamic temperature conversion at nu0
    conv_fac = convfac(bandpass) / convfac(nu0)

    # Combine scaling model and thermodynamic conversion
    scale_fac = scale_fac / conv_fac

    return scale_fac


def make_tophat(p, precise=False):

    if precise:
        nu=np.linspace(p[2],p[3],num=100)
    else:
        nu=np.linspace(p[0]-p[1]/2.,p[0]+p[1]/2.,num=100)
    bandpass=np.array([nu,nu**-2]).T

    return bandpass

