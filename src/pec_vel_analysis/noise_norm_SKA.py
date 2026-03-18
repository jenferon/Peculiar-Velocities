#!/usr/bin/python
# python noise_norm_SKA.py
# normalise any image rms to expected rms of noise in Kelvin
import numpy as np
import sys
import math
from astropy.io import fits

path = '/gpfs01/home/ppxjf3/peculiar_vel/data/NOISE/' #location of noise images
path_out = '/gpfs01/home/ppxjf3/peculiar_vel/data/NOISE/' #location for output
fov=1.0 # field of view in degrees
Dim = 512 # pixels per side
t_int_h = 1000 # number of hours to normalise noise to
n_f = 101 # number of frequencies
del_nu = 0.1e6 # frequency interval / Hz
f1 = 140.0 # first frequency map / MHz

#################### Advanced Use #################### 
#SKA

n_tiles = 1 # per station
n_dipoles = 256 # per tile
tel = 'SKA_all_stations'
Di = 8 # from baseline design.
const_PSF = 0
if (tel=='SKA_core_area'):
    n_stations = 224
    D =  1000  
elif (tel=='SKA_central_area'):
    n_stations = 296
    D =  3400
elif (tel=='SKA_all_stations'):
    n_stations = 512
    D = 3400

c = 2.998e8
eta_a = 1.0 #antenna efficiency
pi = 3.14159
k_b = 1.38e-23
eta_s=0.9 #system efficiency

t_int = t_int_h*3600. # integration time / s

FWHM_arcmin = 4.0
sz=[Dim,Dim]

W = 1.3 #weighting factor (1.3-2) #WE SET THIS AS 1 FOR SUSSEX but also didnt have factor of 2. with equation as they are now, set to 2 to get 67 mK at 150 MHz for 48 stattions. 1.3 89 mK for 24 stations . set to 1.3 to match website.

slice_norm = np.zeros(shape=(sz[0],sz[1],n_f))
for kk in range(0,n_f):
    freq = f1 + (del_nu*1e-6)*kk
    hdulist = fits.open(path+'Noise_N512_FOV1.000_'+str(np.round(freq,1))+'MHz_SKA_SKA_all_stations_Rev3_EOR0_1.0_0512_natural_I_I.fits')
    cube = hdulist[0].data#
    cube = np.squeeze(cube)
    lamb = c/(freq*1e6)
    print('Noise_N512_FOV1.000_'+str(np.round(freq,1))+'MHz_SKA_SKA_all_stations_Rev3_EOR0_1.0_0512_natural_I_I.fits')
    # Calculate T_sys
    T_sys = 140.0 + 1.1*60.0*(freq*1e6/c)**(-2.55) # from Table 3 https://astronomers.skatelescope.org/wp-content/uploads/2016/12/SKA-TEL-SKO-DD-001-1_BaselineDesign1.pdf
    # Calculate effective area at each frequency
    eta_rad = (0.056 * freq +82.2)/100.0 # worked out since 85 at 50 MHz and 99 at 300. asssumed linear.
    A_di=lamb**2/(4*pi)*eta_rad*Di #from Table 3 https://astronomers.skatelescope.org/wp-content/uploads/2016/12/SKA-TEL-SKO-DD-001-1_BaselineDesign1.pdf foot note
    if (A_di > 3.2):
        A_di=3.2 # Lower frequencies limited by antenna area
    A=A_di*n_dipoles*n_tiles # Area per station.
 # Calculate SEFD at each frequency.
    K=(A)/(2.0*k_b)
    SEFD=T_sys/K # in J m^-2=W m^-2 Hz^-1
    SEFD=SEFD*1.0e26 # in Jy
 # Calculate noise sensitivities in Jy 
    noise_Jy=(W/eta_s)*(SEFD/math.sqrt(2*n_stations*(n_stations-1)*del_nu))*math.sqrt(1000/(6*10))

    rms =0.
    for ii in range(0,sz[0]):
        for jj in range(0,sz[1]):
            rms = rms + cube[ii,jj]**2
    rms = math.sqrt(rms/(sz[0]*sz[1])) 

    if (const_PSF == 1):
        FWHM = (FWHM_arcmin/60.0)*(pi/180.0)
    else:
        FWHM = 1.22 * lamb / D  # radians

    beamarea = pi * FWHM**2 / (4.0*math.log(2)) # beam solid angle of Gaussian beam

    slice_norm[:,:,kk] = cube[:,:] * (noise_Jy/rms) # the slice is now in Janskys
    slice_norm[:,:,kk] = slice_norm[:,:,kk]*1e-26*(c/(freq*1.0e6))**2*(1.0/(2.0*k_b*beamarea)) # converting from Jy to K.
    print(kk,freq,noise_Jy,10e-3*1e-26*(c/(freq*1.0e6))**2*(1.0/(2.0*k_b*beamarea)))

fits.writeto(path_out+'Noise_N512_FOV1.000_140-150_1000hrs_SKA_all_stations_Rev3_EOR0_1.0_0512_natural_I_test.fits',slice_norm,overwrite=True)
