import numpy as np
import sys
import math
from astropy.io import fits
import logging


logger = logging.getLogger(__name__)
logging.basicConfig(filename='example.log', encoding='utf-8', level=logging.DEBUG)

def noise_normalisation(base, fov, dim, t_int_h, nfreq, del_nu, nu1):
    """_summary_

    Args:
        base (_type_): file location
        fov (_type_): field of view of the observation
        dim (_type_): pixel dimensiom of the image
        t_int_h (_type_): intergration time of the observation in hours
        nfreq (_type_): number of frequencies in the observation    
        del_nu (_type_): frequency spacing of the observation
        nu1 (_type_): initial frequency of the observation
    """

    #################### Advanced Use #################### 
    #Below is defined as for the SKAO-Low interferometer
    n_tiles = 1 # per station
    n_dipoles = 256 # per tile
    tel = 'SKA_central_area' #change between SKA_central_area and SKA_core_area to define now many stations are being used
    Di = 8 # from baseline design.
    const_PSF = 0
    if (tel=='SKA_core_area'):
        n_stations = 224
        D =  1000  
    elif (tel=='SKA_central_area'):
        n_stations = 296
        D =  3400

    c = 2.998e8 #speed of light
    eta_a = 1.0 #antenna efficiency
    pi = 3.14159
    k_b = 1.38e-23 #boltzman constant
    eta_s=0.9 #system efficiency

    t_int = t_int_h*3600. # integration time / s

    FWHM_arcmin = 4.0
    sz=[dim,dim]

    W = 1.3 #weighting factor (1.3-2) 

    slice_norm = np.zeros(shape=(sz[0],sz[1],nfreq))
    
    #open data cube 
    hdulist = fits.open(base+'NOISE_I_75-95_SKA_central_area_Rev3_EOR0_1.0degree.fits')
    cube = hdulist[0].data#
    cube = np.squeeze(cube)
    for kk in range(0,nfreq):
        nu = nu1 + (del_nu*1e-6)*kk
        lamb = c/(nu*1e6)
        
        # Calculate T_sys
        T_sys = 40.0 + 1.1*60.0*(nu*1e6/c)**(-2.55) # from Table 3 https://astronomers.skatelescope.org/wp-content/uploads/2016/12/SKA-TEL-SKO-DD-001-1_BaselineDesign1.pdf
        # Calculate effective area at each frequency
        eta_rad = (0.056 * nu +82.2)/100.0 # worked out since 85 at 50 MHz and 99 at 300. asssumed linear.
        A_di=lamb**2/(4*pi)*eta_rad*Di #from Table 3 https://astronomers.skatelescope.org/wp-content/uploads/2016/12/SKA-TEL-SKO-DD-001-1_BaselineDesign1.pdf foot note
        if (A_di > 3.2):
            A_di=3.2 # Lower frequencies limited by antenna area
        A=A_di*n_dipoles*n_tiles # Area per station.
        # Calculate SEFD at each frequency.
        K=(A)/(2.0*k_b)
        SEFD=T_sys/K # in J m^-2=W m^-2 Hz^-1
        SEFD=SEFD*1.0e26 # in Jy
        # Calculate noise sensitivities in Jy 
        noise_Jy=(W/eta_s)*(SEFD/math.sqrt(2*n_stations*(n_stations-1)*del_nu*t_int))

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

        slice_norm[:,:,kk] = cube[:,:,kk] * (noise_Jy/rms) # the slice is now in Janskys
        slice_norm [:,:,kk]= slice_norm[:,:,kk] *1e-26*(c/(nu*1.0e6))**2*(1.0/(2.0*k_b*beamarea)) # converting from Jy to K.
        logger.debug('kk; {}, frequency: {}, noise [Jy] {}, conversio factor: {}'.format(kk,nu,noise_Jy,10e-3*1e-26*(c/(nu*1.0e6))**2*(1.0/(2.0*k_b*beamarea))))

    return slice_norm