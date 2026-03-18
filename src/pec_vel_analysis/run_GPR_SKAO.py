import numpy as np
import astropy.io.fits as fits
import GPy
from gpr4im import fg_tools as fg
from gpr4im import pk_tools as pk
from gpr4im import obs_tools as obs

def GPR(data, frequency):
    
    #choose kernel
    # kernel for the smooth foreground:
    kern_sfg = GPy.kern.RBF(1)
    #mixing kernel
    #kern_mix = GPy.kern.Matern32(1)
    #ex kernel
    kern_ex = GPy.kern.Matern52(1)
    # kernel for the HI cosmological signal:
    kern_21 = GPy.kern.Exponential(1)
    #set lengthscales to ensure the kernals fit to the correcrt part of the signal based on the data in Mertens et al 2020

    kern_sfg.lengthscale.constrain_bounded(10,100)
    kern_21.lengthscale.constrain_bounded(0.1,1.2)
    #kern_mix.lengthscale.constrain_bounded(1,10)
    kern_ex.lengthscale.constrain_bounded(0.2,8)
    kern_fg = kern_sfg + kern_ex # kern_mix +

    gpr_result = fg.GPRclean(data, frequency, kern_fg, kern_21, NprePCA=0, num_restarts=10,
                                              noise_data=None, heteroscedastic=False, zero_noise=True, invert=False)

    model_gpr = gpr_result.fgfit 
    resids_gpr = data - model_gpr

    return model_gpr, resids_gpr

base = '/gpfs01/home/ppxjf3/peculiar_vel/'

name_LC = 'EoR_H21cm_v9.fits'

with fits.open(base + name_LC, memmap=True) as hdu:
    LC = np.array(hdu[0].data) #10MHz centered on 166MHz
    hdu.info()

LC = np.transpose(LC, (1,2,0))
LC = np.flip(LC, axis=2)
LC = LC[:,:,340:741]

LC_1 = LC[:,:,0:201]
LC_2 = LC[:,:,100:301]
LC_3 = LC[:,:,200:401]

with fits.open(base + 'Fg_N1024_FOV1.0_140-180MHz_0.1MHz.fits', memmap=True) as hdu:
    Fgs = np.array(hdu[0].data) #10MHz centered on 166MHz
    hdu.info()

LC_Fg = LC + Fgs[256:768,256:768,:]

LC_Fg_1 = LC_Fg[:,:,0:201]
LC_Fg_2 = LC_Fg[:,:,100:301]
LC_Fg_3 = LC_Fg[:,:,200:401]

freq = np.arange(140, 180.1, 0.1)
model_GPR, resids_GPR = GPR(LC_Fg, freq)
model_GPR_1, resids_GPR_1 = GPR(LC_Fg_1, freq[:201])
model_GPR_2, resids_GPR_2 = GPR(LC_Fg_2, freq[100:301])
model_GPR_3, resids_GPR_3 = GPR(LC_Fg_3, freq[200:401])

fits.writeto(base + 'residuals_fg_GPR_140_180_LC_SKAO.fits', resids_GPR, overwrite=True)
fits.writeto(base + 'residuals_fg_GPR_140_160_LC_SKAO.fits', resids_GPR_1, overwrite=True)
fits.writeto(base + 'residuals_fg_GPR_150_170_LC_SKAO.fits', resids_GPR_2, overwrite=True)
fits.writeto(base + 'residuals_fg_GPR_160_180_LC_SKAO.fits', resids_GPR_3, overwrite=True)



"""
fits.writeto(base + 'model_fg_GPR_140_160_LC_1024.fits', model_GPR_1, overwrite=True)
fits.writeto(base + 'model_fg_GPR_150_170_LC_1024.fits', model_GPR_2, overwrite=True)
fits.writeto(base + 'model_fg_GPR_160_180_LC_1024.fits', model_GPR_3, overwrite=True)
"""
