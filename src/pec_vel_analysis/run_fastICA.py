import numpy as np
import astropy.io.fits as fits
from sklearn.decomposition import FastICA

def do_fastica(data, comps):
    shape = data.shape
    f_ica = FastICA(n_components=comps)
    #generate the 4 componets
    S = f_ica.fit_transform(data.reshape((shape[0]*shape[1],shape[2])))
    
    #get mixing matrix
    A = f_ica.mixing_
    
    #make model
    model_fICA = (np.matmul(A,S.T).T).reshape((shape[0],shape[1],shape[2]))
    
    #get resids
    resids_fICA = data - model_fICA #residuals 
    
    return model_fICA, resids_fICA
base = '/gpfs01/home/ppxjf3/peculiar_vel/'

name_LC = 'LightconeRSD_N1024_FOV1.0000_dnu0.10MHz_180.00MHz_140.00MHz_ds0.016905_div00.00_pv1_oneevent0_evo1_lcon1_fresid_on.fits'

with fits.open(base + name_LC, memmap=True) as hdu:
    LC = np.array(hdu[0].data) #10MHz centered on 166MHz
    hdu.info()

LC = np.transpose(LC, (1,2,0))

LC_1 = LC[:,:,0:201]
LC_2 = LC[:,:,100:301]
LC_3 = LC[:,:,200:401]

with fits.open(base + 'Fg_N1024_FOV1.0_140-180MHz_0.1MHz.fits', memmap=True) as hdu:
    Fgs = np.array(hdu[0].data) #10MHz centered on 166MHz
    hdu.info()
Fgs = Fgs*(10e3)
LC_Fg = LC + Fgs

LC_Fg_1 = LC_Fg[:,:,0:201]
LC_Fg_2 = LC_Fg[:,:,100:301]
LC_Fg_3 = LC_Fg[:,:,200:401]

c = 6
model_LC_1, resids_LC_1 = do_fastica(LC_Fg_1, c)
model_LC_2, resids_LC_2 = do_fastica(LC_Fg_2, c)
model_LC_3, resids_LC_3 = do_fastica(LC_Fg_3, c)

fits.writeto(base + 'residuals_x10e3_fg_6c_140_160_LC_1024.fits', resids_LC_1, overwrite=True)
fits.writeto(base + 'residuals_x10e3_fg_6c_150_170_LC_1024.fits', resids_LC_2, overwrite=True)
fits.writeto(base + 'residuals_x10e3_fg_6c_160_180_LC_1024.fits', resids_LC_3, overwrite=True)
"""
fits.writeto(base + 'model_fg_x10e3_4c_140_160_LC_1024.fits', model_LC_1, overwrite=True)
fits.writeto(base + 'model_fg_x10e3_4c_150_170_LC_1024.fits', model_LC_2, overwrite=True)
fits.writeto(base + 'model_fg_x10e3_4c_160_180_LC_1024.fits', model_LC_3, overwrite=True)
"""
