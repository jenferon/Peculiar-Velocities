import numpy as np
import astropy.io.fits as fits

nu = np.arange(140.0,150.1,0.1)
base = '/gpfs01/home/ppxjf3/peculiar_vel/data/NOISE/'
LC = np.empty([101,512,512])

for ii in range(0,len(nu)):
	with fits.open(base + 'Noise_N512_FOV1.000_'+str(round(nu[ii],1))+'MHz_SKA_SKA_all_stations_Rev3_EOR0_1.0_0512_natural_I_I.fits', memmap=True) as hdu:
		LC[ii,:,:] = np.array(hdu[0].data) #10MHz centered on 166MHz

fits.writeto(base + 'NOISE_I_140-180_SKA_all_stations_Rev3_EOR0_1.0degree'+str(np.random.randint(1000000000,size=1)[0])+'_uvmax=1000m.fits', LC, overwrite=True)
"""
LC = np.empty([201,512,512])

for ii in range(0,len(nu)):
        with fits.open(base + 'Noise_N512_FOV1.000_'+str(round(nu[ii],1))+'MHz_SKA_SKA_central_area_Rev3_EOR0_1.0_0512_natural_PSF_PSF.fits', memmap=True) as hdu:
                LC[ii,:,:] = np.array(hdu[0].data) #10MHz centered on 166MHz

fits.writeto(base + 'NOISE_PSF_140-160_SKA_central_area_Rev3_EOR0_1.0.fits', LC, overwrite=True)

"""
