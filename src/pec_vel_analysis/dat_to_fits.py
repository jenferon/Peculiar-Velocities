import numpy as np
import matplotlib.pyplot as plt
import struct
import astropy.io.fits as fits

base = '/gpfs01/home/ppxjf3/peculiar_vel/data/deltaTb_RSD/'
name = 'Lightcone_N512_FOV1.0000_dnu0.10MHz_165.00MHz_155.00MHz_ds0.016905_div00.00_pv1_oneevent1_evo0_lcon0_test'
# Open the file in binary read mode
with open(base + name +'.dat', "rb") as fid:
    # Read the binary data from the file
    data = fid.read()

    # Unpack the binary data into a tuple of floats
    LC = struct.unpack('f' * (len(data) // struct.calcsize('f')), data)
    
npix = 512 
nfreq = 101

# Convert the tuple to a NumPy array
LC = np.array(LC)
LC = LC.reshape(npix,npix,nfreq)
fits.writeto(base + name + '.fits',  np.transpose(LC, (2,0,1)), overwrite=True)

