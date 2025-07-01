import numpy as np
import matplotlib.pyplot as plt
import struct
import astropy.io.fits as fits

base = '/home/ppxjf3/RSD_LC/comparison/ds_2/'
name = 'Oneevent_N512_FOV1.0000_dnu0.10MHz_180.00MHz_140.00MHz_ds0.006762_div00.00_pv1_oneevent0_evo1_lcon0_dz_000.50_ds_2'
# Open the file in binary read mode
with open(base + name +'.dat', "rb") as fid:
    # Read the binary data from the file
    data = fid.read()

    # Unpack the binary data into a tuple of floats
    LC = struct.unpack('i' * (len(data) // struct.calcsize('i')), data)
    
npix = 512
nfreq = 401

# Convert the tuple to a NumPy array
LC = np.array(LC)
LC = LC.reshape(npix,npix,nfreq)
#fits.writeto(base + name + '.fits',  np.transpose(LC, (2,0,1)), overwrite=True)
fits.writeto(base + name + '.fits',  LC, overwrite=True)

