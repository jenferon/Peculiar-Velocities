import numpy as np
import matplotlib.pyplot as plt
import struct
import astropy.io.fits as fits

base = '/home/ppxjf3/RSD_LC/comparison/ds_2/'
name = 'Oneevent_N512_FOV1.0000_dnu0.10MHz_180.00MHz_140.00MHz_ds0.006762_div00.00_pv1_oneevent0_evo1_lcon0_dz_000.50_ds_2'
        
npix = 512
nfreq = 401
# Open the file in binary read mode

def dat_to_fits(base, name, npix, nfreq):
    """Function to change data type from a .dat file to a .fits file 

    Args:
        base (string): location where the file is saved
        name (string): name of the file
        npix (int): x and y dimension of the image cube
        nfreq (int): z dimension of the image cube
    """
    with open(base + name +'.dat', "rb") as fid:
        # Read the binary data from the file
        data = fid.read()

        # Unpack the binary data into a tuple of floats
        LC = struct.unpack('i' * (len(data) // struct.calcsize('i')), data)


    # Convert the tuple to a NumPy array
    LC = np.array(LC)
    LC = LC.reshape(npix,npix,nfreq)
    fits.writeto(base + name + '.fits',  LC, overwrite=True)

