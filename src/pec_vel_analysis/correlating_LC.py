import numpy as np
import struct
import matplotlib.pyplot as plt
def pearson_correl(x,y):
	return (np.sum((x-np.mean(x))*(y-np.mean(y))))/(np.sqrt(np.sum((x-np.mean(x))**2)*np.sum((y-np.mean(y))**2)))
 
"""base = '/gpfs01/home/ppxjf3/peculiar_vel/data/deltaTb_RSD/'
name_LC = 'LightconeRSD_N512_FOV1.0000_dnu0.50MHz_165.00MHz_155.00MHz_ds0.016905_div00.00_pv1_oneevent1_evo0_test.dat'
 
# Open the file in binary read mode
with open(base + name_LC, "rb") as fid:
	# Read the binary data from the file
	data = fid.read()

	# Unpack the binary data into a tuple of floats
	LC = struct.unpack('f' * (len(data) // struct.calcsize('f')), data)
    
npix = 512 
nfreq =21
# Convert the tuple to a NumPy array
LC = np.array(LC)
LC = LC.reshape(npix,npix,nfreq)
"""
z = np.arange(6.0, 26.0, 1.0)
xHI = np.array([])

base_delta = '/gpfs01/home/ppxjf3/peculiar_vel/data/Ionization/'
for ii in z:
	with open(base_delta + 'xHII_z'+str(np.round(ii,1))+'00_N600_L200.0.dat', "rb") as fid:
		# Read the binary data from the file
		data = fid.read()

	# Unpack the binary data into a tuple of floats
	delta = struct.unpack('f' * (len(data) // struct.calcsize('f')), data)
	    
	delta = np.array(delta)
	delta = delta.reshape(600,600,600)
	
	xHI = np.append(xHI, np.mean(delta))
		
np.savetxt('/gpfs01/home/ppxjf3/peculiar_vel/xHII_paper.txt', [z,xHI])
