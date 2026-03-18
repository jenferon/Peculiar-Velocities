import numpy as np
import matplotib.pyplot as plt

base = '/gpfs01/home/ppxjf3/peculiar_vel/data/Ionization/'
# name = 'xHII_z17.500_N600_L200.0.dat'
z = np.linspace(6.0, 25.0, num=39)
xHII = np.empty([len(z)])
N = 600

for i in range(0,len(z)):
	fname = base + f'xHII_z{z[i]:.3}_N{N:d}_L200.0.dat'
	#open file
	with open(fname, "rb") as fid:
		# Read the binary data from the file
		# data = fid.read()
		# Unpack the binary data into a tuple of floats
		# xHIIi = struct.unpack('f' * (len(data) // struct.calcsize('f')), data)
                xHIIi = np.fromfile(fid, count=N*N*N, dtype=np.float32)
                assert(np.all(xHIIi <= 1.0))
                assert(np.all(xHIIi >= 0.0))
	print('z='+str(i))
	print(np.mean(xHIIi))
	xHII[i] = np.mean(xHIIi)

#save as a txt file
np.savetxt(base + 'xHII_evolution.txt', np.array([z, xHII]).T, header='z <x_HII>_V')
