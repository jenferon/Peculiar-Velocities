import numpy as np
import struct
import astropy.io.fits as fits
import tools21cm as t2c
from sklearn.decomposition import FastICA

from astropy.cosmology import Planck13 as cosmoP
from astropy.cosmology import FlatLambdaCDM, LambdaCDM

import astropy.units as u
cosmo = FlatLambdaCDM(H0=71 * u.km / u.s / u.Mpc, Om0=0.27)

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

def zFromNu(nu):
    """
    Convert frequency of 21cm line to redshift
    
    Input: nu [MHz]
    """
    nu21 = 1.420405e3  #MHz
    return nu21/nu - 1.0

def NuFromz(z):
    """
    Convert frequency of 21cm line to redshift
    
    Input: nu [MHz]
    """
    nu21 = 1.420405e3  #MHz
    return nu21/(z + 1.0)

def get_lengths(nu_low, nu_hi, nu_mid, theta_FOV):
    z_lo = zFromNu(nu_low)
    z_mid = zFromNu(nu_mid)
    z_hi = zFromNu(nu_hi)
    
    L_para = cosmo.comoving_distance(z_lo) - cosmo.comoving_distance(z_hi)
    L_perp = cosmo.comoving_distance(z_mid) * (np.pi * theta_FOV / 180.0)
    return L_para/u.Mpc, L_perp/u.Mpc

base = '/gpfs01/home/ppxjf3/peculiar_vel/data/'
npix=512
nfreq=201

with open(base + 'LightconeRSD_N512_FOV1.0000_dnu0.10MHz_095.00MHz_075.00MHz_ds0.003099_div00.00_pv1_oneevent0_evo1_lcon0_dz_000.10.dat', "rb") as fid:
    # Read the binary data from the file
    data = fid.read()

    # Unpack the binary data into a tuple of floats
    RSD = struct.unpack('f' * (len(data) // struct.calcsize('f')), data)

RSD = np.array(RSD)*1e3
RSD = RSD.reshape(npix,npix,nfreq)
with open(base + 'LightconeRSD_N512_FOV1.0000_dnu0.10MHz_095.00MHz_075.00MHz_ds0.007747_div00.00_pv1_oneevent0_evo1_lcon1_dz_000.10.dat', "rb") as fid:
# Read the binary data from the file
    data = fid.read()

    # Unpack the binary data into a tuple of floats
    LC = struct.unpack('f' * (len(data) // struct.calcsize('f')), data)

LC = np.array(LC)*1e3
LC = LC.reshape(npix,npix,nfreq)

with fits.open(base + 'Fg_N512_FOV1.0_75_95.0MHz_0.1MHz.fits', memmap=True) as hdu:
    FG = np.array(hdu[0].data)

FG=FG*1e3

FG_RSD = FG + RSD
FG_LC = FG + LC

names = np.loadtxt(base + 'filenames.txt', dtype='str')

kbins = 12
mubins = 15
ps_dist_array = np.zeros([kbins,len(names),2])
ps_dist_array_mu = np.zeros([mubins,kbins,len(names),2])
ps_array_1d_NOISE = np.zeros([kbins,len(names)])
ps_array_mu_NOISE = np.zeros([mubins,kbins,len(names)])

i0 = 50
i1 = 151

freq = np.arange(75.0, 95.1, 0.1)

L_para, L_perp = get_lengths(freq[i0], freq[i1], (freq[i0]+freq[i1])/2, 1.0)
box_dims = [L_perp, L_perp, L_para]
V_eor = L_perp*L_perp*L_para

for nn in range(0,len(names)):
    print("opening: "+base+"normalised_CD_1000hrs/"+names[nn])
    hdulist = fits.open(base+"normalised_CD_1000hrs/"+names[nn])
    NOISE = hdulist[0].data
    print("opened: "+base+"normalised_CD_1000hrs/"+names[nn])

    NOISE = NOISE *1e3
    NOISE = np.transpose(NOISE, (1,2,0))

    FG_RSD_NOISE = FG_RSD + NOISE
    FG_LC_NOISE = FG_LC + NOISE

    model, resids_RSD_FG_NOISE = do_fastica(FG_RSD_NOISE, 4)
    model, resids_LC_FG_NOISE = do_fastica(FG_LC_NOISE, 4)

    p, k = t2c.power_spectrum_1d(resids_RSD_FG_NOISE[:,:,i0:i1], kbins=kbins, box_dims=box_dims, binning ='log', return_n_modes=False)
    ps_dist_array[:,nn,0] = (p*V_eor*k**3)/(2*np.pi**2)
    p, k = t2c.power_spectrum_1d(resids_LC_FG_NOISE[:,:,i0:i1], kbins=kbins, box_dims=box_dims, binning ='log', return_n_modes=False)
    ps_dist_array[:,nn,1] = (p*V_eor*k**3)/(2*np.pi**2)

    p, mu, k = t2c.power_spectrum_mu(resids_RSD_FG_NOISE[:,:,i0:i1], los_axis = 2, box_dims=box_dims, mubins=mubins,kbins=kbins, exclude_zero_modes=True,return_n_modes=False,absolute_mus=False)
    ps_dist_array_mu[:,:,nn,0] = (p*V_eor*k**3)/(2*np.pi**2)
    p, mu, k = t2c.power_spectrum_mu(resids_LC_FG_NOISE[:,:,i0:i1], los_axis = 2, box_dims=box_dims, mubins=mubins,kbins=kbins, exclude_zero_modes=True,return_n_modes=False,absolute_mus=False)
    ps_dist_array_mu[:,:,nn,1] = (p*V_eor*k**3)/(2*np.pi**2)

    p, k = t2c.power_spectrum_1d(NOISE[:,:,i0:i1], kbins=kbins, box_dims=box_dims, binning ='log', return_n_modes=False)
    ps_array_1d_NOISE[:,nn] = (p*V_eor*k**3)/(2*np.pi**2)

    p, mu, k = t2c.power_spectrum_mu(NOISE[:,:,i0:i1], los_axis = 2, box_dims=box_dims, mubins=mubins,kbins=kbins, exclude_zero_modes=True,return_n_modes=False,absolute_mus=False)
    ps_array_mu_NOISE[:,:,nn] = (p*V_eor*k**3)/(2*np.pi**2)


errors = np.zeros([kbins,3])

for ii in range(0,kbins):

    errors[ii,0] = np.std(ps_dist_array[ii,:,0])
    errors[ii,1] = np.std(ps_dist_array[ii,:,1])
    errors[ii,2] = np.std(ps_array_1d_NOISE[ii,:])

np.savetxt(base+"errors_on_FastICA_CD_1Dps.txt",errors)

errors_mu = np.zeros([mubins,kbins,3])
for ii in range(0,kbins):
    for jj in range(0,mubins):
        errors_mu[jj,ii,0] = np.std(ps_dist_array_mu[jj,ii,:,0])
        errors_mu[jj,ii,1] = np.std(ps_dist_array_mu[jj,ii,:,1])
        errors_mu[jj,ii,2] = np.std(ps_array_mu_NOISE[jj,ii,:])

np.savetxt(base+"errors_on_FastICA_CD_mups.txt",errors_mu.reshape((mubins*kbins*3)))

