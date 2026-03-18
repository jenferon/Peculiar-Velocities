import numpy as np
import struct
import astropy.io.fits as fits
import tools21cm as t2c
from sklearn.decomposition import FastICA
from scipy import interpolate
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

bins = 8

ps_dist_array = np.zeros([bins,bins,len(names),2])

i0 = 50
i1 = 151

freq = np.arange(75, 95.1, 0.1)

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

    pp_NOISE_LC, kper_NOISE_LC, kpar_NOISE_LC, n_NOISE_LC = t2c.power_spectrum_2d(resids_LC_FG_NOISE[:,:,i0:i1], kbins=bins, box_dims= box_dims, return_modes=True)
    for ii in range(0,len(kper_NOISE_LC)):
        for jj in range(0,len(kpar_NOISE_LC)):
            pp_NOISE_LC[ii,jj] = (pp_NOISE_LC[ii,jj]*V_eor*np.sqrt(kper_NOISE_LC[ii]**2 + kpar_NOISE_LC[jj]**2)**3)/(2*np.pi**2) 
    fp_NOISE_LC = interpolate.interp2d(kper_NOISE_LC, kpar_NOISE_LC, pp_NOISE_LC.T, kind='linear')
    X_NOISE, Y_NOISE = kper_NOISE_LC, kpar_NOISE_LC
    ps_dist_array[:,:,nn,0] = fp_NOISE_LC(X_NOISE,Y_NOISE)

    pp_NOISE_RSD, kper_NOISE_RSD, kpar_NOISE_RSD, n_NOISE_RSD = t2c.power_spectrum_2d(resids_RSD_FG_NOISE[:,:,i0:i1], kbins=bins, box_dims= box_dims, return_modes=True)
    for ii in range(0,len(kper_NOISE_RSD)):
        for jj in range(0,len(kpar_NOISE_RSD)):
            pp_NOISE_RSD[ii,jj] = (pp_NOISE_RSD[ii,jj]*V_eor*np.sqrt(kper_NOISE_RSD[ii]**2 + kpar_NOISE_RSD[jj]**2)**3)/(2*np.pi**2) 
    fp_NOISE_RSD = interpolate.interp2d(kper_NOISE_RSD, kpar_NOISE_RSD, pp_NOISE_RSD.T, kind='linear')
    X_NOISE, Y_NOISE = kper_NOISE_RSD, kpar_NOISE_RSD
    ps_dist_array[:,:,nn,1] = fp_NOISE_RSD(X_NOISE,Y_NOISE)

errors = np.zeros([bins,bins,2])

for ii in range(0,bins):
    for jj in range(0,bins):
        errors[ii,jj,0] = np.std(ps_dist_array[ii,jj,:,0])
        errors[ii,jj,1] = np.std(ps_dist_array[ii,jj,:,1])

np.savetxt(base+"errors_on_FastICA_CD_cyclindricalps.txt",errors.reshape((bins*bins*2)))
