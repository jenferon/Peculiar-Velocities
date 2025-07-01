import py21cmsense as p21sense
from py21cmsense import GaussianBeam, Observation, Observatory, PowerSpectrum, beam
import numpy as np
from astropy.cosmology import Planck13 as cosmoP
from astropy.cosmology import FlatLambdaCDM, LambdaCDM

import astropy.units as u
cosmo = FlatLambdaCDM(H0=71 * u.km / u.s / u.Mpc, Om0=0.27)

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

nu = 145
c = 2.998e8
eta_a = 1.0 #antenna efficiency
pi = 3.14159
k_b = 1.38e-23
eta_s=0.9 #system efficiency
T_recv = 40.0 # http://arxiv.org/pdf/0901.3359v1.pdf (140-180K)
zs = zFromNu(nu)

baseline_data = np.loadtxt("/home/ppxjf3/RSD_LC/SKA_all_stations_Rev3.tm/layout_ecef.txt")
"""aa4 =  Observatory(antpos=baseline_data* u.m, beam=beam.GaussianBeam(
           frequency=nu*u.MHz, dish_size=35*u.m), latitude = -27.0*u.deg, Trcv=T_recv*u.K) #check beam and dish size

#observation
observation_params = {}
observation_params["ndays"] = 166.7
observation_params["cosmo"] = cosmo
observation_params["h"] = cosmo.H0.value / 100.0
observation_params["freq_bands"] = nu
observation_params["redshifts"] = zs
observation_params["time_per_day_hrs"] = 6.0
observation_params["bandwidth"] = 10e6

observation = p21sense.Observation(
    observatory=aa4,
    lst_bin_size=observation_params["time_per_day_hrs"] * u.hour,
    time_per_day=observation_params["time_per_day_hrs"] * u.hour,
    n_days=observation_params["ndays"],
    bandwidth=observation_params["bandwidth"]*u.Hz,
    coherent=False, #add baselines coherantly or not
    cosmo=observation_params["cosmo"],
    tsky_amplitude = 17.1*u.K,
    tsky_ref_freq = 408.0*u.MHz,
    #baseline_filters=p21sense.BaselineRange(bl_max=1000* u.m) 
    )
xx=np.logspace(-1., 2.5, 30)
print(xx)
kperp_edges= u.Quantity(xx, "littleh/Mpc")
ska_aa4_senses1 =PowerSpectrum(foreground_model="foreground_free", horizon_buffer=0.0/u.Mpc,
            observation=observation,
        ).at_frequency(nu * u.MHz)
sense1d_both = ska_aa4_senses1.calculate_sensitivity_1d_binned(thermal=True, sample=False, k=kperp_edges)

array = np.zeros([2,len(xx)])
array[0,:] = xx*0.71
array[1,:] = sense1d_both
np.savetxt("/home/ppxjf3/RSD_LC/eor_sensitivity_full_baselines.txt", array)"""


aa4 =  Observatory(antpos=baseline_data* u.m, beam=beam.GaussianBeam(
           frequency=85*u.MHz, dish_size=35*u.m), latitude = -27.0*u.deg, Trcv=T_recv*u.K) #check beam and dish size

#observation
observation_params = {}
observation_params["ndays"] = 166.7
observation_params["cosmo"] = cosmo
observation_params["h"] = cosmo.H0.value / 100.0
observation_params["freq_bands"] = 85
observation_params["redshifts"] = zFromNu(85)
observation_params["time_per_day_hrs"] = 6.0
observation_params["bandwidth"] = 10e6

observation = p21sense.Observation(
    observatory=aa4,
    lst_bin_size=observation_params["time_per_day_hrs"] * u.hour,
    time_per_day=observation_params["time_per_day_hrs"] * u.hour,
    n_days=observation_params["ndays"],
    bandwidth=observation_params["bandwidth"]*u.Hz,
    coherent=False, #add baselines coherantly or not
    cosmo=observation_params["cosmo"],
    tsky_amplitude = 17.1*u.K,
    tsky_ref_freq = 408.0*u.MHz,
    #baseline_filters=p21sense.BaselineRange(bl_max=8000* u.m) 
    )
xx=np.logspace(-1., 2.5, 30)
print(xx)
kperp_edges= u.Quantity(xx, "littleh/Mpc")
ska_aa4_senses1 =PowerSpectrum(foreground_model="foreground_free", horizon_buffer=0.0/u.Mpc,
            observation=observation,
        ).at_frequency(85 * u.MHz)
sense1d_both = ska_aa4_senses1.calculate_sensitivity_1d_binned(thermal=True, sample=False, k=kperp_edges)

array = np.zeros([2,len(xx)])
array[0,:] = xx*0.71
array[1,:] = sense1d_both
np.savetxt("/home/ppxjf3/RSD_LC/cd_sensitivity_full_baselines.fits", array)