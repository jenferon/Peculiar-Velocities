import numpy as np
import foreground_sim as fs
import astropy.io.fits as fits
import matplotlib.pyplot as plt

base = '/gpfs01/home/ppxjf3/peculiar_vel/'
numin = 140
nfreq = 401
ngrid= 512
FoV = 1.0
nustep = 0.1

sim_test = fs.ForeGsim(nfreq=nfreq, numin=numin, nustep=nustep, ngrid=ngrid, imgsize=FoV, nexgal=1477487)

syn_fg = sim_test.gen_sync_map() #generates synchrotron emission foreground map

ff_fg = sim_test.gen_freefree_map() #free-free emission due to bremsstrahllung raditation in diffuse ionised galatic gas

egfg_test = sim_test.gen_exgal_map()

fg = syn_fg + ff_fg + egfg_test
#fg_egfg = syn_fg + ff_fg + egfg_test

fits.writeto(base + 'Fg_egfg_N'+str(ngrid)+'_FOV'+str(FoV)+'_'+str(numin)+'_'+str(np.round(numin+nustep*nfreq))+'MHz_'+str(nustep)+'MHz.fits', fg, overwrite=True)
#fits.writeto(base + 'Fg_egfg_N512_FOV1.0_55_95MHz_0.1MHz.fits', fg_egfg, overwrite=True)

#global_temp_fg = np.empty([fg.shape[2]]) 
#global_temp_fg_egfg = np.empty([fg_egfg.shape[2]]) 
#for ii in range(0,Fgs.shape[2]):
#    global_temp_fg[ii] = np.mean(fg[:,:,ii])
#    global_temp_fg_egfg[ii] = np.mean(fg_egfg[:,:,ii])

#freq = np.arange(140, 180.1, 0.1)
  

#plt.plot(freq, global_temp_fg, color='k')
#plt.plot(freq, global_temp_fg_egfg, color='b')
#plt.savefig(base+'foreground_check.pdf')
