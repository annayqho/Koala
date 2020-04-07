""" put a limit on the amount of nickel that could be present """

from astropy.cosmology import Planck15
import numpy as np

c = 3E10
dcm = Planck15.luminosity_distance(z=0.2714).cgs.value
z = 0.2714

gmag = 21.51
egmag = 0.21
dt = 10.07
dt_rest = dt/(1+z)
lam = 4722.7E-8 # cm
lam_rest = lam/(1+z)

fnu = 10**((gmag+48.6)/(-2.5))
# In rest-frame
fnu_rest = fnu/(1+z)
flam = fnu_rest * c / lam_rest**2 
Flam = flam * 4 * np.pi * dcm**2
lum = lam_rest * Flam
print(lum)

tni = 8.8
tco = 113.6
mni = lum / (2E43) / (
        3.9*np.exp(-dt_rest/tni)+0.678*(np.exp(-dt_rest/tco)-np.exp(-dt_rest/tni)))
print(mni)

# from Khatami's equation
beta = 4/3
mni = (lum*beta**2*8**2) / (2*2E43*tni**2) / (1-(1+beta*8/tni)*np.exp(-beta*8/tni))
