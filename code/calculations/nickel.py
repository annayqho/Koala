""" put a limit on the amount of nickel that could be present """

from astropy.cosmology import Planck15
import numpy as np

c = 3E10
dcm = Planck15.luminosity_distance(z=0.2714).cgs.value

rmag = 21.55
#rmag = 19.4

# r-band lambda is 6231 (Sloan r)
# but in rest-frame 4950
lam = 4950E-8 # cm

# g-band lambda (rest-frame)
#lam = 3820E-8

# r-band mag = 21.55 +/- 0.23

# -19.5 mag is 
fnu = 10**((rmag+48.6)/(-2.5))
print(fnu)

# which is
flam = fnu * c / lam**2 

# at the distance of the source
Flam = flam * 4 * np.pi * dcm**2

# luminosity
lum = lam * Flam
print(lum)

tni = 8.8
tco = 113.6
t = 8
mni = lum / (2E43) / (3.9*np.exp(-t/tni)+0.678*(np.exp(-t/tco)-np.exp(-t/tni)))
print(mni)

# from Khatami's equation
beta = 4/3
mni = (lum*beta**2*8**2) / (2*2E43*tni**2) / (1-(1+beta*8/tni)*np.exp(-beta*8/tni))
