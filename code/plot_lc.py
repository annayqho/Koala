import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.cosmology import Planck15

dat = ascii.read("../data/ZTF18abvkwla_lct.csv")
filts = np.array([val[2:3] for val in dat['filter']])

fig,axarr = plt.subplots(2,1, figsize=(5,6))

# Plot the g-band light curve
ax = axarr[0]
choose = np.logical_and(filts=='g', dat['mag'] > 0)
jd = dat['jdobs'][choose]
t0 = jd[0]
dt = jd-t0
rmag = dat['mag'][choose]
ermag = dat['mag_unc'][choose]
ax.errorbar(dt, rmag, yerr=ermag, fmt='s', c='k', label="g-band")

# Plot the r-band light curve
choose = np.logical_and(filts=='r', dat['mag'] > 0)
jd = dat['jdobs'][choose]
dt = jd-t0
gmag = dat['mag'][choose]
egmag = dat['mag_unc'][choose]
ax.errorbar(
        dt, gmag, yerr=egmag, fmt='o', c='grey', mec='k', 
        mfc='white', label="r-band")

# Format the top panel
ax.set_ylabel("Days since first g-band detection")
ax.set_ylabel("Apparent magnitude")
ax.invert_yaxis()
ax.legend()

# Plot the g-r color over time
ax = axarr[1]
dt_gr = []
abs_g = []
gr = []
for dt_val in np.arange(0,10):
    choose = np.logical_and(dt > dt_val-0.1, dt < dt_val+0.1)
    ng = sum(gmag[choose])
    g = gmag[dt_val]
    r = rmag[dt_val]
    if ng > 1:
        g = np.average(gmag[choose])
    nr = sum(rmag[choose])
    if nr > 1:
        r = np.average(rmag[choose])
    if np.logical_and(nr>0, ng>0):
        gr.append(-(g-r))
        dt_gr.append(dt_val)
        abs_g.append(g-Planck15.distmod(z=0.27).value)
cax = ax.scatter(gr, abs_g, marker='o', c=dt_gr)
fig.colorbar(cax, ax=ax)

# Format the bottom panel
ax.set_xlabel("$g-r$")
ax.set_ylabel("Absolute $g$-band mag")
ax.invert_yaxis()

plt.tight_layout()
plt.show()
