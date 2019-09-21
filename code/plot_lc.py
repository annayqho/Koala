import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from astropy.io import ascii
from astropy.cosmology import Planck15

z = 0.2714

dat = ascii.read("../data/ZTF18abvkwla_lct.csv")
filts = np.array([val[2:3] for val in dat['filter']])

fig,ax = plt.subplots(1,1, figsize=(7,5))

# Plot the g-band light curve
choose = np.logical_and(filts=='g', dat['mag'] > 0)
jd = dat['jdobs'][choose]
t0 = jd[0]
dt = (jd-t0)/(1+z)
rmag = dat['mag'][choose]
ermag = dat['mag_unc'][choose]
ax.errorbar(
        dt, rmag, yerr=ermag, 
        fmt='s', c='k', lw=0.5, label=r"ZTF $g$")

# Plot the g-band upper limits
choose = np.logical_and(filts=='g', dat['mag'] == 0)
jd = dat['jdobs'][choose]
dt = (jd-t0)/(1+z)
lims = dat['limmag'][choose]
ax.scatter(dt, lims, c='k', marker='s', label=None)
for ii,dt_val in enumerate(dt):
    ax.arrow(dt_val, lims[ii], 0, +0.1, color='k', length_includes_head=True,
            head_width=0.2, head_length=0.05, label=None)
 
# Plot the r-band light curve
choose = np.logical_and(filts=='r', dat['mag'] > 0)
jd = dat['jdobs'][choose]
dt = (jd-t0)/(1+z)
gmag = dat['mag'][choose]
egmag = dat['mag_unc'][choose]
ax.errorbar(
        dt, gmag, yerr=egmag, fmt='o', c='grey', mec='k', 
        mfc='white', lw=0.5, label=r"ZTF $r$")

# Plot the r-band upper limits
choose = np.logical_and(filts=='r', dat['mag'] == 0)
jd = dat['jdobs'][choose]
dt = (jd-t0)/(1+z)
lims = dat['limmag'][choose]
ax.scatter(dt, lims, facecolor='white', edgecolor='k', marker='o', label=None)
for ii,dt_val in enumerate(dt):
    ax.arrow(dt_val, lims[ii], 0, +0.1, edgecolor='k', lw=0.5,
            facecolor='white', length_includes_head=True,
            head_width=0.2, head_length=0.05, label=None)

# Format the panel
ax2 = ax.twinx()
ax2.set_ylabel(
        r"Absolute Magnitude",
        fontsize=16, rotation=270, labelpad=15.0)
y_f = lambda y_i: y_i-Planck15.distmod(z=0.2714).value
ymin, ymax = ax.get_ylim()
print(ymin)
print(y_f(ymin))
ax2.set_ylim((y_f(ymin), y_f(ymax)))
ax2.plot([],[])
ax2.tick_params(axis='both', labelsize=16)
ax2.invert_yaxis()

ax.set_xlabel("Rest-frame days since first detection", fontsize=16)
ax.set_ylabel("Apparent magnitude", fontsize=16)
ax.tick_params(axis='both', labelsize=16)
ax.invert_yaxis()
ax.legend(fontsize=16)
ax.set_xlim(-2.2,8.1)
ax.set_ylim(22,19.4)

plt.tight_layout()
plt.show()
#plt.savefig("lc.png", dpi=100)
