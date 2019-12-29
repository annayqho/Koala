import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from astropy.io import ascii
from astropy.cosmology import Planck15

z = 0.2714

dat = ascii.read("../../data/ZTF18abvkwla_lct.csv")
filts = np.array([val[2:3] for val in dat['filter']])

fig,ax = plt.subplots(1,1, figsize=(7,5))

# Plot the g-band light curve
choose = np.logical_and(filts=='g', dat['mag'] > 0)
jd = dat['jdobs'][choose]
t0 = jd[0]
dt = (jd-t0)/(1+z)
rmag = dat['mag'][choose] - 0.167 #extinction in SDSS g
ermag = dat['mag_unc'][choose]
ax.errorbar(
        dt, rmag, yerr=ermag, 
        fmt='s', c='k', lw=0.5, 
        label=r"ZTF $g$ ($\lambda_\mathrm{rest}=3820$\AA)")

# Plot the g-band upper limits
choose = np.logical_and(filts=='g', dat['mag'] == 0)
jd = dat['jdobs'][choose]
dt = (jd-t0)/(1+z)
lims = dat['limmag'][choose] - 0.167 # extinction in SDSS g
ax.scatter(dt, lims, c='k', marker='s', label=None)
for ii,dt_val in enumerate(dt):
    ax.arrow(dt_val, lims[ii], 0, +0.1, color='k', length_includes_head=True,
            head_width=0.2, head_length=0.05, label=None)
 
# Plot the r-band light curve
choose = np.logical_and(filts=='r', dat['mag'] > 0)
jd = dat['jdobs'][choose]
dt = (jd-t0)/(1+z)
rmag = dat['mag'][choose] - 0.115
ermag = dat['mag_unc'][choose]
ax.errorbar(
        dt, rmag, yerr=ermag, fmt='o', c='grey', mec='k', 
        mfc='white', lw=0.5, 
        label=r"ZTF $r$ ($\lambda_\mathrm{rest}=4950$\AA)")

# Plot the r-band upper limits
choose = np.logical_and(filts=='r', dat['mag'] == 0)
jd = dat['jdobs'][choose]
dt = (jd-t0)/(1+z)
lims = dat['limmag'][choose] - 0.115
ax.scatter(dt, lims, facecolor='white', edgecolor='k', marker='o', label=None)
for ii,dt_val in enumerate(dt):
    ax.arrow(dt_val, lims[ii], 0, +0.1, edgecolor='k', lw=0.5,
            facecolor='white', length_includes_head=True,
            head_width=0.2, head_length=0.05, label=None)

# Plot the 18cow light curves
# scale to the distance of Koala
dmag = Planck15.distmod(z=0.2714).value-Planck15.distmod(0.0141).value

# what filter is closest to 4950 rest-frame?
# for LT, it's in g-band
t_g = np.array([58287.1500, 58288.2020, 58289.0234, 58289.1889, 58289.9147,
    58290.1082, 58290.2189, 58290.2922, 58290.9029, 58291.0145, 58291.1742])
t0_cow = t_g[0] - 2
dt_g = t_g - t0_cow
mag_g = np.array([13.40, 13.65, 14.10, 14.18, 14.48, 
    14.57, 14.63, 14.70, 14.94, 14.99, 15.09]) - 0.287 # extinction in g
ax.plot(dt_g, mag_g + dmag, c='k', lw=0.5)
ax.text(4.5, 21.1, "18cow $g$", fontsize=11, rotation=-55)

# Connect the LC to the last upper limit 
# doesn't look good to do this...
#ax.plot(
#        [58284.1300-t0_cow, dt_g[0]], 
#        [18.90 + dmag, mag_g[0]+dmag], c='k', lw=0.5, ls='--')
# instead, just indicate the last non-detection
ax.plot([-1, -1], [21.6, 21.8], c='k', lw=0.5)
ax.text(-0.9, 21.68, 'time of last 18cow upper limit', fontsize=11)

# what filter is closest to 3820 rest-frame? I think u-band
# but it's sort of between u-band and b-band
# for LT, it's pretty squarely in u-band
t_u = np.array([58289.0227, 58289.1963, 58289.6434, 58289.9154, 58290.2243,
    58290.2267, 58290.6697, 58290.9036, 58291.1809])
mag_u = np.array(
        [13.97, 14.03, 14.24, 14.31, 14.43, 14.46, 14.63, 14.66, 14.83])
ax.plot(t_u-t0_cow, mag_u - 0.368 + dmag, c='k', lw=0.5, ls='-')
ax.text(5, 21.05, "18cow $u$", fontsize=11, rotation=-55)

# Format the panel
ax.set_ylim(21.7, 19.2)
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

ax.text(0.9, 19.32, 'S', fontsize=14)

ax.set_xlabel("Rest-frame days since first detection", fontsize=16)
ax.set_ylabel("Apparent magnitude", fontsize=16)
ax.tick_params(axis='both', labelsize=16)
ax.legend(fontsize=13)
ax.set_xlim(-2.2,8.1)

plt.tight_layout()
#plt.show()
plt.savefig("lc.png", dpi=300)
