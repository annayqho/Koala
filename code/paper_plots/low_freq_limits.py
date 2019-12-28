""" Get low-frequency VLASS limits """

import numpy as np
from matplotlib import rc
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/AT2018cow/code")
rc("font", family="serif")
rc("text", usetex=True)
import matplotlib.pyplot as plt
import sys
from matplotlib.ticker import MaxNLocator
from astropy.io import ascii
from astropy.cosmology import Planck15
from get_radio import get_data_all
from synchrotron_fit import self_abs,fit_self_abs


def plot_18cow_atca():
    """ construct a 3 GHz and 9 GHz light curve based on the ATCA data """
    dcm = Planck15.luminosity_distance(z=0.014).cgs.value

    # Step 1: get all of the days with ATCA data
    tel, freq, days, flux, eflux_form, eflux_sys = get_data_all()
    ef = np.sqrt(eflux_form**2 + eflux_sys**2)
    # there are 7 days with ATCA data
    # of those 7 days, the first three have multiple points
    # and the last four only have 34 GHz

    # Plot the 9 GHz data
    choose = np.logical_and(tel=='ATCA', freq==9)
    lum = 9E9 * flux[choose] * 1E-3 * 1E-23 * 4 * np.pi * dcm**2
    elum = 9E9 * ef[choose] * 1E-3 * 1E-23 * 4 * np.pi * dcm**2
    plt.errorbar(days[choose], lum, yerr=elum, c='blue', fmt='s')

    # Step 2: for each of those days, fit nu squared and extrapolate
    # down to 3 GHz, to generate a 3 GHz light curve
    dt = []
    lc = []
    udays = np.unique(days[tel=='ATCA'])
    for day in udays:
        dt.append(day)
        choose = np.logical_and(days == day, eflux_form != 99) # no upper limits
        popt, pcov = fit_self_abs(freq[choose], flux[choose], ef[choose])
        lc.append(self_abs(3, popt[0]))
    dt = np.array(dt)
    lc = np.array(lc)

    # Step 3: Plot the 3 GHz light curve as luminosity
    # give them all 20% uncertainties
    lum = 3E9 * lc * 1E-3 * 1E-23 * 4 * np.pi * dcm**2
    # but only t <= 40 days or so
    choose = dt <= 40
    plt.errorbar(dt[choose], lum[choose], yerr=0.2*lum[choose], 
            c='k', fmt='.', ms=10, 
            label='AT2018cow, 3 GHz')


def plot_18cow_margutti():
    # Now, plot the Margutti (2018) low-freqeuncy data as a luminosity...
    # choose the 3.5 GHz point, which is at 83.52 days
    dcm = Planck15.luminosity_distance(z=0.014).cgs.value
    lum = 3.5E9 * 3.24 * 1E-3 * 1E-23 * 4 * np.pi * dcm**2
    elum = 3.5E9 * 0.06 * 1E-3 * 1E-23 * 4 * np.pi * dcm**2
    plt.errorbar(83.52, lum, yerr=elum, c='k', fmt='.', ms=10)

    # now their 9 GHz data
    lum = 9E9 * 9.10 * 1E-3 * 1E-23 * 4 * np.pi * dcm**2
    elum = 9E9 * 0.30 * 1E-3 * 1E-23 * 4 * np.pi * dcm**2
    plt.errorbar(
            83.51, lum, yerr=elum, c='blue', fmt='s', 
            label="AT2018cow, 9 GHz")


def plot_koala():
    # Plot the X-band Koala detections
    f = np.array([8.3E39,5.85E38])
    plt.errorbar([81,313], f, 0.0006*f, c='blue', fmt='*', ms=20)
    plt.plot([81,313], f, c='blue', lw=2)
    plt.text(
            70, 1E40, "ZTF18abvkwla", fontsize=14, 
            horizontalalignment='right', verticalalignment='center')


def plot_98bw():
    dat = np.loadtxt("/Users/annaho/Dropbox/Projects/Research/IcBL/data/radio_compilations/98bw_low_freq.txt", delimiter=',')
    dcm = Planck15.luminosity_distance(z=0.0085).cgs.value
    lum = dat[:,1] * 1E-3 * 1E-23 * 4 * np.pi * dcm**2 * 1.38E9
    plt.plot(dat[:,0], lum, c='grey', label='98bw, 1.38 GHz')


plot_18cow_atca()
plot_18cow_margutti()
plot_koala()
plot_98bw()

# Step 4: Add in the VLASS limits
dat = ascii.read("../data/radio_lims.dat")
flim_raw = dat['limit']
z_raw = dat['z']
bad = np.logical_or(flim_raw.mask, z_raw <= 0)
dt = dat['dt'][~bad]
flim = flim_raw[~bad]
z = z_raw[~bad]

dcm = Planck15.luminosity_distance(z=z).cgs.value
llim = 3E9*flim*1E-6*1E-23*4*np.pi*dcm**2
plt.scatter(dt, llim, marker='v', s=30, c='k',
        label='VLASS limits, 2--4 GHz')

# Formatting
plt.tick_params(axis='both', labelsize=14)
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r"$\nu L_{\nu}$ [erg/s]", fontsize=16)
plt.xlabel("$\Delta t$ [days]", fontsize=16)
plt.legend(fontsize=11, loc='lower right', ncol=2)
plt.tight_layout()

#plt.show()
plt.savefig("vlass_lims.png", dpi=200)
