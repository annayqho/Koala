""" Use ZTFQuery to get basic info about the nightly cadence survey """

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from astropy.time import Time
from sqlalchemy import create_engine
import pandas as pd


def cadence_hist(mjd, fields, filt):
    g_dt_all = []
    r_dt_all = []

    for fieldid in np.unique(fields):
        # Use g-band 
        choose = np.logical_and(fields == fieldid, filt=='g')
        g_obsdates = mjd[choose].values
        g_dt = g_obsdates[1:]-g_obsdates[0:-1]
        [g_dt_all.append(val) for val in g_dt]
        start = Time(min(g_obsdates), format='mjd').iso
        end = Time(max(g_obsdates), format='mjd').iso

        choose = np.logical_and(fields == fieldid, filt=='r')
        r_obsdates = mjd[choose].values
        r_dt = r_obsdates[1:]-r_obsdates[0:-1]
        [r_dt_all.append(val) for val in r_dt]
        start = Time(min(r_obsdates), format='mjd').iso
        end = Time(max(r_obsdates), format='mjd').iso

    g_dt_all = np.array(g_dt_all)
    r_dt_all = np.array(r_dt_all)

    plt.figure(figsize=(6,4))
    weights = np.ones_like(g_dt_all)/float(len(g_dt_all))
    plt.hist(
            g_dt_all, histtype='step', range=(0,10),
            color='#140b34', weights=weights, label="$g$-band", lw=2)
    weights = np.ones_like(r_dt_all)/float(len(r_dt_all))
    plt.hist(
            r_dt_all, histtype='step', range=(0,10),
            color='#e55c30', weights=weights, label="$r$-band", lw=2)
    plt.xlabel("Nights Between Observations (Days)", fontsize=16)
    plt.ylabel("Fraction of Observations", fontsize=16)
    plt.tick_params(labelsize=14)
    plt.legend(fontsize=14, loc='upper right')
    plt.tight_layout()
    #plt.show()
    plt.savefig("cadence_hist.eps", format='eps', dpi=500)


# These are all the 1DC observations
path = '/Users/annaho/Dropbox/Projects/Research/Koala/docs/anna.db'
engine = create_engine('sqlite:///%s' %path)

df = pd.read_sql('Summary', engine)

# Eric says: note that the limiting mags are not meaningful

# Load in the relevant info
mjd = df['expMJD']
print(Time(min(mjd), format='mjd').iso)
print(Time(max(mjd), format='mjd').iso)
fields = df['fieldID']
filt = df['filter']

fieldids = np.unique(fields)
np.savetxt("fieldlist.txt", fieldids, fmt='%s')

# Make the cadence histogram
# cadence_hist(mjd, fields, filt)

# Number of epochs as a function of position in the sky
