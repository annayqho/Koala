""" Use ZTFQuery to get basic info about the nightly cadence survey """

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from astropy.time import Time
from sqlalchemy import create_engine
import pandas as pd
import requests


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
            g_dt_all, histtype='step', range=(-0.5,10.5), bins=11,
            color='#140b34', weights=weights, label="$g$-band", lw=2)
    weights = np.ones_like(r_dt_all)/float(len(r_dt_all))
    plt.hist(
            r_dt_all, histtype='step', range=(-0.5,10.5), bins=11,
            color='#e55c30', weights=weights, label="$r$-band", lw=2)
    plt.xlabel("Nights Between Observations (Days)", fontsize=16)
    plt.ylabel("Fraction of Observations", fontsize=16)
    plt.tick_params(labelsize=14)
    plt.legend(fontsize=14, loc='upper right')
    plt.tight_layout()
    #plt.show()
    plt.savefig("cadence_hist.eps", format='eps', dpi=500)


def calc_efficiency(dt):
    n_prevday = np.sum(np.logical_and(dt > 0.5, dt < 1.5), axis=0) > 0
    n_2d = np.sum(np.logical_and(dt > 1.5, dt < 2.5), axis=0) > 0
    n_3in5 = np.sum(np.logical_and(dt > 0, dt < 5.5), axis=0) >= 3
    good_time = np.logical_and.reduce((n_prevday, n_2d, n_3in5))
    eff = sum(good_time)#/len(good_time)
    return eff


def get_efficiency(mjd, fields, filt, fieldid):
    choose = np.logical_and(fields == fieldid, filt=='g')
    g_obsdates = mjd[choose].values
    g_dt = g_obsdates[None,:]-g_obsdates[:,None]
    g_eff = calc_efficiency(g_dt)

    choose = np.logical_and(fields == fieldid, filt=='r')
    r_obsdates = mjd[choose].values
    r_dt = r_obsdates[None,:]-r_obsdates[:,None]
    r_eff = calc_efficiency(r_dt)
    return g_eff, r_eff


def get_nobs(fields, filt, fieldid):
    choose = np.logical_and(fields == fieldid, filt=='g')
    g_nobs = sum(choose)

    choose = np.logical_and(fields == fieldid, filt=='r')
    r_nobs = sum(choose)
    return g_nobs, r_nobs


def plot_efficiencies():
    nflds = len(fieldids)
    geff = np.zeros(nflds)
    reff = np.zeros(nflds)

    for ii,fid in enumerate(fieldids):
        gval, rval = get_efficiency(mjd, fields, filt, fid)
    #    gval, rval= get_nobs(fields, filt, fid)
        if fid==554:
            print(gval, rval)
        geff[ii] = gval
        reff[ii] = rval

    print(sum(geff))
    print(sum(reff))

    plt.hist(
            reff, histtype='step', color='#e55c30', lw=2, 
            bins=11)
            #range=(-0.05, 1.05), bins=11, label="$r$")

    plt.hist(
            geff, histtype='step', color='#140b34', lw=2, 
            bins=11)
            #range=(-0.05, 1.05), bins=11, label="$g$")

    plt.tick_params(axis='both', labelsize=16)
    plt.xlabel("Number of Obs", fontsize=16)
    plt.ylabel("Number of Fields", fontsize=16)
    plt.legend(fontsize=14)
    plt.tight_layout()
    #plt.savefig("nobs.eps", format='eps', dpi=500)
    plt.show()


def get_footprint(fields,nights,night):
    """ for a given night, get the footprint """
    choose = nights == night
    fids = list(np.unique(fields[choose]))
    fids = [str(val) for val in fids]
    fstr = "%2C%0D%0A+".join(fids)
    url = "http://yupana.caltech.edu/cgi-bin/ptf/tb//zoc?skycoverage=1&fieldlist=%s&submitskycoverage=SUBMIT+Sky+Coverage" %fstr
    out = requests.get(url)
    toparse = out.text.split('\n')
    footprint = float((toparse[6].strip()).split("=")[-1])
    return footprint


# These are all the 1DC observations
path = '/Users/annaho/Dropbox/Projects/Research/Koala/docs/anna.db'
engine = create_engine('sqlite:////Users/annaho/Dropbox/Projects/Research/Koala/data/anna.db')

df = pd.read_sql('Summary', engine)

# Eric says: note that the limiting mags are not meaningful

# Load in the relevant info
mjd = df['expMJD']
fields = df['fieldID']
filt = df['filter']
nights = df['night']
fieldids = np.unique(fields)
unights = np.unique(nights)

coverage = []

for night in unights:
    print(night)
    area = get_footprint(fields,nights,night)
    print(area)
    coverage.append(area)

# np.savetxt("fieldlist.txt", fieldids, fmt='%s')

# Make the cadence histogram
# cadence_hist(mjd, fields, filt)
