""" Radio analysis of rapidly evolving transients
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", family="serif")
rc("text", usetex=True)
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck15
from astropy.time import Time
import time

def get_transients():
    """ Get data on the rapidly evolving transients """
    data_dir = "/Users/annaho/Dropbox/Projects/Research/AT2018cow/data"
    drout = Table.read(
        "%s/drout_transients.txt" %data_dir, 
        format="ascii.no_header", delimiter="&")
    des = Table.read(
        "%s/des_transients.txt" %data_dir, 
        format="ascii.no_header", delimiter="&")
    des_dates_gs = Table.read(
        "%s/des_transients_dates.txt" %data_dir, 
        format="ascii.no_header", delimiter="&")
    des_dates_b = Table.read(
        "%s/des_transients_bronze_dates.txt" %data_dir, 
        format="ascii.no_header", delimiter=" ")
    ksn = np.array(["13:31:51.64", "-10:44:09.48", 0.090])

    ra_raw = np.hstack((drout['col3'], des['col2'], ksn[0]))
    dec_raw = np.hstack((drout['col4'], des['col3'], ksn[1]))
    z_raw = np.hstack((drout['col5'], des['col4'], ksn[2]))
    names_raw = np.hstack((drout['col1'], des['col1'], 'KSN'))

    names_mjd = np.hstack((des_dates_gs['col1'], des_dates_b['col1']))
    mjd = np.hstack((des_dates_gs['col2'], des_dates_b['col2']))

    dates = []
    for date_raw in drout['col2']:
        split = date_raw.split()
        yyyy = split[0]
        mm = split[1].strip('\.')
        mm = mm.strip('\\')
        # get month as a number
        try:
            t = str(time.strptime(mm, '%B').tm_mon).zfill(2)
        except:
            t = str(time.strptime(mm, '%b').tm_mon).zfill(2)
        dd = split[2]
        date = Time("%s-%s-%s" %(yyyy,t,dd), format='iso')
        dates.append(date.mjd)
    for name in des['col1']:
        ind = np.where(names_mjd==name)[0][0]
        dates.append((Time(mjd[ind], format='mjd')).mjd)
    dates.append((Time(2457232.70, format='jd')).mjd)
    dates = np.array(dates)

    ra = np.array([val for val in ra_raw])

    # curate the RA and Dec list
    c = SkyCoord(ra_raw, dec_raw, frame='icrs', unit=(u.hourangle, u.deg))

    # curate the redshift list
    z = []
    for ii, zval in enumerate(z_raw):
        try:
            z.append(float(zval))
        except:
            z.append(-1)
    z = np.array(z)
    return names_raw, ra_raw, dec_raw, dates, z


if __name__=="__main__":
    # Retrieve data
    names, ra, dec, dates, z = get_transients()

    # initialize plot
    fig,axarr = plt.subplots(2, 1, figsize=(5,6))

    # redshift distribution of all of them
    choose = z > 0
    axarr[0].hist(z[choose], histtype='step', color='k')
    ax = axarr[0]
    ax.set_xlabel("Redshift", fontsize=14)
    ax.set_ylabel("\# Transients", fontsize=14)
    ax.tick_params(axis='both', labelsize=14)
    ax.axvline(x=0.014, c='k', ls='--', lw=2.0)

    # distribution of peak flux at 230 GHz
    dist = Planck15.luminosity_distance(z=z[choose]).cgs.value
    ref = Planck15.luminosity_distance(z=0.014).cgs.value
    flux = 50e3 * (ref / dist)**2 
    ax = axarr[1]
    ax.hist(
            flux, histtype='step', color='k', 
            bins=np.logspace(np.log10(1.1), np.log10(5000), 10))
    ax.set_xscale('log')
    ax.set_xlabel("Peak Flux ($\\mu$Jy) at $230\,$GHz", fontsize=14)
    ax.set_ylabel("\# Transients", fontsize=14)
    ax.tick_params(axis='both', labelsize=14)
    ax.axvline(x=50e3, c='k', ls='--', lw=2.0)
    ax.text(
            0.94, 0.9, "AT2018cow", transform=ax.transAxes,
            horizontalalignment='right', fontsize=14)
    ax.text(
            0.94, 0.8, "$50\,$mJy", transform=ax.transAxes,
            horizontalalignment='right', fontsize=14)

    # formatting
    plt.tight_layout()
    plt.show()
