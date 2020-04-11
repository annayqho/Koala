""" Calculate rise and fade timescales for light curves """

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import pandas as pd
from astropy.cosmology import Planck15
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF_fast_transient_search/code")
from forced_phot.run_forced_phot import get_forced_phot

# corrections
def correct(mpeak, trise, etrise, tfade, etfade, z, mw_ext):
    Mpeak = mpeak-mw_ext-Planck15.distmod(z=z).value
    trise_corr = trise / (1+z)
    etrise_corr = etrise / (1+z)
    tfade_corr = tfade / (1+z)
    etfade_corr = etfade / (1+z)
    return Mpeak, trise_corr, etrise_corr, tfade_corr, etfade_corr


def iptf15ul():
    z = 0.066
    mw_ext = 0

    dat = ascii.read("../../data/iptf15ul/iptf15ul.txt")
    mjd = dat['col2']
    filt = dat['col3']
    mag = dat['col4']
    emag = dat['col5']
    choose = np.logical_and(filt == 'g', mag != '>')

    xall = mjd[choose].astype(float)
    yall = mag[choose].astype(float)
    yerrall = emag[choose].astype(float)

    tpeak = xall[np.argmin(yall)]
    mpeak = np.min(yall)
    empeak = yerrall[np.argmin(yall)]

    x = xall[xall<=tpeak].data
    y = yall[xall<=tpeak].data
    ey = yerrall[xall<=tpeak].data
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order]

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = tpeak-np.interp(mpeak+0.75, ysamples[ii], x)

    rise = np.mean(trise)
    erise = np.std(trise)

    x = xall[xall>=tpeak].data 
    y = yall[xall>=tpeak].data 
    ey = yerrall[xall>=tpeak].data 
    order = np.argsort(y)
    x = x[order]
    y = y[order] 
    ey = ey[order] 

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(mpeak+0.75, ysamples[ii], x)-tpeak

    fade = np.mean(tfade)
    efade = np.std(tfade)

    Mpeak, trise_corr, etrise_corr, tfade_corr, etfade_corr = correct(
        mpeak, rise, erise, fade, efade, z, mw_ext)

    print("mpeak is %s +/- %s" %(Mpeak, empeak))
    print("fade time is %s +/- %s" %(tfade_corr,etfade_corr))
    print("rise time is %s +/- %s" %(trise_corr,etrise_corr))


def at2018cow():
    z = 0.0141
    dat = pd.read_fwf("../../data/at2018cow_photometry_table.dat")
    mjd = dat['MJD']
    filt = dat['Filt']
    mag = dat['ABMag']
    emag = dat['Emag']
    choose = filt == 'g'

    xall = mjd[choose][1:].values
    yall = mag[choose][1:].values.astype(float)-0.287
    yerrall = emag[choose][1:].values.astype(float)

    # Add the detection prior to peak
    # o-band dominated by r-band flux
    # r-band extinction: 0.198
    xall = np.insert(xall, 0, 58285.44)
    yall = np.insert(yall, 0, 14.70-0.198-0.4)
    yerrall = np.insert(yerrall, 0, 0.10)

    # Add the i-band detection prior to peak
    xall = np.insert(xall, 0, 58286.1950)
    yall = np.insert(yall, 0, 14.32-0.147-0.7)
    yerrall = np.insert(yerrall, 0, 0.03)

    tpeak = xall[np.argmin(yall)]
    empeak = yerrall[np.argmin(yall)]
    mpeak = np.min(yall)
    print(mpeak-Planck15.distmod(z=z).value)
    print(empeak)

    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order] 
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(mpeak+0.75, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise),np.std(trise)))

    x = xall[xall>=tpeak]
    y = yall[xall>=tpeak]
    ey = yerrall[xall>=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order] 
    ey = ey[order] 

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(mpeak+0.75, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade),np.std(tfade)))


def dougie():
    redshift = 0.19
    mw_ext = 0.031

    dat = ascii.read("../../data/dougie.txt")
    xall = dat['col1']
    yall = dat['col2']
    yerrall = np.sqrt(dat['col3']**2+0.1**2)

    ind = np.argmin(yall)
    tpeak = xall[ind]
    mpeak = yall[ind]
    empeak = yerrall[ind]

    # The shape of the LC is weird so have to do this manually
    x = xall[xall<=tpeak].data[2:4]
    y = yall[xall<=tpeak].data[2:4]
    ey = yerrall[xall<=tpeak].data[2:4]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order]

    nsim = 1000
    trise= np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = tpeak-np.interp(mpeak+0.75, ysamples[ii], x)

    rise = np.mean(trise)
    erise = np.std(trise)

    x = xall[xall>=tpeak].data[4:6]
    y = yall[xall>=tpeak].data[4:6]
    ey = yerrall[xall>=tpeak].data[4:6]
    order = np.argsort(y)
    x = x[order]
    y = y[order] 
    ey = ey[order] 

    nsim = 1000
    tfade= np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(mpeak+0.75, ysamples[ii], x)-tpeak

    fade = np.mean(tfade)
    efade = np.std(tfade)

    Mpeak, trise_corr, etrise_corr, tfade_corr, etfade_corr = correct(
        mpeak, rise, erise, fade, efade, redshift, mw_ext)

    print("peak is %s+/-%s" %(Mpeak,empeak))
    print("fade time is %s +/- %s" %(tfade_corr,etfade_corr))
    print("rise time is %s +/- %s" %(trise_corr,etrise_corr))


def iptf16asu():
    z = 0.187
    # no need to correct for extinction
    mw_ext = 0

    dat = pd.read_fwf("../../data/iptf16asu.txt")
    mjd = dat['MJD']
    filt = dat['filt']
    mag = dat['mag']
    emag = dat['emag']

    # Rising behavior
    islim = np.array(['>' in m for m in mag])
    choose = np.logical_and(filt == 'g', islim==False)

    xall = mjd[choose].values
    yall = mag[choose].values.astype(float)
    yerrall = emag[choose].values.astype(float)

    ind = np.argmin(yall)
    tpeak = xall[ind]
    mpeak = yall[ind]
    empeak = yerrall[ind]
    mpeak = np.min(yall)
    Mpeak = mpeak-Planck15.distmod(z=z).value
    print("peak is %s +/- %s" %(Mpeak, empeak))

    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order] 
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(mpeak+0.75, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise)/(1+z),np.std(trise)/(1+z)))


    # Fading behavior
    choose = filt == 'r'

    xall = mjd[choose][0:3].values
    yall = mag[choose][0:3].values.astype(float)
    yerrall = emag[choose][0:3].values.astype(float)

    ind = np.argmin(yall)
    tpeak = xall[ind]
    mpeak = yall[ind]
    empeak = yerrall[ind]
    mpeak = np.min(yall)
    Mpeak = mpeak-Planck15.distmod(z=z).value
    print("peak is %s +/- %s" %(Mpeak, empeak))

    x = xall[xall>=tpeak]
    y = yall[xall>=tpeak]
    ey = yerrall[xall>=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order] 
    ey = ey[order] 

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(mpeak+0.75, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade)/(1+z),np.std(tfade)/(1+z)))


def sn2011kl():
    # this is already shifted to rest-frame and corrected for extinction
    dat = ascii.read("../../data/2011kl.txt")
    xall = dat['col1']
    yall = dat['col2']
    yerrall = dat['col3']

    ind = np.argmin(yall)
    tpeak = xall[ind]
    mpeak = yall[ind]
    empeak = yerrall[ind]

    # The shape of the LC is weird so have to do this manually
    x = xall[xall<=tpeak].data
    y = yall[xall<=tpeak].data
    ey = yerrall[xall<=tpeak].data
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order]

    nsim = 1000
    trise= np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = tpeak-np.interp(mpeak+0.75, ysamples[ii], x)

    x = xall[xall>=tpeak].data
    y = yall[xall>=tpeak].data
    ey = yerrall[xall>=tpeak].data
    order = np.argsort(y)
    x = x[order]
    y = y[order] 
    ey = ey[order] 

    nsim = 1000
    tfade= np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(mpeak+0.75, ysamples[ii], x)-tpeak

    print("peak is %s+/-%s" %(mpeak,empeak)) # already abs mag
    print("rise time is %s +/- %s" %(np.mean(trise), np.std(trise)))
    print("fade time is %s +/- %s" %(np.mean(tfade), np.std(tfade)))


def koala():
    dat = ascii.read("../../data/koala_lc.txt")
    xall = dat['col1']
    order = np.argsort(xall)
    xall = xall[order]
    yall = dat['col2'][order]
    yerrall = dat['col3'][order]
    filt = dat['col4'][order]
    choose = filt == 'r'
    xall = xall[choose].data
    yall = yall[choose].data
    yerrall = yerrall[choose].data

    ind = np.argmax(yall)
    tpeak = xall[ind]
    mpeak = yall[ind]
    empeak = yerrall[ind]

    print("peak is %s+/-%s" %(mpeak,empeak))

    # The shape of the LC is weird so have to do this manually
    x = xall[3:ind+1]
    y = yall[3:ind+1]
    ey = yerrall[3:ind+1]
    # from averaging
    y[-1] = 103.4
    ey[-1] = 9.33
    y[0] = 2.95
    ey[0] = 8.33

    nsim = 1000
    trise= np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = tpeak-np.interp(mpeak/2, ysamples[ii], x)

    mean_rise = np.mean(trise)/1.2714
    err_rise = np.std(trise)/1.2714

    print("rise time is %s +/- %s" %(mean_rise,err_rise))

    x = xall[ind+1:14]
    y = yall[ind+1:14]
    ey = yerrall[ind+1:14]
    y[0] = 103.4
    ey[0] = 9.33
    y[3] = 60.05
    ey[3] = 8.07
    y[4] = 60.05
    ey[4] = 8.07
    x = x[4:6]
    y = y[4:6]
    ey = ey[4:6]

    nsim = 1000
    tfade= np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(mpeak/2, ysamples[ii], x)-tpeak

    mean_rise = np.mean(tfade)/1.2714
    err_rise = np.std(tfade)/1.2714

    print("fade time is %s +/- %s" %(mean_rise, err_rise))


def SNLS04D4ec():
    z = 0.593
    mw_ext = 0.047

    dat = pd.read_fwf("../../data/SNLS04D4ec.txt")
    jd = dat['JD']
    filt = dat['F']
    mag = dat['mag']
    emag = dat['emag']

    choose = filt == 'i'

    xall = jd[choose].values
    yall = mag[choose].values
    yerrall = emag[choose].values
    yall[-1] = 25.998
    yerrall[-1] = 0.1
    xall = xall.astype(float)
    yall = yall.astype(float)
    yerrall = yerrall.astype(float)

    ind = np.argmin(yall)
    tpeak = xall[ind]
    mpeak = yall[ind]
    empeak = yerrall[ind]
    mpeak = np.min(yall)

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall>=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order] 
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(mpeak+0.75, ysamples[ii], x)-tpeak

    rise = np.mean(trise)
    erise = np.std(trise)

    # Fading behavior
    x = xall[xall>=tpeak][0:4]
    y = yall[xall>=tpeak][0:4]
    ey = yerrall[xall>=tpeak][0:4]
    order = np.argsort(y)
    x = x[order]
    y = y[order] 
    ey = ey[order] 

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(mpeak+0.75, ysamples[ii], x)-tpeak

    fade = np.mean(tfade)
    efade = np.std(tfade)

    Mpeak, trise_corr, etrise_corr, tfade_corr, etfade_corr = correct(
        mpeak, rise, erise, fade, efade, z, mw_ext)

    print("mpeak is %s +/- %s" %(Mpeak, empeak))
    print("fade time is %s +/- %s" %(tfade_corr,etfade_corr))
    print("rise time is %s +/- %s" %(trise_corr,etrise_corr))


def SNLS05D2bk():
    z = 0.699
    mw_ext = 0.029

    dat = pd.read_fwf("../../data/SNLS05D2bk.txt")
    jd = dat['JD']
    filt = dat['F']
    mag = dat['mag']
    emag = dat['emag']

    choose = filt == 'i'

    xall = jd[choose].values[0:-2].astype(float)
    yall = mag[choose].values[0:-2].astype(float)
    yerrall = emag[choose].values[0:-2].astype(float)

    ind = np.argmin(yall)
    tpeak = xall[ind]
    mpeak = yall[ind]
    empeak = yerrall[ind]
    mpeak = np.min(yall)

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall>=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order] 
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(mpeak+0.75, ysamples[ii], x)-tpeak

    rise = np.mean(trise)
    erise = np.std(trise)

    # Fading behavior
    x = xall[xall>=tpeak][0:6]
    y = yall[xall>=tpeak][0:6]
    ey = yerrall[xall>=tpeak][0:6]
    order = np.argsort(y)
    x = x[order]
    y = y[order] 
    ey = ey[order] 

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(mpeak+0.75, ysamples[ii], x)-tpeak

    fade = np.mean(tfade)
    efade = np.std(tfade)

    Mpeak, trise_corr, etrise_corr, tfade_corr, etfade_corr = correct(
        mpeak, rise, erise, fade, efade, z, mw_ext)

    print("mpeak is %s +/- %s" %(Mpeak, empeak))
    print("fade time is %s +/- %s" %(tfade_corr,etfade_corr))
    print("rise time is %s +/- %s" %(trise_corr,etrise_corr))


def SNLS06D1hc():
    z = 0.555
    mw_ext = 0.048

    dat = pd.read_fwf("../../data/SNLS06D1hc.txt")
    jd = dat['JD']
    filt = dat['F']
    mag = dat['mag']
    emag = dat['emag']

    choose = filt == 'i'

    xall = jd[choose].values[0:-1].astype(float)
    yall = mag[choose].values[0:-1].astype(float)
    yerrall = emag[choose].values[0:-1].astype(float)

    ind = np.argmin(yall)
    tpeak = xall[ind]
    mpeak = yall[ind]
    empeak = yerrall[ind]
    mpeak = np.min(yall)

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall>=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order] 
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(mpeak+0.75, ysamples[ii], x)-tpeak

    rise = np.mean(trise)
    erise = np.std(trise)

    # Fading behavior
    x = xall[xall>=tpeak][0:5]
    y = yall[xall>=tpeak][0:5]
    ey = yerrall[xall>=tpeak][0:5]
    order = np.argsort(y)
    x = x[order]
    y = y[order] 
    ey = ey[order] 

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(mpeak+0.75, ysamples[ii], x)-tpeak

    fade = np.mean(tfade)
    efade = np.std(tfade)

    Mpeak, trise_corr, etrise_corr, tfade_corr, etfade_corr = correct(
        mpeak, rise, erise, fade, efade, z, mw_ext)

    print("mpeak is %s +/- %s" %(Mpeak, empeak))
    print("fade time is %s +/- %s" %(tfade_corr,etfade_corr))
    print("rise time is %s +/- %s" %(trise_corr,etrise_corr))


def DESX1eho():
    z = 0.76
    dat = pd.read_table("../../data/DES16X1eho_i.dat")
    xall = dat['# t'][1:-1].values.astype(float)
    yall= dat['M'][1:-1].values.astype(float)
    yerrall = dat['M_e'][1:-1].values.astype(float)
    # this is absolute magnitude in observer frame,
    # correctd for Galactic extinction

    ind = np.argmin(yall)
    tpeak = xall[ind]
    mpeak = yall[ind]
    empeak = yerrall[ind]

    # Fading behavior
    x = xall[xall>=tpeak][0:5]
    y = yall[xall>=tpeak][0:5]
    ey = yerrall[xall>=tpeak][0:5]
    order = np.argsort(y)
    x = x[order]
    y = y[order] 
    ey = ey[order] 

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(mpeak+0.75, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade)/(1+z),np.std(tfade)/(1+z)))


def ZTF18aalrxas():
    """ rest-frame g-band is just observer-frame (z<0.15) """
    mw_ext = 0.073
    z = 0.058832

    name = 'ZTF18aalrxas'
    ra = 237.298480
    dec = 32.288010
    jdobs = 2458229.8193
    zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[5,5])

    #dat = pd.read_table("ZTF18aalrxas_forced_phot.txt", sep=' ')
    #jd = dat['jd']
    #filt = dat['filt']
    #flux = dat['flux']
    #eflux = dat['eflux']

    choose = filt == 'g'

    xall = jd[choose]
    yall = flux[choose]
    yerrall = eflux[choose]

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]
    mpeak = mag[choose][ind]-0.073
    empeak = emag[choose][ind]
    Mpeak = mpeak-Planck15.distmod(z=z).value
    print("peak is %s +/- %s" %(Mpeak, empeak))

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise)/(1+z),np.std(trise)/(1+z)))

    # Fading behavior
    x = xall[xall>=tpeak]
    y = yall[xall>=tpeak]
    ey = yerrall[xall>=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order] 
    ey = ey[order] 

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade)/(1+z),np.std(tfade)/(1+z)))
     

def ZTF19aavbjfp():
    """ rest-frame g-band is just observer-frame (z<0.15) """
    mw_ext = 0.054
    z = 0.028

    name = 'ZTF19aavbjfp'
    ra = 186.753502 
    dec = 62.163866
    jdobs = 2458622.6656
    zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[7,70])

    #dat = pd.read_table("ZTF18aalrxas_forced_phot.txt", sep=' ')
    #jd = dat['jd']
    #filt = dat['filt']
    #flux = dat['flux']
    #eflux = dat['eflux']

    choose = filt == 'g'

    xall = jd[choose]
    yall = flux[choose]
    yerrall = eflux[choose]

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]
    mpeak = mag[choose][ind]-mw_ext
    empeak = emag[choose][ind]
    Mpeak = mpeak-Planck15.distmod(z=z).value
    print("peak is %s +/- %s" %(Mpeak, empeak))

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise)/(1+z),np.std(trise)/(1+z)))

    # Fading behavior
    x = xall[xall>=tpeak]
    y = yall[xall>=tpeak]
    ey = yerrall[xall>=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order] 
    ey = ey[order] 

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade)/(1+z),np.std(tfade)/(1+z)))
    

def ZTF19aapfmki():
    """ rest-frame g-band is just observer-frame (z<0.15) """
    mw_ext = 0.096
    z = 0.05469

    name = 'ZTF19aapfmki'
    ra = 211.431544
    dec = 9.515708
    jdobs = 2458586.7728
    zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[5,15])

    #dat = pd.read_table("ZTF18aalrxas_forced_phot.txt", sep=' ')
    #jd = dat['jd']
    #filt = dat['filt']
    #flux = dat['flux']
    #eflux = dat['eflux']

    choose = filt == 'g'

    xall = jd[choose]
    yall = flux[choose]
    yerrall = eflux[choose]

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]
    mpeak = mag[choose][ind]-mw_ext
    empeak = emag[choose][ind]
    Mpeak = mpeak-Planck15.distmod(z=z).value
    print("peak is %s +/- %s" %(Mpeak, empeak))

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise)/(1+z),np.std(trise)/(1+z)))

    # Fading behavior
    x = xall[xall>=tpeak]
    y = yall[xall>=tpeak]
    ey = yerrall[xall>=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order] 
    ey = ey[order] 

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade)/(1+z),np.std(tfade)/(1+z)))


def ZTF19abuvqgw():
    """ rest-frame g-band is just observer-frame (z<0.15) """
    mw_ext = 0.456
    z = 0.087

    name = 'ZTF19abuvqgw'
    ra = 297.526517
    dec = 66.082482
    jdobs = 2458730.7971

    # not enough 1DC data to make these measurements. need to use MSIP.
    # this is photometry from the marshal
    xall = np.array([2458726.7984, 2458727.8199, 2458727.8930, 2458730.7971, 
                     2458732.7954, 2458732.8111, 2458733.7948, 2458733.8285, 
                     2458734.7960, 2458734.8835, 2458735.9134])
    yall = np.array([19.94, 19.47, 19.36, 19.15, 19.37, 19.31, 19.44, 19.52,
                     19.64, 19.58, 19.93])
    yerrall = np.array([0.15, 0.12, 0.12, 0.11, 0.10, 0.12, 0.13, 0.12, 0.14, 
                        0.18, 0.16])

    ind = np.argmin(yall)
    tpeak = xall[ind]
    mpeak = yall[ind]
    empeak = yerrall[ind]
    Mpeak = mpeak-Planck15.distmod(z=z).value
    print("peak is %s +/- %s" %(Mpeak, empeak))

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(mpeak+0.75, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise)/(1+z),np.std(trise)/(1+z)))

    # Fading behavior
    x = xall[xall>=tpeak]
    y = yall[xall>=tpeak]
    ey = yerrall[xall>=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order] 
    ey = ey[order] 

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(mpeak+0.75, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade)/(1+z),np.std(tfade)/(1+z)))


def ZTF19abgbdcp():
    """ rest-frame g-band is just observer-frame (z<0.15) """
    mw_ext = 0.092
    z = 0.0318

    name = 'ZTF19abgbdcp'
    ra = 193.173301
    dec = 22.997260
    jdobs = 2458678.7610
    zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[10,20])

    choose = filt == 'g'

    xall = jd[choose]
    yall = flux[choose]
    yerrall = eflux[choose]
    # need to add one point for the preceding limit
    xall = np.append(xall, 2458674.7011)
    yall = np.append(yall, 0)
    yerrall = np.append(yerrall, 0.1)

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]
    mpeak = mag[choose][ind]-mw_ext
    empeak = emag[choose][ind]
    Mpeak = mpeak-Planck15.distmod(z=z).value
    print("peak is %s +/- %s" %(Mpeak, empeak))

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise)/(1+z),np.std(trise)/(1+z)))

    # Fading behavior
    # can't get forced photometry for public data here. 
    # just use the photometry from the marshal.
    # mpeak = 17.43, so half-max is 18.28
    x = [2458690.6979, 2458694.7128]
    y = [18.08, 18.34]
    ey = [0.10, 0.08]

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(18.28, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade)/(1+z),np.std(tfade)/(1+z)))


def ZTF18abskrix():
    """ rest-frame g-band is just observer-frame (z<0.15) """
    mw_ext = 0.585
    # this is a CV

    name = 'ZTF18abskrix'
    ra = 328.697754 
    dec = 34.722499
    jdobs = 2458364.7671
    zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[5,5])

    choose = filt == 'g'

    xall = jd[choose]
    yall = flux[choose]
    yerrall = eflux[choose]

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]
    mpeak = mag[choose][ind]-mw_ext
    empeak = emag[choose][ind]
    print("peak mag is %s +/- %s" %(mpeak,empeak))

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise),np.std(trise)))

    x = xall[xall>=tpeak]
    y = yall[xall>=tpeak]
    ey = yerrall[xall>=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade),np.std(tfade)))


def ZTF18absrffm():
    """ rest-frame g-band is just observer-frame (z<0.15) """
    # this is a CV

    name = 'ZTF18absrffm'
    ra = 260.235923 
    dec = 7.809376
    jdobs = 2458368.7287
    zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[5,5])

    choose = filt == 'g'
    xall = jd[choose]
    yall = flux[choose]
    yerrall = eflux[choose]

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]

    # peak mag
    mag[choose][ind]-0.545

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise),np.std(trise)))

    # Fading behavior. There aren't enough g-band points on the decline,
    # so let's use r-band instead.
    choose = filt == 'r'
    xall = jd[choose]
    yall = flux[choose]
    yerrall = eflux[choose]

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]

    x = xall[xall>=tpeak]
    y = yall[xall>=tpeak]
    ey = yerrall[xall>=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade),np.std(tfade)))


def ZTF18abyzkeq():
    """ rest-frame g-band is just observer-frame (z<0.15) """
    # this is a CV

    name = 'ZTF18abyzkeq'
    ra = 342.392602 
    dec = 37.931542
    jdobs = 2458387.7838
    zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[5,5])

    choose = filt == 'g'
    xall = jd[choose]
    yall = flux[choose]
    yerrall = eflux[choose]

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]

    # Peak mag
    mag[ind]-0.443

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise),np.std(trise)))

    # fading behavior
    x = xall[xall>=tpeak]
    y = yall[xall>=tpeak]
    ey = yerrall[xall>=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade),np.std(tfade)))


def ZTF18aasaiyp():
    """ rest-frame g-band is just observer-frame (z<0.15) """
    mw_ext = 0.053
    z = 0.104

    name = 'ZTF18aasaiyp'
    ra = 205.583613
    dec = 60.565342
    jdobs = 2458248.7201
    zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[10,20])

    choose = filt == 'g'

    xall = jd[choose]
    yall = flux[choose]
    yerrall = eflux[choose]

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]
    mpeak = mag[choose][ind]-mw_ext
    empeak = emag[choose][ind]
    Mpeak = mpeak-Planck15.distmod(z=z).value
    print("peak is %s +/- %s" %(Mpeak, empeak))

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise)/(1+z),np.std(trise)/(1+z)))

    # Fading behavior
    x = xall[xall>=tpeak]
    y = yall[xall>=tpeak]
    ey = yerrall[xall>=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade)/(1+z),np.std(tfade)/(1+z)))


def ZTF18ablxawt():
    """ rest-frame g-band is just observer-frame (z<0.15) """
    # this is a CV

    name = 'ZTF18ablxawt'
    ra = 270.155281 
    dec = 10.805404
    jdobs = 2458340.6722
    zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[10,10])

    choose = np.logical_and(filt == 'r', ~np.isnan(flux)) # no g-band data
    xall = jd[choose]
    yall = flux[choose]
    yerrall = eflux[choose]

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]

    # peak mag
    mag[choose][ind]-0.362

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order][2:]
    y = y[order][2:]
    ey = ey[order][2:]

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise),np.std(trise)))

    # Fading behavior. 
    x = xall[xall>=tpeak]
    y = yall[xall>=tpeak]
    ey = yerrall[xall>=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade),np.std(tfade)))


def ZTF18abuvqgo():
    """ rest-frame g-band is r-band (z>0.15) """
    mw_ext = 0.079
    z = 0.156

    name = 'ZTF18abuvqgo'
    ra = 23.953878 
    dec = 0.393257
    jdobs = 2458372.9142
    zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[10,20])

    choose = filt == 'r'

    xall = jd[choose]
    order = np.argsort(xall)
    xall = xall[order]
    yall = flux[choose][order]
    yerrall = eflux[choose][order]

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]
    mpeak = mag[choose][order][ind]-mw_ext
    empeak = emag[choose][order][ind]
    Mpeak = mpeak-Planck15.distmod(z=z).value
    print("peak is %s +/- %s" %(Mpeak, empeak))

    # Rising behavior
    x = xall[xall<=tpeak][11:13]
    y = yall[xall<=tpeak][11:13]
    ey = yerrall[xall<=tpeak][11:13]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise)/(1+z),np.std(trise)/(1+z)))

    # Fading behavior
    x = xall[xall>=tpeak][9:]
    y = yall[xall>=tpeak][9:]
    ey = yerrall[xall>=tpeak][9:]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order]

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade)/(1+z),np.std(tfade)/(1+z)))


def ZTF18abydmfv():
    """ rest-frame g-band is g-band (z<0.15) """
    mw_ext = 0.324
    z = 0.042

    name = 'ZTF18abydmfv'
    ra = 81.489304 
    dec = 2.082508
    jdobs = 2458384.9822
    zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[5,15])

    choose = filt == 'g'

    xall = jd[choose]
    order = np.argsort(xall)
    xall = xall[order]
    yall = flux[choose][order]
    yerrall = eflux[choose][order]

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]
    mpeak = mag[choose][order][ind]-mw_ext
    empeak = emag[choose][order][ind]
    Mpeak = mpeak-Planck15.distmod(z=z).value
    print("peak is %s +/- %s" %(Mpeak, empeak))

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise)/(1+z),np.std(trise)/(1+z)))

    # Fading behavior
    x = xall[xall>=tpeak][4:]
    y = yall[xall>=tpeak][4:]
    ey = yerrall[xall>=tpeak][4:]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order]

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade)/(1+z),np.std(tfade)/(1+z)))


def ZTF18acepuyx():
    """ rest-frame g-band is g-band (z<0.15) """
    mw_ext = 0.053
    z = 0.0711

    name = 'ZTF18acepuyx'
    ra = 160.035059 
    dec = 31.162211
    jdobs = 2458468.9449
    zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[5,20])

    choose = filt == 'g'

    xall = jd[choose]
    order = np.argsort(xall)
    xall = xall[order]
    yall = flux[choose][order]
    yerrall = eflux[choose][order]

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]
    mpeak = mag[choose][order][ind]-mw_ext
    empeak = emag[choose][order][ind]
    Mpeak = mpeak-Planck15.distmod(z=z).value
    print("peak is %s +/- %s" %(Mpeak, empeak))

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise)/(1+z),np.std(trise)/(1+z)))

    # Fading behavior
    x = xall[xall>=tpeak]
    y = yall[xall>=tpeak]
    ey = yerrall[xall>=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order]

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade)/(1+z),np.std(tfade)/(1+z)))


def ZTF19aatoboa():
    """ rest-frame g-band is g-band (z<0.15) """
    mw_ext = 0.044
    z = 0.0758

    name = 'ZTF19aatoboa'
    ra = 186.418990 
    dec = 44.746862
    jdobs = 2458609.7189
    zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[5,20])

    choose = filt == 'g'

    xall = jd[choose]
    order = np.argsort(xall)
    xall = xall[order]
    yall = flux[choose][order]
    yerrall = eflux[choose][order]

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]
    mpeak = mag[choose][order][ind]-mw_ext
    empeak = emag[choose][order][ind]
    Mpeak = mpeak-Planck15.distmod(z=z).value
    print("peak is %s +/- %s" %(Mpeak, empeak))

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise)/(1+z),np.std(trise)/(1+z)))

    # Fading behavior
    x = xall[xall>=tpeak][0:5]
    y = yall[xall>=tpeak][0:5]
    ey = yerrall[xall>=tpeak][0:5]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order]

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade)/(1+z),np.std(tfade)/(1+z)))


def ZTF19abyggpx():
    """ rest-frame g-band is g-band (z<0.15) """
    mw_ext = 0.048
    z = 0.0625

    name = 'ZTF19abgbbpx'
    ra = 192.165824 
    dec = 30.577645
    jdobs = 2458683.6760
    zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[10,20])

    choose = filt == 'g'

    xall = jd[choose]
    order = np.argsort(xall)
    xall = xall[order]
    yall = flux[choose][order]
    yerrall = eflux[choose][order]

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]
    mpeak = mag[choose][order][ind]-mw_ext
    empeak = emag[choose][order][ind]
    Mpeak = mpeak-Planck15.distmod(z=z).value
    print("peak is %s +/- %s" %(Mpeak, empeak))

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise)/(1+z),np.std(trise)/(1+z)))

    # Fading behavior
    x = xall[xall>=tpeak]
    y = yall[xall>=tpeak]
    ey = yerrall[xall>=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order]

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade)/(1+z),np.std(tfade)/(1+z)))


def ZTF19abiyyhd():
    """ rest-frame g-band is g-band (z<0.15) """
    mw_ext = 0.380
    z = 0.07

    name = 'ZTF19abiyyhd'
    ra = 357.534252 
    dec = 36.825260
    jdobs = 2458693.9109
    zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[5,10])

    choose = np.logical_and.reduce((filt == 'g', flux<140, flux>-20))

    xall = jd[choose]
    order = np.argsort(xall)
    xall = xall[order]
    yall = flux[choose][order]
    yerrall = eflux[choose][order]

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]
    mpeak = mag[choose][order][ind]-mw_ext
    empeak = emag[choose][order][ind]
    Mpeak = mpeak-Planck15.distmod(z=z).value
    print("peak is %s +/- %s" %(Mpeak, empeak))

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise)/(1+z),np.std(trise)/(1+z)))

    # Fading behavior
    x = xall[xall>=tpeak][:-3]
    y = yall[xall>=tpeak][:-3]
    ey = yerrall[xall>=tpeak][:-3]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order]

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade)/(1+z),np.std(tfade)/(1+z)))


def ZTF18acnbfel():
    """ rest-frame g-band is g-band (z<0.15) """
    mw_ext = 0.160
    z = 0.084

    name = 'ZTF18acnbfel'
    ra = 124.602467 
    dec = 64.413799
    jdobs = 2458442.9944
    zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[30,30])

    choose = np.logical_and(filt == 'g', flux<400)

    xall = jd[choose]
    order = np.argsort(xall)
    xall = xall[order]
    yall = flux[choose][order]
    yerrall = eflux[choose][order]

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]
    mpeak = mag[choose][order][ind]-mw_ext
    empeak = emag[choose][order][ind]
    Mpeak = mpeak-Planck15.distmod(z=z).value
    print("peak is %s +/- %s" %(Mpeak, empeak))

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise)/(1+z),np.std(trise)/(1+z)))

    # Fading behavior
    x = xall[xall>=tpeak][:-3]
    y = yall[xall>=tpeak][:-3]
    ey = yerrall[xall>=tpeak][:-3]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order]

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade)/(1+z),np.std(tfade)/(1+z)))


def ZTF19aaadfcp():
    """ rest-frame g-band is g-band (z<0.15) """
    mw_ext = 0.185
    z = 0.08

    name = 'ZTF19aaadfcp'
    ra = 120.019249 
    dec = 58.761185
    jdobs = 2458490.7184
    zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[10,10])

    choose = np.logical_and(filt == 'g', flux<250)

    xall = jd[choose]
    order = np.argsort(xall)
    xall = xall[order]
    yall = flux[choose][order]
    yerrall = eflux[choose][order]

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]
    mpeak = mag[choose][order][ind]-mw_ext
    empeak = emag[choose][order][ind]
    Mpeak = mpeak-Planck15.distmod(z=z).value
    print("peak is %s +/- %s" %(mpeak, empeak))

    # Rising behavior
    x = xall[xall<=tpeak][2:]
    y = yall[xall<=tpeak][2:]
    ey = yerrall[xall<=tpeak][2:]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise)/(1+z),np.std(trise)/(1+z)))

    # Fading behavior
    x = xall[xall>=tpeak]
    y = yall[xall>=tpeak]
    ey = yerrall[xall>=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order]

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade)/(1+z),np.std(tfade)/(1+z)))


def ZTF19aanvhyc():
    """ rest-frame g-band is g-band (z<0.15) """
    mw_ext = 0.343
    z = 0.065

    name = 'ZTF19aanvhyc'
    ra = 270.597724 
    dec = 19.058252
    jdobs = 2458577.9966
    zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[10,25])

    choose = np.logical_and(filt == 'g', flux<400)

    xall = jd[choose]
    order = np.argsort(xall)
    xall = xall[order]
    yall = flux[choose][order]
    yerrall = eflux[choose][order]

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]
    mpeak = mag[choose][order][ind]-mw_ext
    empeak = emag[choose][order][ind]
    Mpeak = mpeak-Planck15.distmod(z=z).value
    print("peak is %s +/- %s" %(mpeak, empeak))

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise)/(1+z),np.std(trise)/(1+z)))

    # Fading behavior
    x = xall[xall>=tpeak]
    y = yall[xall>=tpeak]
    ey = yerrall[xall>=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order]

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade)/(1+z),np.std(tfade)/(1+z)))


def ZTF19abqiacu():
    """ rest-frame g-band is g-band (z<0.15) """
    mw_ext = 0.163
    z = 0.1

    name = 'ZTF19abqiacu'
    ra = 269.132757 
    dec = 67.987523
    jdobs = 2458722.7978
    zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[20,25])

    choose = np.logical_and(filt == 'g', flux<400)

    xall = jd[choose]
    order = np.argsort(xall)
    xall = xall[order]
    yall = flux[choose][order]
    yerrall = eflux[choose][order]

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]
    mpeak = mag[choose][order][ind]-mw_ext
    empeak = emag[choose][order][ind]
    Mpeak = mpeak-Planck15.distmod(z=z).value
    print("peak is %s +/- %s" %(Mpeak, empeak))

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise)/(1+z),np.std(trise)/(1+z)))

    # Fading behavior
    x = xall[xall>=tpeak]
    y = yall[xall>=tpeak]
    ey = yerrall[xall>=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order]

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade)/(1+z),np.std(tfade)/(1+z)))


def ZTF18abxxeai():
    """ rest-frame g-band is just observer-frame (z<0.15) """
    # this is a CV

    name = 'ZTF18abxxeai'
    ra = 264.675448 
    dec = 15.201592
    jdobs = 2458381.6688
    zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[5,10])

    choose = np.logical_and(filt == 'g', ~np.isnan(flux))
    xall = jd[choose]
    yall = flux[choose]
    yerrall = eflux[choose]

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]

    # peak mag
    mag[choose][ind]-0.246

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise),np.std(trise)))

    # Fading behavior. 
    x = xall[xall>=tpeak]
    y = yall[xall>=tpeak]
    ey = yerrall[xall>=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade),np.std(tfade)))


def ZTF18abyjgaa():
    """ rest-frame g-band is just observer-frame (z<0.15) """
    # this is a CV

    name = 'ZTF18abyjgaa'
    ra = 93.617640 
    dec = 60.243844
    jdobs = 2458384.9898
    zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[5,10])

    choose = np.logical_and(filt == 'g', ~np.isnan(flux))
    xall = jd[choose]
    yall = flux[choose]
    yerrall = eflux[choose]

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]

    # Peak mag
    mag[choose][ind]-0.460

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise),np.std(trise)))

    # Fading behavior. 
    x = xall[xall>=tpeak]
    y = yall[xall>=tpeak]
    ey = yerrall[xall>=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade),np.std(tfade)))


def ZTF18acgnwpo():
    """ rest-frame g-band is just observer-frame (z<0.15) """
    # this is a CV

    name = 'ZTF18acgnwpo'
    ra = 104.305144 
    dec = 49.834241
    jdobs = 2458431.8840
    zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[5,10])

    choose = np.logical_and(filt == 'g', flux>-15)
    xall = jd[choose]
    yall = flux[choose]
    yerrall = eflux[choose]

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]

    # peak mag
    mag[choose][ind]-0.309

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise),np.std(trise)))

    # Fading behavior. 
    x = xall[xall>=tpeak]
    y = yall[xall>=tpeak]
    ey = yerrall[xall>=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade),np.std(tfade)))


def ZTF19aanqqzb():
    """ rest-frame g-band is just observer-frame (z<0.15) """
    # this is a CV

    name = 'ZTF19aanqqzb'
    ra = 67.603414 
    dec = 72.268338
    jdobs = 2458572.6601
    zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[5,5])

    choose = np.logical_and(filt == 'g', flux>-15)
    xall = jd[choose]
    yall = flux[choose]
    yerrall = eflux[choose]

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]

    # Peak mag
    mag[choose][ind]-0.824

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise),np.std(trise)))

    # Fading behavior. 
    x = xall[xall>=tpeak]
    y = yall[xall>=tpeak]
    ey = yerrall[xall>=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade),np.std(tfade)))


def ZTF19aaqfdvu():
    """ rest-frame g-band is just observer-frame (z<0.15) """
    # this is a CV

    name = 'ZTF19aaqfdvu'
    ra = 267.642963 
    dec = 9.682246
    jdobs = 2458587.9628
    zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[5,5])

    choose = np.logical_and.reduce((filt == 'g', ~np.isnan(flux), flux<200))
    xall = jd[choose]
    yall = flux[choose]
    yerrall = eflux[choose]

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]

    # Peak mag
    mag[choose][ind]-0.528

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise),np.std(trise)))

    # Fading behavior. 
    x = xall[xall>=tpeak]
    y = yall[xall>=tpeak]
    ey = yerrall[xall>=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade),np.std(tfade)))


def ZTF19aaxfqyx():
    """ rest-frame g-band is just observer-frame (z<0.15) """
    # this is a CV

    name = 'ZTF19aaxfqyx'
    ra = 277.900171 
    dec = 23.793139
    jdobs = 2458638.8728
    zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[5,5])

    choose = filt == 'g'
    xall = jd[choose]
    yall = flux[choose]
    yerrall = eflux[choose]

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]

    # peak mag
    mag[choose][ind]-0.344

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise),np.std(trise)))

    # Fading behavior. 
    x = xall[xall>=tpeak]
    y = yall[xall>=tpeak]
    ey = yerrall[xall>=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade),np.std(tfade)))


def ZTF19abfzfbs():
    """ rest-frame g-band is just observer-frame (z<0.15) """
    # this is a CV

    name = 'ZTF19abfzfbs'
    ra = 7.458464 
    dec = 48.103267
    jdobs = 2458681.9257
    zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[10,20])

    choose = filt == 'g'
    xall = jd[choose]
    yall = flux[choose]
    yerrall = eflux[choose]

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]

    # peak mag
    mag[choose][ind]-0.509

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise),np.std(trise)))

    # Fading behavior. 
    x = xall[xall>=tpeak]
    y = yall[xall>=tpeak]
    ey = yerrall[xall>=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade),np.std(tfade)))


def ZTF19abpwygn():
    """ rest-frame g-band is just observer-frame (z<0.15) """
    # this is a CV

    name = 'ZTF19abpwygn'
    ra = 64.363065
    dec = 73.895042
    jdobs = 2458711.9343
    zp,filt,jd,flux,eflux,mag,emag = get_forced_phot(name,ra,dec,jdobs,[5,5])

    choose = filt == 'g'
    xall = jd[choose]
    yall = flux[choose]
    yerrall = eflux[choose]

    ind = np.argmax(yall)
    tpeak = xall[ind]
    fpeak = yall[ind]
    efpeak = yerrall[ind]

    # peak mag
    mag[choose][ind]-0.477

    # Rising behavior
    x = xall[xall<=tpeak]
    y = yall[xall<=tpeak]
    ey = yerrall[xall<=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    trise = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        trise[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("rise time is %s +/- %s" %(np.mean(trise),np.std(trise)))

    # Fading behavior. 
    x = xall[xall>=tpeak]
    y = yall[xall>=tpeak]
    ey = yerrall[xall>=tpeak]
    order = np.argsort(y)
    x = x[order]
    y = y[order]
    ey = ey[order] 

    nsim = 1000
    tfade = np.zeros(nsim)

    ysamples = np.zeros((nsim,len(x)))
    for ii,val in enumerate(y):
        ysamples[:,ii] = np.random.normal(loc=val,scale=ey[ii],size=nsim)

    for ii in np.arange(nsim):
        tfade[ii] = np.interp(fpeak/2, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade),np.std(tfade)))
