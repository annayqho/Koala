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


def ztf_transients():
    # ZTF18aalrxas
    # rest-frame g-band is just observer-frame (z<0.15)
    ext_g = 0.073
    

if __name__=="__main__":
    at2018cow()
