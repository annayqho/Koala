""" Calculate rise and fade timescales for light curves """

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import pandas as pd


def iptf15ul():
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

    x = xall[xall<=tpeak].data
    y = yall[xall<=tpeak].data
    ey = yerrall[xall<=tpeak].data
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
        tfade[ii] = tpeak-np.interp(mpeak+0.75, ysamples[ii], x)

    print("rise time is %s +/- %s" %(np.mean(tfade),np.std(tfade)))

    x = xall[xall>=tpeak].data 
    y = yall[xall>=tpeak].data 
    ey = yerrall[xall>=tpeak].data 
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
        tfade[ii] = np.interp(mpeak+0.75, ysamples[ii], x)-tpeak

    print("fade time is %s +/- %s" %(np.mean(tfade),np.std(tfade)))


def at2018cow():
    dat = pd.read_fwf("../../data/at2018cow_photometry_table.dat")
    mjd = dat['MJD']
    filt = dat['Filt']
    mag = dat['ABMag']
    emag = dat['Emag']
    choose = filt == 'g'

    xall = mjd[choose][1:].values
    yall = mag[choose][1:].values.astype(float)
    yerrall = emag[choose][1:].values.astype(float)

    tpeak = xall[np.argmin(yall)]
    mpeak = np.min(yall)

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
    dat = ascii.read("../../data/dougie.txt")
    xall = dat['col1']
    yall = dat['col2']+0.031
    yerrall = np.sqrt(dat['col3']**2+0.1**2)

    ind = np.argmin(yall)
    tpeak = xall[ind]
    mpeak = yall[ind]
    empeak = yerrall[ind]

    print("peak is %s+/-%s" %(mpeak,empeak))

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

    print("rise time is %s +/- %s" %(np.mean(trise),np.std(trise)))

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

    print("fade time is %s +/- %s" %(np.mean(tfade),np.std(tfade)))


def iptf16asu():
    dat = pd.read_fwf("../../data/iptf16asu.txt")
    mjd = dat['MJD']
    filt = dat['filt']
    mag = dat['mag']
    emag = dat['emag']

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


def sn2011kl():
    dat = ascii.read("../../data/2011kl.txt")
    xall = dat['col1']
    yall = dat['col2']
    yerrall = dat['col3']

    ind = np.argmin(yall)
    tpeak = xall[ind]
    mpeak = yall[ind]
    empeak = yerrall[ind]

    print("peak is %s+/-%s" %(mpeak,empeak))

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

    print("rise time is %s +/- %s" %(np.mean(trise),np.std(trise)))

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

    print("fade time is %s +/- %s" %(np.mean(tfade),np.std(tfade)))


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

    print("rise time is %s +/- %s" %(np.mean(trise),np.std(trise)))

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

    print("fade time is %s +/- %s" %(np.mean(tfade),np.std(tfade)))


def SNLS04D4ec():
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

    print("rise time is %s +/- %s" %(np.mean(trise),np.std(trise)))

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

    print("fade time is %s +/- %s" %(np.mean(tfade),np.std(tfade)))


def SNLS05D2bk():
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

    print("rise time is %s +/- %s" %(np.mean(trise),np.std(trise)))

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

    print("fade time is %s +/- %s" %(np.mean(tfade),np.std(tfade)))


def SNLS06D1hc():
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

    print("rise time is %s +/- %s" %(np.mean(trise),np.std(trise)))

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

    print("fade time is %s +/- %s" %(np.mean(tfade),np.std(tfade)))
