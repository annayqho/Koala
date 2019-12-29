""" Scan the partnership and Caltech 1DC streams.
This was previously run as a Jupyter notebook,
but converted to a Python script for easier use.


To do items:
    - Make the cutouts the same size
    - Print out the distance from the host (distnr, distpsnr, etc)
    - Make cutouts live links to webpages
    - Print the dRB score
    - Say whether the source has already been saved to the marshal


"""

import numpy as np
import warnings
import os
import requests
from PIL import Image
from astropy.cosmology import Planck15
from bson.json_util import loads, dumps
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
import subprocess
import pandas as pd
import glob
#import progressbar
from penquins import Kowalski
import psycopg2
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import ICRS
from astropy.io import fits
from ztfquery.utils import stamps
import gzip
import io
import json
from helpers import *


def logon():
    """ Log onto Kowalski """
    username = 'ah'
    password = 'TetraodonInsists'
    s = Kowalski(
        protocol='https', host='kowalski.caltech.edu', port=443,
        verbose=False, username=username, password=password)
    return s


def count_obs(obs, program):
    """ Count the number of observations taken
    in a given field, in each filter """
    choose = obs['qcomment'] == program
    fields = np.unique(obs['field'][choose])
    counts = np.zeros((len(fields), 2))

    for ii,fieldid in enumerate(fields):
        filts = [1,2]
        counts[ii,:] = np.array(
            [np.sum(np.logical_and(obs['field'][choose]==fieldid, obs['fil'][choose]==f)) for f in filts])
        
    return counts


def initial_search(s, obs):
    """ Run an initial broad query """

    # Set search window
    obst = obs['UT_START'].values
    start = Time(obst[0], format='isot').jd
    end = Time(obst[-1], format='isot').jd

    q = {"query_type": "find",
         "query": {
             "catalog": "ZTF_alerts",
             "filter": {
                     'candidate.jd': {'$gt': start, '$lt': end},
                     'candidate.magpsf': {'$lt': 20},
                     'candidate.isdiffpos': {'$in': ['1', 't']},
                     'candidate.programid': {'$gt': 0},
                     'candidate.ssdistnr': {'$lt': -1},
                     'candidate.drb': {'$gt': 0.5},
             },
             "projection": {
                     "objectId": 1,
                     "candidate.distpsnr1": 1,
                     "candidate.sgscore1": 1,
                     "candidate.srmag1": 1,
                     "candidate.sgmag1": 1,
                     "candidate.simag1": 1,
                     "candidate.szmag1": 1,
             }
         }
         }

    query_result = s.query(query=q)
    out = query_result['result_data']['query_result']

    names_all = np.array([val['objectId'] for val in out])
    names = np.unique(names_all)
    print("There are %s unique cands from this initial filter." %len(names))

    dist = np.array([val['candidate']['distpsnr1'] for val in out]) 
    sg = np.array([val['candidate']['sgscore1'] for val in out]) 
    rmag = np.array([val['candidate']['srmag1'] for val in out])
    gmag = np.array([val['candidate']['sgmag1'] for val in out])
    imag = np.array([val['candidate']['simag1'] for val in out])
    zmag = np.array([val['candidate']['szmag1'] for val in out])

    pointunderneath = np.zeros(len(names_all), dtype=bool)
    crit = np.logical_and.reduce((dist>0, dist<2, sg>0.76))
    pointunderneath[crit] = True

    crit = np.logical_and(
            np.logical_and.reduce((sg==0.5, dist>0, dist<0.5)),
            np.logical_or.reduce((
                np.logical_and(rmag>0,rmag<17),
                np.logical_and(gmag>0,gmag<17),
                np.logical_and(imag>0,imag<17),
                np.logical_and(zmag>0,zmag<17),
            )))
    pointunderneath[crit] = True

    crit = np.logical_or.reduce((
        np.logical_and.reduce((rmag > 0, rmag < 12, sg > 0.49, dist < 20)),
        np.logical_and.reduce((rmag > 0, rmag < 15, sg > 0.49, dist < 5)),
        np.logical_and.reduce((rmag > 0, rmag < 15, sg > 0.8, dist < 20)),
        np.logical_and.reduce((gmag > 0, gmag < 12, sg > 0.49, dist < 20)),
        np.logical_and.reduce((gmag > 0, gmag < 15, sg > 0.8, dist < 20)),
        np.logical_and.reduce((imag > 0, imag < 12, sg > 0.49, dist < 20)),
        np.logical_and.reduce((imag > 0, imag < 11.5, sg > 0.49, dist < 20)),
        np.logical_and.reduce((imag > 0, imag < 14.5, sg > 0.49, dist < 5)),
        np.logical_and.reduce((imag > 0, imag < 15, sg > 0.8, dist < 20)),
        np.logical_and.reduce((zmag > 0, zmag < 11.5, sg > 0.49, dist < 10)),
        np.logical_and.reduce((zmag > 0, zmag < 14.0, sg > 0.49, dist < 2.5)),
        np.logical_and.reduce((zmag > 0, zmag < 15, sg > 0.8, dist < 20))))
    pointunderneath[crit] = True

    names_no_star = np.unique(names_all[pointunderneath==False])
    print("There are %s that are not stars." %len(names_no_star))
    return np.array(names_no_star)


def get_dets(s, name):
    q = {"query_type": "find",
         "query": {
             "catalog": "ZTF_alerts",
             "filter": {
                     'objectId': {'$eq': name},
                     'candidate.isdiffpos': {'$in': ['1', 't']},
             },
             "projection": {
                     "_id": 0,
                     "candidate.jd": 1,
                     "candidate.magpsf": 1,
                     "candidate.sigmapsf": 1,
                     "candidate.fid": 1,
                     "candidate.programid": 1,
             }
         }  
         }  
    query_result = s.query(query=q)
    out = query_result['result_data']['query_result']
    return out


def get_prv_dets(s, name):
    q = {"query_type": "find",
         "query": {
             "catalog": "ZTF_alerts_aux",
             "filter": {
                     '_id': {'$eq': name},
             },
             "projection": {
                     "_id": 0,
                     "prv_candidates": 1,
             }
         }  
         }  
    query_result = s.query(query=q)
    out = query_result['result_data']['query_result'][0]['prv_candidates']
    return out


def plot_band(ax, dt,mag,emag,mask,col):
    indpeak = np.argmin(mag[mask])
    tpeak = dt[mask][indpeak]
    ax.errorbar(dt[mask], mag[mask], emag[mask], c=col, fmt='o')
    ax.axhline(y=mag[mask][indpeak], c=col, lw=0.5)
    ax.axvline(x=tpeak, c=col, lw=0.5)
    ax.axhline(y=mag[mask][indpeak]+0.75, c=col, lw=0.5)


def plot_lc(ax,name,dt,mag,emag,filt):
    if sum(filt==1) > 0:
        plot_band(ax, dt, mag, emag, filt==1, 'green')
    if sum(filt==2) > 0:
        plot_band(ax, dt, mag, emag, filt==2, 'red')
    if sum(filt==3) > 0:
        plot_band(ax, dt, mag, emag, filt==3, 'grey')

    ax.set_xlabel("Days", fontsize=14)
    ax.set_ylabel("App Mag", fontsize=14)
    ax.invert_yaxis()
    ax.tick_params(axis='both', labelsize=14)


def plotlims_band(ax, dt,mag,mask,col):
    ax.scatter(
            dt[mask], mag[mask], edgecolor=col, facecolor='white', marker='v')


def plot_lims(ax, limdt,limmag,limfilt):
    plotlims_band(ax, limdt, limmag, limfilt==1, 'green')
    plotlims_band(ax, limdt, limmag, limfilt==2, 'red')


def get_nondet(s, name):
    q = {"query_type": "find",
         "query": {
             "catalog": "ZTF_alerts",
             "filter": {
                     'objectId': {'$eq': name}
             },
             "projection": {
                     "_id": 0,
                     "prv_candidates": 1,
             }
         }  
         }  
    query_result = s.query(query=q)
    out = query_result['result_data']['query_result']
    
    limjds = []
    limmags = []
    limfilts = []
    for res in out:
        if res['prv_candidates'] is not None:
            for det in res['prv_candidates']:
                if det['magpsf'] is None:
                    limjds.append(det['jd'])
                    limmags.append(det['diffmaglim'])
                    limfilts.append(det['fid'])
    
    limjds = np.array(limjds)
    limmags = np.array(limmags)
    limfilts = np.array(limfilts)    
    
    order = np.argsort(np.array(limjds))
    limjds = limjds[order]
    limmags = limmags[order]
    limfilts = limfilts[order]
    
    return limjds,limmags,limfilts


def calc_trise(dt,mag,mask):
    if len(dt[mask]) > 0:
        indpeak = np.argmin(mag[mask])
        tpeak = dt[mask][indpeak]
        mpeak = mag[mask][indpeak]
        dt_fit = dt[mask][0:indpeak+1]
        mag_fit = mag[mask][0:indpeak+1]
        order = np.argsort(mag_fit)
        dt_fit = dt_fit[order]
        mag_fit = mag_fit[order]
        tval = np.interp(mpeak+0.75, mag_fit, dt_fit, left=-999, right=-999)
        return tpeak-tval
    else:
        return 1000
    
    
def calc_tfade(dt,mag,mask):
    if len(dt[mask]) > 0:
        indpeak = np.argmin(mag[mask])
        tpeak = dt[mask][indpeak]
        mpeak = mag[mask][indpeak]
        dt_fit = dt[indpeak:]
        mag_fit = mag[indpeak:]
        order = np.argsort(mag_fit)
        dt_fit = dt_fit[order]
        mag_fit = mag_fit[order]
        tval = np.interp(mpeak+0.75, mag_fit, dt_fit, left=-999, right=-999)
        return tpeak-tval
    else:
        return 1000
    
    
def calc_dmdt(dt,mag,emag):
    dm = np.array(mag[1:]-mag[0:-1])
    dt = np.array(dt[1:]-dt[0:-1])
    dmdt = np.abs(dm/dt)
    sig = np.array(np.abs(dm)/np.sqrt(emag[1:]**2 + emag[0:-1]**2))
    # only pay attention to dm/dt values that are significant
    choose = sig > 1
    if sum(choose) > 0:
        return np.max(dmdt[choose])
    return 0


def make_triplet(alert, normalize: bool = False):
    """
        Feed in alert packet
    """
    cutout_dict = dict()

    for cutout in ('science', 'template', 'difference'):
        cutout_data = loads(
                dumps([alert[f'cutout{cutout.capitalize()}']['stampData']]))[0]

        # unzip
        with gzip.open(io.BytesIO(cutout_data), 'rb') as f:
            with fits.open(io.BytesIO(f.read())) as hdu:
                data = hdu[0].data
                # replace nans with zeros
                cutout_dict[cutout] = np.nan_to_num(data)
                # normalize
                if normalize:
                    cutout_dict[cutout] /= np.linalg.norm(cutout_dict[cutout])

        # pad to 63x63 if smaller
        shape = cutout_dict[cutout].shape
        if shape != (63, 63):
            cutout_dict[cutout] = np.pad(cutout_dict[cutout], [(0, 63 - shape[0]), (0, 63 - shape[1])],
                                         mode='constant', constant_values=1e-9)

    triplet = np.zeros((63, 63, 3))
    triplet[:, :, 0] = cutout_dict['science']
    triplet[:, :, 1] = cutout_dict['template']
    triplet[:, :, 2] = cutout_dict['difference']

    return triplet


def plot_triplet(tr):
    """ From Dima's Kowalski tutorial """
    fig,axarr = plt.subplots(1,3,figsize=(8, 2), dpi=120)
    titles = ['Science', 'Reference', 'Difference']
    u_scale_factor = [40, 40, 10]
    l_scale_factor = [30, 40, 1]
    for ii,ax in enumerate(axarr):
        ax.axis('off')
        data = tr[:,:,ii]
        dat = data.flatten()
        sig = np.median(np.abs(dat-np.median(dat)))
        median = np.median(data)
        ax.imshow(
            data, origin='upper', cmap=plt.cm.bone, 
            vmin=median-l_scale_factor[ii]*sig, 
            vmax=median+u_scale_factor[ii]*sig)
        #norm=LogNorm())
        ax.title.set_text(titles[ii])
    fig.subplots_adjust(wspace=0)
    return fig


def plot_source_summary(s, name):
    """ Plot one png file per source """

    #################################################
    # Only run any of this if the LC duration > 0.02d.
    #################################################

    # Get LC
    jd,mag,emag,fid,program,limjds,limmags,limfilts = get_lc(s, name)
    if jd[-1]-jd[0] > 0.02:
        # Set up the plot
        fig = plt.figure(figsize=(8,10))
        gs = gridspec.GridSpec(5, 5, figure=fig, hspace=0.4, wspace=0.4)

        #################################################
        # Top row: basic information about the candidate
        #################################################

        # Position
        query = """db['ZTF_alerts'].find_one(
                {'objectId': {'$eq': '%s'}},
                {'_id': 0})""" %name
        query_result = s.query(
                {'query_type': 'general_search', 'query': query})
        out = query_result['result_data']['query_result']
        ra = out['candidate']['ra']
        dec = out['candidate']['dec']
        rastr = out['coordinates']['radec_str'][0]
        decstr = out['coordinates']['radec_str'][1]
        
        ax = plt.subplot(gs[0,0:2])
        ax.text(0,0.9,"Candidate: %s" %name, transform=ax.transAxes)
        pos = SkyCoord(ra, dec, unit='deg')
        ax.text(0,0.7,"RA, Dec: %s,%s" %(np.round(ra,5),np.round(dec,5)), 
                transform=ax.transAxes)
        ax.text(0, 0.5, "%s,%s" %(rastr,decstr))

        # In Galactic coordinates as well (to see latitude)
        c = ICRS(ra*u.degree, dec*u.degree)
        ax.text(0, 0.3, "Galactic (l,b): %s,%s" %(
            np.round(
                pos.galactic.l.deg,5), np.round(pos.galactic.b.deg, 5)),
                transform=ax.transAxes)

        # Formatting
        ax.set_axis_off()

        # Now show cutouts 
        tr = make_triplet(out)
        ax = plt.subplot(gs[0,2])
        ax.axis('off')
        ax.imshow(tr[:, :, 0], origin='upper', cmap=plt.cm.bone, norm=LogNorm())
        ax.title.set_text('Science')
        ax = plt.subplot(gs[0,3])
        ax.axis('off')
        ax.imshow(tr[:, :, 1], origin='upper', cmap=plt.cm.bone, norm=LogNorm())
        ax.title.set_text('Reference')
        ax = plt.subplot(gs[0,4])
        ax.axis('off')
        ax.imshow(tr[:, :, 2], origin='upper', cmap=plt.cm.bone)
        ax.title.set_text('Difference')

        #################################################
        # Second row: added information from other surveys
        #################################################

        # PS1
        ax = plt.subplot(gs[1,3])
        img = stamps.get_ps_stamp(ra, dec, size=240, color=["y","g","i"])
        ax.imshow(np.asarray(img))
        ax.set_title("PS1 (y/g/i)")
        ax.axis('off')

        # LegacySurvey
        ax = plt.subplot(gs[1,4])
        url = "http://legacysurvey.org/viewer/jpeg-cutout?ra=%s&dec=%s&layer=dr8&pixscale=0.27&bands=grz" %(ra,dec)
        r = requests.get(url)
        ax.imshow(Image.open(io.BytesIO(r.content)))
        ax.set_title("LegSurv DR8")
        ax.axis('off')

        # Cross-match to get photo-z
        ax = plt.subplot(gs[1,0])
        conn = psycopg2.connect(
                host="private.caltech.edu", user="ztfuser", 
                password="86R03DpSnTvtKsM&VVo9", database="ztfimages", 
                port="6666")
        cur = conn.cursor()
        cur.execute('SELECT z_phot_median,"TYPE" FROM dr8_north \
                     WHERE q3c_radial_query("RA", "DEC", %s, %s, 0.0002777 * 3);' %(ra,dec))
        out_north = cur.fetchall()
        cur.execute('SELECT z_phot_median,"TYPE" FROM dr8_south \
                     WHERE q3c_radial_query("RA", "DEC", %s, %s, 0.0002777 * 3);' %(ra,dec))
        out_south = cur.fetchall()
        out = [val for val in out_north]
        for val in out_south:
            out.append(val)
        nsrc = len(out)
        ax.text(
                0,0.9,"LegacySurvey: %s sources in 3arcsec" %nsrc, 
                transform=ax.transAxes)
        if nsrc > 0:
            photoz = out[0][0]
            stype = out[0][1]
            ax.text(
                    0,0.7,"Closest source: photoz=%s, type=%s" %(
                        np.round(photoz,2), stype), 
                    transform=ax.transAxes)
            dm = Planck15.distmod(z=photoz).value
            ax.text(
                    0,0.5,"distmod=%s" %(np.round(dm,2)),
                    transform=ax.transAxes)
        # Age: time since first detection
        age = np.round(Time.now().jd-jd[0],1)
        ax.text(0,0.1,"Age (since first det): %s days" %(age))
        ax.axis('off')

        #################################################
        #################################################
        #################################################
        # Second row: full light curve
        #################################################

        # Decide on a t0 for the plot
        t0 = jd[np.argmin(mag)]

        # Plot the FULL LC
        ax = plt.subplot(gs[2,:]) # second row, all columns
        plot_lc(ax,name,jd-t0,mag,emag,fid)
        plot_lims(ax,limjds-t0,limmags,limfilts)  

        # Calculate some basic numbers for each filter
        for fidx in [1, 2]:
            choose = fid == fidx
            if sum(choose) > 1:
                jd_fit = jd[choose]
                dt_fit = jd_fit-t0
                mag_fit = mag[choose]
                emag_fit = emag[choose]
                ndet = len(jd_fit)

                choose_lim = limfilts == fidx
                limdt_fit = limjds[choose_lim]-t0
                limmags_fit = limmags[choose_lim]
                limfilts_fit = limfilts[choose_lim]

                ax = plt.subplot(gs[fidx+2,0:2])
                plot_lc(
                        ax,name,dt_fit,mag_fit,emag_fit,
                        np.array([fidx]*ndet))
                plot_lims(ax,limdt_fit,limmags_fit,limfilts_fit)  
                ax.set_xlim(-20,20)

                # Calculate the duration
                src_duration = dt_fit[-1]-dt_fit[0]
                ax.text(
                        0.9, 0.15, 
                        "Duration [d]: %s" %np.round(src_duration, 2),
                        transform=ax.transAxes, horizontalalignment='right')

                # Calculate the max dmag/dt of the detections
                dmdt_det = calc_dmdt(dt_fit, mag_fit, emag_fit)

                # If it exists,
                # calculate the dmag/dt of the last non-detection
                # to the first detection
                choose_lim = limfilts == fidx
                before_first_det = limjds[choose_lim] < jd_fit[0]
                if sum(before_first_det) > 0:
                    lastlim_jd = limjds[choose_lim][before_first_det][-1] 
                    lastlim_dt = lastlim_jd - t0
                    lastlim_mag = limmags[choose_lim][before_first_det][-1] 
                    dmdt_lim = calc_dmdt(
                            np.array([lastlim_dt, dt_fit[0]]), 
                            np.array([lastlim_mag, mag_fit[0]]),
                            np.array([0, emag_fit[0]]))

                    # Keep whichever is bigger,
                    # assuming that the limiting mag is actually
                    # fainter than the detection
                    if np.logical_and(dmdt_lim>dmdt_det,
                            lastlim_mag > mag_fit[0]):
                        dmdt = dmdt_lim
                    else:
                        dmdt = dmdt_det
                else:
                    dmdt = dmdt_det
                ax.text(
                        0.90, 0.05, "dm/dt [mag/d]: %s" %np.round(dmdt, 2),
                        transform=ax.transAxes, horizontalalignment='right')

                # Check whether there is a peak
                indpeak = np.argmin(mag_fit)
                mpeak = mag_fit[indpeak]
                empeak = emag_fit[indpeak]
                dmrise = mag_fit[0]-mpeak
                dmfade = mag_fit[-1]-mpeak
                rise_sig = dmrise/np.sqrt(empeak**2+emag_fit[0]**2)
                fade_sig = dmfade/np.sqrt(empeak**2+emag_fit[-1]**2)

                has_peak = np.logical_and.reduce((
                        mag_fit[0] > mpeak, mag_fit[-1] > mpeak,
                        rise_sig > 1, fade_sig > 1))

                # If the source has peaked, check the rise time
                ax = plt.subplot(gs[fidx+2,3:])
                tpeak = dt_fit[indpeak]
                dt_rise = dt_fit[0:indpeak+1]
                mag_rise = mag_fit[0:indpeak+1]
                emag_rise = emag_fit[0:indpeak+1]
                ndet = len(dt_rise)
                plot_lc(
                        ax,name,dt_rise,mag_rise,emag_rise,
                        np.array([fidx]*ndet))
                ax.set_xlim(min(dt_rise), max(dt_rise))
                if has_peak:
                    order = np.argsort(mag_rise)
                    dt_rise_sorted = dt_rise[order]
                    mag_rise_sorted = mag_rise[order]
                    tval = np.interp(
                            mpeak+0.75, mag_rise_sorted, dt_rise_sorted, 
                            left=-999, right=-999)
                    rise_time = tpeak-tval
                    if rise_time > 998:
                        dur_rise = dt_rise[-1]-dt_rise[0]
                        ax.text(
                                0.90, 0.05, 
                                "$t_{1/2}$ [d] > %s" %np.round(dur_rise,2),
                                transform=ax.transAxes, 
                                horizontalalignment='right')
                    else:
                        ax.text(
                                0.90, 0.05, 
                                "$t_{1/2}$ [d]: %s" %np.round(rise_time,2),
                                transform=ax.transAxes, 
                                horizontalalignment='right')
                else:
                    ax.text(
                            0.90, 0.05, 
                            "No significant peak", 
                            transform=ax.transAxes, 
                            horizontalalignment='right')


        # Formatting
        #gs.tight_layout(fig)

        # Save figure
        #plt.show()
        plt.savefig("%s.png" %(name))
        plt.close()


def print_source_line(outputf, s, cur, name):
    """ Print one line to html file per source """
    # Get LC
    isalert,jd,mag,emag,fid,program,limjds,limmags,limfilts,limprog = get_lc(s, name)

    # Get position
    query = """db['ZTF_alerts'].find_one(
            {'objectId': {'$eq': '%s'}},
            {'_id': 0})""" %name
    query_result = s.query(
            {'query_type': 'general_search', 'query': query})
    out = query_result['result_data']['query_result']
    ra = out['candidate']['ra']
    dec = out['candidate']['dec']
    rastr = out['coordinates']['radec_str'][0]
    decstr = out['coordinates']['radec_str'][1]

    # In Galactic coordinates as well (to see latitude)
    pos = SkyCoord(ra, dec, unit='deg')
    gal_l = pos.galactic.l.deg
    gal_b = pos.galactic.b.deg
    
    # Only do this if duration is between 1 and 100 days
    isdet = mag < 99
    dur = max(jd[isdet])-min(jd[isdet])
    if np.logical_and(dur>1,dur<100):
        # Print name
        outputf.write(name)
        outputf.write('&nbsp;&nbsp;&nbsp;&nbsp;')
        if sum(program==2) > 0:
            outputf.write('[<a href="http://skipper.caltech.edu:8080/cgi-bin/growth/ingest_avro_id.cgi?programidx=1&name=%s">Save to Marshal (Partnership)</a>]' %name)
            outputf.write('&nbsp;&nbsp;&nbsp;&nbsp;')
        if sum(program==3) > 0:
            outputf.write('[<a href="http://skipper.caltech.edu:8080/cgi-bin/growth/ingest_avro_id.cgi?programidx=11&name=%s">Save to Marshal (Caltech)</a>]' %name)
            outputf.write('&nbsp;&nbsp;&nbsp;&nbsp;')
        outputf.write('[<a href="http://skipper.caltech.edu:8080/cgi-bin/growth/view_source.cgi?name=%s">Marshal</a>]' %name)
        outputf.write('&nbsp;&nbsp;&nbsp;&nbsp;')
        outputf.write('[<a href="https://lasair.roe.ac.uk/object/%s">Lasair</a>]' %name)
        outputf.write('</br>')

        outputf.write("RA, Dec: %s, %s" %(np.round(ra,5),np.round(dec,5)))
        outputf.write('</br>')
        outputf.write("%s, %s" %(rastr,decstr))
        outputf.write('</br>')

        #c = ICRS(ra*u.degree, dec*u.degree)
        outputf.write(
            "Galactic (l, b): %s, %s" %(np.round(gal_l,5), np.round(gal_b,5)))
        outputf.write('</br>')
        outputf.write('</br>')

        # Now show cutouts of science image
        tr = make_triplet(out)
        fig = plot_triplet(tr) 
        plt.savefig("%s_triplet.png" %name, bbox_inches='tight', pad_inches=0)
        plt.close()
        outputf.write('<img src="%s_triplet.png">'%(name))

        # PS1
        img = stamps.get_ps_stamp(ra, dec, size=240, color=["y","g","i"])
        plt.figure(figsize=(2,2))
        plt.imshow(np.asarray(img))
        plt.title("PS1 (y/g/i)")
        plt.axis('off')
        plt.savefig("%s_ps1.png" %name, bbox_inches='tight', pad_inches=0)
        plt.close()
        outputf.write('<img src="%s_ps1.png">'%(name))

        #LegacySurvey
        url = "http://legacysurvey.org/viewer/jpeg-cutout?ra=%s&dec=%s&layer=dr8&pixscale=0.27&bands=grz" %(ra,dec)
        r = requests.get(url)
        plt.figure(figsize=(2,2))
        plt.imshow(Image.open(io.BytesIO(r.content)))
        plt.title("LegSurv DR8")
        plt.axis('off')
        plt.savefig("%s_ls.png" %name)
        plt.close()
        outputf.write('<img src="%s_ls.png">'%(name))
        outputf.write('</br>')

        # Cross-match to get photo-z
        if dec > 32:
            cur.execute('SELECT z_phot_median,z_phot_l68,z_phot_u68,"TYPE" FROM dr8_north \
                         WHERE q3c_radial_query("RA", "DEC", %s, %s, 0.0002777);' %(ra,dec))
            out = cur.fetchall()
        else:
            cur.execute('SELECT z_phot_median,z_phot_l68,z_phot_u68,"TYPE" FROM dr8_south \
                         WHERE q3c_radial_query("RA", "DEC", %s, %s, 0.0002777);' %(ra,dec))
            out = cur.fetchall()
        nsrc = len(out)
        outputf.write("LegacySurvey: %s sources in 1arcsec" %nsrc)
        outputf.write('</br>')
        if nsrc > 0:
            photoz = out[0][0]
            lphotoz = out[0][1]
            uphotoz = out[0][2]
            stype = out[0][3]
            print(lphotoz,uphotoz)
            outputf.write(
                    "Closest source: photoz=%s (68 percentile bounds %s,%s), type=%s" %(
                    np.round(photoz,2), np.round(lphotoz,2), 
                    np.round(uphotoz,2), stype))
            outputf.write('</br>')
            dm = Planck15.distmod(z=photoz).value
            outputf.write("distmod=%s" %(np.round(dm,2)))
            outputf.write('</br>')

        # Calculate the duration
        outputf.write("Duration [d]: %s" %np.round(dur, 2))
        outputf.write('</br>')

        # Decide on a t0 for the plot
        t0 = jd[np.argmin(mag)]

        # Plot the FULL LC
        fig,ax = plt.subplots(1,1,figsize=(8,2))
        plot_lc(ax,name,jd-t0,mag,emag,fid)
        plot_lims(ax,limjds-t0,limmags,limfilts)  
        plt.tight_layout()
        plt.savefig("%s_lc.png" %name)
        plt.close()
        outputf.write('<img src="%s_lc.png">'%(name))

        # Plot a zoomed-in version
        fig,ax = plt.subplots(1,1,figsize=(3,2))
        plot_lc(ax,name,jd-t0,mag,emag,fid)
        plot_lims(ax,limjds-t0,limmags,limfilts)  
        plt.xlim(-5,5)
        choose = np.logical_and(jd-t0>-5,jd-t0<5)
        plt.ylim(max(mag[choose]+0.5), min(mag[choose])-0.5)
        plt.tight_layout()
        plt.savefig("%s_zoom_lc.png" %name)
        plt.close()
        outputf.write('<img src="%s_zoom_lc.png">'%(name))
        outputf.write('</br>')

        # Calculate some basic numbers for each filter
        for fidx in [1, 2]:
            filtname = {1: 'g', 2: 'r'}
            choose = fid == fidx
            if sum(choose) > 1:
                jd_fit = jd[choose]
                dt_fit = jd_fit-t0
                mag_fit = mag[choose]
                emag_fit = emag[choose]
                ndet = len(jd_fit)

                choose_lim = limfilts == fidx
                limdt_fit = limjds[choose_lim]-t0
                limmags_fit = limmags[choose_lim]
                limfilts_fit = limfilts[choose_lim]


                # Calculate the max dmag/dt of the detections
                dmdt_det = calc_dmdt(dt_fit, mag_fit, emag_fit)

                # If it exists,
                # calculate the dmag/dt of the last non-detection
                # to the first detection
                choose_lim = limfilts == fidx
                before_first_det = limjds[choose_lim] < jd_fit[0]
                if sum(before_first_det) > 0:
                    lastlim_jd = limjds[choose_lim][before_first_det][-1] 
                    lastlim_dt = lastlim_jd - t0
                    lastlim_mag = limmags[choose_lim][before_first_det][-1] 
                    dmdt_lim = calc_dmdt(
                            np.array([lastlim_dt, dt_fit[0]]), 
                            np.array([lastlim_mag, mag_fit[0]]),
                            np.array([0, emag_fit[0]]))

                    # Keep whichever is bigger,
                    # assuming that the limiting mag is actually
                    # fainter than the detection
                    if np.logical_and(dmdt_lim>dmdt_det,
                            lastlim_mag > mag_fit[0]):
                        dmdt = dmdt_lim
                    else:
                        dmdt = dmdt_det
                else:
                    dmdt = dmdt_det

                # Check whether there is a peak
                indpeak = np.argmin(mag_fit)
                mpeak = mag_fit[indpeak]
                empeak = emag_fit[indpeak]
                dmrise = mag_fit[0]-mpeak
                dmfade = mag_fit[-1]-mpeak
                rise_sig = dmrise/np.sqrt(empeak**2+emag_fit[0]**2)
                fade_sig = dmfade/np.sqrt(empeak**2+emag_fit[-1]**2)

                has_peak = np.logical_and.reduce((
                        mag_fit[0] > mpeak, mag_fit[-1] > mpeak,
                        rise_sig > 1, fade_sig > 1))

                # If the source has peaked, check the rise time
                tpeak = dt_fit[indpeak]
                dt_rise = dt_fit[0:indpeak+1]
                mag_rise = mag_fit[0:indpeak+1]
                emag_rise = emag_fit[0:indpeak+1]
                ndet = len(dt_rise)
                fig,ax = plt.subplots(1,1,figsize=(3,2))
                plot_lc(
                        ax,name,dt_rise,mag_rise,emag_rise,
                        np.array([fidx]*ndet))
                ax.set_xlim(min(dt_rise), max(dt_rise))
                plt.tight_layout()
                plt.savefig("%s_%s_rise.png" %(name,fidx))
                plt.close()
                outputf.write('<img src="%s_%s_rise.png">'%(name,fidx))
                outputf.write('</br>')
                if has_peak:
                    order = np.argsort(mag_rise)
                    dt_rise_sorted = dt_rise[order]
                    mag_rise_sorted = mag_rise[order]
                    tval = np.interp(
                            mpeak+0.75, mag_rise_sorted, dt_rise_sorted, 
                            left=-999, right=-999)
                    rise_time = tpeak-tval
                    if rise_time > 998:
                        dur_rise = dt_rise[-1]-dt_rise[0]
                        outputf.write(
                                "%s-band: t1/2 [d] > %s" %(
                                    filtname[fidx], np.round(dur_rise,2)))
                        outputf.write('</br>')
                    else:
                        outputf.write(
                                "%s-band: t1/2 [d]: %s" %(
                                    filtname[fidx], np.round(rise_time,2)))
                        outputf.write('</br>')
                else:
                    outputf.write(
                            "%s-band: no significant peak" %filtname[fidx])
                    outputf.write('</br>')
                outputf.write("dm/dt [mag/d]: %s" %np.round(dmdt, 2))
                outputf.write("</br>")


        # # Formatting
        # #gs.tight_layout(fig)

        # # Save figure
        # #plt.show()
        # plt.savefig("%s.png" %(name))
        # plt.close()
        outputf.write('</br>')
        outputf.write('<hr>')


def make_marshal(names):
    s = logon()

    # Connect to LegacySurvey photo-z table
    conn = psycopg2.connect(
            host="private.caltech.edu", user="ztfuser", 
            password="86R03DpSnTvtKsM&VVo9", database="ztfimages", 
            port="6666")
    cur = conn.cursor()

    # Set up the mini marshal
    outputf = open("1DC_candidates.html", "w")
    outputf.write('<!doctype html>')
    outputf.write('<html>')
    outputf.write('<head>')
    outputf.write('<title>SOURCE MARSHAL</title>')
    outputf.write('</head>')
    outputf.write('<body>')

    outputf.write('<font size="5" color="green">Sources from the 1DC survey</font>')
    outputf.write('</br>')
    outputf.write('</br>')
    outputf.write('</br>')
    outputf.write('<hr>')
    outputf.write('</br>') 

    for ii,name in enumerate(names):
        print("%s/%s" %(ii, len(names)))
        # For debugging
        print(name)
        print_source_line(outputf, s, cur, name)

    outputf.close()


if __name__=="__main__":
    files = glob.glob("output/*duration.txt")
    names = []
    for f in files:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            dat = np.loadtxt(f, dtype=str, ndmin=1)
            if len(dat) > 0:
                [names.append(val) for val in dat]
    names = np.array(names)
    names = np.unique(names)
    make_marshal(names)
    #make_marshal(['ZTF18abvkwla', 'ZTF18aalrxas'])
