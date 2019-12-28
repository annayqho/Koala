""" Create a mini marshal for a given observation date 
and list of candidate names """

import matplotlib.pyplot as plt
import numpy as np
import io
import glob
from penquins import Kowalski
import psycopg2
from bson.json_util import loads, dumps
import gzip
from ztfquery.utils import stamps
from ztfquery import query as zq
from helpers import *
from astropy.io import fits
import requests
from PIL import Image
from astropy.cosmology import Planck15
import warnings


def logon():
    """ Log onto Kowalski """
    username = 'ah'
    password = 'TetraodonInsists'
    s = Kowalski(
        protocol='https', host='kowalski.caltech.edu', port=443,
        verbose=False, username=username, password=password)
    return s


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


def print_source_line(outputf, s, cur, name):
    """ Print one line to html file per source """
    conn = psycopg2.connect(
            host="private.caltech.edu", user="ztfuser",
            password="86R03DpSnTvtKsM&VVo9", database="ztfimages",
            port="6666")
    cur = conn.cursor()

    # Get LC
    isalert,jd,mag,emag,fid,program,limjds,limmags,limfilts,limprogram = get_lc(s, name)
    choose = program==3
    isalert = isalert[choose]
    jd = jd[choose]
    mag = mag[choose]
    emag = emag[choose]
    fid = fid[choose]
    choose = limprogram==3
    limjds = limjds[choose]
    limmags = limmags[choose]
    limfilts = limfilts[choose]
    
    # Query for the alert
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
    field = out['candidate']['field']
    drb = out['classifications']['braai']
    rb = out['candidate']['rb']
    distnr = out['candidate']['distnr']
    magnr = out['candidate']['magnr']
    pos = Position(ra,dec,rastr,decstr)
    
    # Print name
    outputf.write(
            '<a href="http://www.astro.caltech.edu/~ah/archival_search_partnership/candidates.html#%s"><b><font size="3">%s</font></b></a>'%(name,name))
    outputf.write( 
            '<div id="%s"></div>'%(name))
    outputf.write('&nbsp;&nbsp;&nbsp;&nbsp;')
    if np.logical_or(sum(program==2) > 0, sum(program==1) > 0): 
        outputf.write('[<a href="http://skipper.caltech.edu:8080/cgi-bin/growth/ingest_avro_id.cgi?programidx=1&name=%s">Save to Marshal (Partnership)</a>]' %name)
        outputf.write('&nbsp;&nbsp;&nbsp;&nbsp;')
    if sum(program==3) > 0:
        outputf.write('[<a href="http://skipper.caltech.edu:8080/cgi-bin/growth/ingest_avro_id.cgi?programidx=11&name=%s">Save to Marshal (Caltech)</a>]' %name)
        outputf.write('&nbsp;&nbsp;&nbsp;&nbsp;')
    outputf.write('[<a href="http://skipper.caltech.edu:8080/cgi-bin/growth/view_source.cgi?name=%s">Marshal</a>]' %name)      
    outputf.write('&nbsp;&nbsp;&nbsp;&nbsp;')
    outputf.write('[<a href="https://lasair.roe.ac.uk/object/%s">Lasair</a>]' %name)
    outputf.write('</br>')

    outputf.write("RA, Dec: %s, %s" %(np.round(pos.ra,5),np.round(pos.dec,5)))
    outputf.write('</br>')
    outputf.write("%s, %s" %(pos.rastr,pos.decstr))
    outputf.write('</br>')

    outputf.write(
        "Galactic (l, b): %s, %s" %(np.round(pos.l,5), np.round(pos.b,5)))
    outputf.write('</br>')
    outputf.write('</br>')
    
    # Now show cutouts of science image
    fname = "%s_triplet.png" %name
    if len(glob.glob(fname))==0:
        tr = make_triplet(out)
        fig = plot_triplet(tr)
        plt.savefig(
                fname, bbox_inches='tight', pad_inches=0)
        plt.close()
    outputf.write('<img src="%s_triplet.png">'%(name))
    
    # PS1
    fname = "%s_ps1.png" %name
    if len(glob.glob(fname))==0:
        img = stamps.get_ps_stamp(
                pos.ra, pos.dec, size=240, color=["y","g","i"])
        plt.figure(figsize=(2,2))
        plt.imshow(np.asarray(img))
        plt.title("PS1 (y/g/i)")
        plt.axis('off')
        plt.savefig(
                "%s_ps1.png" %name, bbox_inches='tight',
                pad_inches=0)
        plt.close()
    outputf.write('<img src="%s_ps1.png">'%(name))

    #LegacySurvey
    fname = "%s_ls.png" %name
    if len(glob.glob(fname))==0:
        url = "http://legacysurvey.org/viewer/jpeg-cutout?ra=%s&dec=%s&layer=dr8&pixscale=0.27&bands=grz" %(ra,dec)
        r = requests.get(url)
        plt.figure(figsize=(2,2))
        plt.imshow(Image.open(io.BytesIO(r.content)))
        plt.title("LegSurv DR8")
        plt.axis('off')
        plt.savefig(fname)
        plt.close()
    outputf.write('<img src="%s_ls.png">'%(name))
    outputf.write('</br>')

    # Cross-match to get photo-z
    cur.execute('SELECT z_phot_median,z_phot_l68,z_phot_u68,"TYPE" FROM dr8_north \
                 WHERE q3c_radial_query("RA", "DEC", %s, %s, 0.0002777*1.5);' %(ra,dec))
    out_north = cur.fetchall()
    cur.execute('SELECT z_phot_median,z_phot_l68,z_phot_u68,"TYPE" FROM dr8_south \
                 WHERE q3c_radial_query("RA", "DEC", %s, %s, 0.0002777*1.5);' %(ra,dec))
    out_south = cur.fetchall()
    out = [val for val in out_north]
    for val in out_south: 
        out.append(val)
    nsrc = len(out)
    outputf.write("LegacySurvey: %s sources in 1.5arcsec" %nsrc)
    outputf.write('</br>')
    if nsrc > 0:
        photoz = out[0][0]
        lphotoz = out[0][1]
        uphotoz = out[0][2]
        stype = out[0][3]
        outputf.write(
                "Closest source: photoz=%s (68 percentile bounds %s,%s), type=%s" %(
                np.round(photoz,2), np.round(lphotoz,2), 
                np.round(uphotoz,2), stype))
        outputf.write('</br>')
        dm = Planck15.distmod(z=photoz).value
        outputf.write("distmod=%s" %(np.round(dm,2)))
        outputf.write('</br>')
        
    # dRB score
    outputf.write("drb score: %s" %(np.round(drb,2)))
    outputf.write('</br>')
    outputf.write("rb scor:: %s" %(np.round(rb,2)))
    outputf.write('</br>')

    # Age: time since first detection
    age = np.round(Time.now().jd-jd[isalert][0],1)
    outputf.write("Age (since first alert): %s days" %(age))
    outputf.write('</br>')
    
    # Calculate the duration
    src_duration = jd[isalert][-1]-jd[isalert][0]
    outputf.write("Duration [d]: %s" %np.round(src_duration, 2))
    outputf.write('</br>')
    
    # Decide on a t0 for the plot
    t0 = jd[np.argmin(mag)]
    
    # Figure out tnow
    tnow = Time.now().jd

    # Plot the FULL LC
    fname = "%s_lc.png" %name
    if len(glob.glob(fname))==0:
        fig,ax = plt.subplots(1,1,figsize=(4,2))
        plot_lc(ax,name,jd-t0,mag,emag,fid)
        plot_lims(ax,limjds-t0,limmags,limfilts)
        plt.axvline(x=tnow-t0, ls='--', c='k')
        plt.tight_layout()
        plt.savefig(fname)
        plt.close()
    outputf.write('<img src=%s>'%(fname))
    
    # Plot a zoomed-in version
    fname = "%s_zoom_lc.png"%(name)
    if len(glob.glob(fname))==0:
        xmin = -2.5
        xmax = 2.5
        fig,ax = plt.subplots(1,1,figsize=(3,2))
        choose = np.logical_and(jd-t0>xmin,jd-t0<xmax)
        plot_lc(ax,name,jd[choose]-t0,mag[choose],emag[choose],fid[choose])
        ymin = min(mag[choose])-0.5
        choose = np.logical_and(limjds-t0>xmin,limjds-t0<xmax)
        plot_lims(ax,limjds[choose]-t0,limmags[choose],limfilts[choose])

        # Get limiting magnitudes from ztfquery for +/- 3 nights
        zquery = zq.ZTFQuery()

        lim_window_min = jd[0]-3
        lim_window_max = jd[0]-0.5
        zquery.load_metadata(
            radec=[pos.ra,pos.dec], size=0.0001, 
            sql_query="obsjd>%s and obsjd<%s" %(lim_window_min,lim_window_max))
        lim_jd = zquery.metatable['obsjd']
        lim_filt = zquery.metatable['fid']
        lim_maglim = zquery.metatable['maglimit']
        plot_lims(ax,lim_jd-t0,lim_maglim,lim_filt)
        opt_a = 0
        if len(lim_maglim) > 0:
            plot_lims(ax,lim_jd-t0,lim_maglim,lim_filt)
            opt_a = max(lim_maglim)

        lim_window_min = jd[-1]+0.5
        lim_window_max = jd[-1]+3
        zquery.load_metadata(
            radec=[pos.ra,pos.dec], size=0.0001, 
            sql_query="obsjd>%s and obsjd<%s" %(lim_window_min,lim_window_max))
        lim_jd = zquery.metatable['obsjd']
        lim_filt = zquery.metatable['fid']
        lim_maglim = zquery.metatable['maglimit']
        opt_b = 0
        if len(lim_maglim) > 0:
            plot_lims(ax,lim_jd-t0,lim_maglim,lim_filt)
            opt_b = max(lim_maglim)

        ymax = np.max([opt_a,opt_b,max(mag)])

        plt.xlim(xmin,xmax)
        plt.ylim(ymin,ymax)
        plt.gca().invert_yaxis()
        plt.axvline(x=tnow-t0, ls='--', c='k')
        plt.tight_layout()
        plt.savefig("%s_zoom_lc.png" %name)

        plt.xlim(min(jd-t0)-0.1,max(jd-t0)+0.1)
        plt.savefig("%s_zoom_zoom_lc.png" %name)
        plt.close()
    outputf.write('<img src="%s_zoom_lc.png">'%(name))
    outputf.write('<img src="%s_zoom_zoom_lc.png">'%(name))
    outputf.write('</br>')
    outputf.write("LC from forced photometry:")
    outputf.write('</br>')
    outputf.write('<img src="forced_phot/%s_forced_phot.png">'%(name))
    outputf.write('</br>')
    outputf.write('</br>')
    outputf.write('<hr>')


def make_marshal(label,toplot):
    # Set up the mini marshal
    s = logon()
    conn = psycopg2.connect(
            host="private.caltech.edu", user="ztfuser",
            password="86R03DpSnTvtKsM&VVo9", database="ztfimages",
            port="6666")
    cur = conn.cursor()

    outputf = open("candidates_%s.html" %(label), "w")
    outputf.write('<!doctype html>')
    outputf.write('<html>')
    outputf.write('<head>')
    outputf.write('<title>SOURCE MARSHAL</title>')
    outputf.write('</head>')
    outputf.write('<body>')
    
    outputf.write('</br>')
    outputf.write('</br>')
    outputf.write('<hr>')
    outputf.write('</br>')
    
    for ii,name in enumerate(toplot):
        print("%s/%s" %(ii, len(toplot)))
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
    make_marshal("1DC", names)
