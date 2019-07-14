""" Print diagnostic plots for the list of sources that survived scanning """

import matplotlib.pyplot as plt
import requests
import numpy as np
import os
import pandas
from astropy.io import fits,ascii
from ztfquery import query
from ztfquery.utils import stamps
from penquins import Kowalski
import glob


names = np.loadtxt("/Users/annaho/Dropbox/Projects/Research/Koala/code/toscan_byeye.txt", dtype=str)

# Log onto Kowalski and ZTFquery
username = 'ah'
password = 'TetraodonInsists'
s = Kowalski(
    protocol='https', host='kowalski.caltech.edu', port=443,
    verbose=False, username=username, password=password)
zquery = query.ZTFQuery()

# Get the position
def get_pos(name):
    query = """db['ZTF_alerts'].find_one(
            {'objectId': {'$eq': '%s'}},
            {'candidate.ra', 'candidate.dec'}
        )""" %name
    query_result = s.query(
            {'query_type': 'general_search', 'query': query})
    out = query_result['result_data']['query_result']
    ra = out['candidate']['ra']
    dec = out['candidate']['dec']
    return ra, dec


# Get the detections
def get_dets(name):
    q = {"query_type": "find",
         "query": {
             "catalog": "ZTF_alerts",
             "filter": {
                     'objectId': {'$eq': name},
                     'candidate.rb': {'$gt': 0.3},
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
    try:
        out = query_result['result_data']['query_result']
        return out
    except:
        return None


def get_lc(name):
    out = get_dets(name)
    jd = []
    dt = []
    mag = []
    emag = []
    filt = []
    pid = []

    if out is not None:
        for det in out:
            cand = det['candidate']

            current_jd = cand['jd']
            current_mag = cand['magpsf']
            current_emag = cand['sigmapsf']
            current_filter = cand['fid']
            current_pid = cand['programid']

            jd.append(current_jd)
            mag.append(current_mag)
            emag.append(current_emag)
            filt.append(current_filter)
            pid.append(current_pid)

        jd = np.array(jd)
        mag = np.array(mag)
        emag = np.array(emag)
        filt = np.array(filt)
        pid = np.array(pid)

        # Sort in order of jd
        order = np.argsort(jd)
        jd = jd[order]
        dt = jd-jd[0]
        mag = mag[order]
        emag = emag[order]
        filt = filt[order] 
        pid = pid[order]
    return jd,dt,mag,emag,filt,pid


def plot_lc(name, ra, dec):
    jd,dt,mag,emag,filt,pid = get_lc(name)
    fname = "%s_lc.png" %name
    plt.figure()

    gband = filt == 1
    plt.errorbar(dt[gband], mag[gband], emag[gband], c='green', fmt='o', label="g-band")
    rband = filt == 2
    plt.errorbar(dt[rband], mag[rband], emag[rband], c='red', fmt='s', label="r-band")
    plt.legend(fontsize=14)
    plt.xlabel("Days since first detection", fontsize=14)
    plt.ylabel("App Mag", fontsize=14)
    plt.gca().invert_yaxis()
    plt.tick_params(axis='both', labelsize=14)
    plt.xlim(dt[0]-10, dt[-1]+10) 
    plt.title("%s, RA=%s, Dec=%s" %(name,ra,dec), fontsize=14)

    # print("searching for non detections around %s, %s" %(ra,dec))
    # zquery.load_metadata(
    #         radec=[ra,dec],size=0.0005,
    #         sql_query="(obsjd BETWEEN %s AND %s) OR (obsjd BETWEEN %s AND %s)" %(
    #             jd[0]-10,jd[0]-0.001,jd[-1]+0.001,jd[-1]+10))
    # obsjd = zquery.metatable['obsjd'].values
    # obsfilt = zquery.metatable['filtercode'].values
    # obslim = zquery.metatable['maglimit'].values
    # obsdt = obsjd-jd[0]

    # gband = obsfilt == 'zg'
    # plt.scatter(obsdt[gband], obslim[gband], c='green', marker='v')
    # rband = obsfilt == 'zr'
    # plt.scatter(obsdt[rband], obslim[rband], c='red', marker='v')

    plt.xlim(-5, 5)
    plt.tight_layout()
    plt.savefig("%s_lc.png" %name)
    plt.close()


def plot_cutout(name):
    fname = "%s_cutout.png" %name
    if len(glob.glob(fname))==0:
        query = """db['ZTF_alerts'].find_one(
                {'objectId': {'$eq': '%s'}},
                {'cutoutDifference'}
            )""" %name
        query_result = s.query(
                {'query_type': 'general_search', 'query': query})
        src = query_result['result_data']['query_result']
        cutout = src['cutoutDifference']
        imgname = cutout['fileName']
        imgdata = cutout['stampData']
        with open(imgname, 'wb') as f:
            f.write(imgdata)
        im = fits.open(imgname)[0].data
        plt.figure()
        plt.imshow(im)
        plt.title("P48 Difference Image")
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(fname, pad_inches=0, bbox_inches='tight')
        plt.close()  


def plot_ps1(name,ra,dec):
    fname = "%s_ps1.png" %name
    if len(glob.glob(fname))==0:
        img = stamps.get_ps_stamp(ra, dec, size=240, color=["y","g","i"])
        plt.figure()
        plt.imshow(np.asarray(img))
        plt.title("Pan-STARRS (y/g/i)")
        plt.axis('off')
        plt.tight_layout()
        
        plt.savefig(fname, bbox_inches='tight')
        plt.close()


def plot_ls(name,ra,dec):
    fname = '%s_ls.jpg' %name
    if len(glob.glob(fname)) == 0:
        url = "http://legacysurvey.org/viewer/jpeg-cutout?ra=%s&dec=%s&layer=decals-dr7&pixscale=0.27&bands=grz" %(ra,dec)
        r = requests.get(url)
        with open(fname, 'wb') as f:
            f.write(r.content)   


def run():
    names = np.loadtxt("toscan_byeye.txt", dtype=str)
    ind = np.where(names=='ZTF18aahjcrq')[0][0]
    for ii,name in enumerate(names[ind:]):
        print(ii,name)
        ra, dec = get_pos(name)
        print("plotting light curves")
        plot_lc(name, ra, dec)
        print("plotting cutouts")
        plot_cutout(name)
        print("plotting ps1 cutout")
        plot_ps1(name,ra,dec)
        print("plotting legacysurvey image")
        plot_ls(name,ra,dec)

if __name__=="__main__":
    run()
