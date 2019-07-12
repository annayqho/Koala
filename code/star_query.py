""" Next step of the query 
- get rid of asteroids
- get rid of stars
- require the LC to be nice
"""

import numpy as np
import glob
import json
import pandas
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF_fast_transient_search/code")
#from print_diagnostics import *
from kowalski_login import logon


def star_check(name):
    query = """db['ZTF_alerts'].find_one(
            {'objectId': {'$eq': '%s'}},
            {'candidate.distpsnr1', 'candidate.sgscore1', 'candidate.srmag1', 'candidate.sgmag1', 'candidate.simag1', 'candidate.szmag1'}
        )""" %name
    query_result = s.query(
            {'query_type': 'general_search', 'query': query})
    out = query_result['result_data']['query_result']
    dist = out['candidate']['distpsnr1']
    sg = out['candidate']['sgscore1']
    rmag = out['candidate']['srmag1']
    gmag = out['candidate']['sgmag1']
    imag = out['candidate']['simag1']
    zmag = out['candidate']['szmag1']

    pointunderneath = False
    if np.logical_and.reduce((dist>0, dist<2, sg>0.76)):
        pointunderneath = True
    if np.logical_and.reduce((sg==0.5, dist>0, dist<0.5)):
        if np.logical_or.reduce((
            np.logical_and(rmag>0,rmag<17),
            np.logical_and(gmag>0,gmag<17),
            np.logical_and(imag>0,imag<17),
            np.logical_and(zmag>0,zmag<17),
            )):
            pointunderneath = True

    brightstar = False
    if np.logical_and.reduce((rmag > 0, rmag < 12, sg > 0.49, dist < 20)):
        brightstar = True
    if np.logical_and.reduce((rmag > 0, rmag < 15, sg > 0.8, dist < 20)):
        brightstar = True
    if np.logical_and.reduce((gmag > 0, gmag < 12, sg > 0.49, dist < 20)):
        brightstar = True
    if np.logical_and.reduce((gmag > 0, gmag < 15, sg > 0.8, dist < 20)):
        brightstar = True
    if np.logical_and.reduce((imag > 0, imag < 12, sg > 0.49, dist < 20)):
        brightstar = True
    if np.logical_and.reduce((imag > 0, imag < 15, sg > 0.8, dist < 20)):
        brightstar = True
    if np.logical_and.reduce((zmag > 0, zmag < 12, sg > 0.49, dist < 20)):
        brightstar = True
    if np.logical_and.reduce((zmag > 0, zmag < 15, sg > 0.8, dist < 20)):
        brightstar = True
    
    return np.logical_or(pointunderneath, brightstar)


if __name__=="__main__":
    s = logon()

    # Read in all of the candidates, and only choose the unique ones 
    files = glob.glob("*filter1.txt")
    for inputf in files[5:]:
        field = inputf.split("_")[0]
        print(field)
        newfname = field+"_nostars.txt"
        if len(glob.glob(newfname)) == 0:
            print("removing stars...")
            cands = np.loadtxt(inputf,dtype=str)
            ucands = np.unique(cands)
            ncands = len(ucands)
            if ncands > 0:
                isstar = np.zeros(ucands.shape, dtype=bool)
                for ii,cand in enumerate(ucands):
                    isstar[ii] = star_check(cand)
                frac_star = np.round(sum(isstar)/ncands * 100, 2)
                print("removed %s percent of sources" %frac_star)
                np.savetxt(field+"_nostars.txt", ucands[~isstar], fmt='%s')
