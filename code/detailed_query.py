""" Query for fast-evolving transients

Criteria:
    - more than 2 detections at a given position with RB score > 0.3
    - brighter than mag of 20
"""

import numpy as np
import json
import pandas
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF_fast_transient_search/code")
#from print_diagnostics import *
from kowalski_login import logon


def get_details(name):
    query = """db['ZTF_alerts'].find_one(
            {'objectId': {'$eq': '%s'}},
            {'candidate.jdstarthist', 'candidate.jd', 'candidate.ssdistnr', 'candidate.distpsnr1', 'candidate.sgscore1', 'candidate.srmag1', 'candidate.fwhm', 'candidate.nbad', 'candidate.chipsf'}
        )""" %name
    query_result = s.query(
            {'query_type': 'general_search', 'query': query})
    out = query_result['result_data']['query_result']
    try:
        ssdistnr = out['candidate']['ssdistnr']
        distnr = out['candidate']['distpsnr1']
        sgscore = out['candidate']['sgscore1']
        magnr = out['candidate']['srmag1']
        fwhm = out['candidate']['fwhm']
        nbad = out['candidate']['nbad']
        chipsf = out['candidate']['chipsf']
    except:
        print("missing crucial fields")
        return None
    return ssdistnr, distnr, sgscore, magnr, fwhm, nbad, chipsf


if __name__=="__main__":
    allfields = np.loadtxt("fieldlist.txt")

    for field in allfields:
        print("running field %s" %field)
        s = logon()

        q = {"query_type": "find",
             "query": {
                 "catalog": "ZTF_alerts",
                 "filter": {
                         'candidate.jd': {'$gt': 2458211, '$lt': 2458585},
                         'candidate.field': {'$eq': field},
                         'candidate.rb': {'$gt': 0.3},
                         'candidate.ndethist': {'$gt': 2},
                         'candidate.magpsf': {'$lt': 20},
                         'candidate.isdiffpos': {'$in': ['1', 't']},
                         'candidate.programid': {'$gt': 2},
                 },
                 "projection": {
                         "_id": 0,
                         "objectId": 1
                 }
             }
             }

        query_result = s.query(query=q)
        out = query_result['result_data']['query_result']

        names = []
        for item in out:
            names.append(item['objectId'])
        names = np.unique(np.array(names))
        print("%s unique cands" %len(names))

        np.savetxt("field%s_filter1.txt" %field, names, fmt='%s')
