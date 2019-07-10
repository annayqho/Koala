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


def get_query(search_start, search_end):
    query = """db['ZTF_alerts'].find(
        {'candidate.jd':
            {'$gt': %(startsearch)f, '$lt': %(endsearch)f},
         '$expr': {'$lt': [{'$subtract': ['$candidate.jd', '$candidate.jdstarthist']}, 100]},
         'candidate.rb': {'$gt': 0.3},
         'candidate.field': {'$eq': 554},
         'candidate.programid': {'$gt': 2},
         'candidate.ndethist': {'$gt': 2},
         'candidate.magpsf': {'$lt': 20},
         'candidate.ssdistnr': {'$lt': -1},
         'candidate.isdiffpos': {'$in': ['1', 't']}},
         {'objectId'}
        )""" %( {'startsearch': search_start, 'endsearch': search_end})
    return query


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
    # The very first observation of this field: 58345.44143477606
    # The last observation of this field before the break: 58386.4729638
    startsearch= 2458211
    endsearch= 2458215
    #endsearch = 2458585.0181908854

    # Koala dates
    startsearch = 2458372
    endsearch = 2458380

    s = logon()
    query_result = s.query({'query_type': 'general_search', 'query': get_query(
        startsearch, endsearch)})
    out = query_result['result_data']['query_result']
    with open('test.json', 'w') as f:
        json.dump(out, f)

    names = []
    dat = json.load(open("test.json"))
    for item in dat:
        names.append(item['objectId'])
    names = np.unique(np.array(names))
    print(len(names))

    keep = []
    for ii,name in enumerate(names):
        print(name)
        vals = get_details(name)
        if vals is not None:
            ssdistnr, distnr, sgscore, magnr, fwhm, nbad, chipsf = vals
            real = np.logical_and.reduce((
                fwhm>0.5, fwhm<8, nbad<5,chipsf<10))
            notrock = ssdistnr<0
            star = np.logical_or(np.logical_and.reduce(
                    (distnr<20, magnr<15, magnr>0, sgscore>0.49)),
                    np.logical_and(distnr<1,sgscore>0.8))
            has_host = distnr < 12
            if np.logical_and.reduce(
                    (real, notrock, ~star, has_host)):
                keep.append(name)

    np.savetxt("test_all.txt", keep, fmt='%s')
