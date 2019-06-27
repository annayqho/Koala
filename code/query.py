""" Query for fast-evolving transients

Criteria:
    - more than 10 detections at a given position with RB score > 0.5
    - rising or fading at a rate >= 1 mag/day
    - jdendhist - jdstarthst < 1 year
"""

import numpy as np
import pandas
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/ZTF_fast_transient_search/code")
#from print_diagnostics import *
from kowalski_login import logon


def get_query(search_start, search_end):
    query = """db['ZTF_alerts'].find_one(
        {'candidate.jd':
            {'$gt': %(startsearch)f, '$lt': %(endsearch)f},
         'candidate.field': {'$eq': 554},
         'candidate.rb': {'$gt': 0.3},
         'candidate.ndethist': {'$gt': 2},
         'candidate.magpsf': {'$lt': 20},
         'candidate.isdiffpos': {'$in': ['1', 't']}},
         {'objectId'}
        )""" %( {'startsearch': search_start, 'endsearch': search_end})
    return query


def get_details(name):
    query = """db['ZTF_alerts'].find_one(
            {'objectId': {'$eq': '%s'}},
            {'candidate.ssdistnr', 'candidate.distpsnr1', 'candidate.sgscore1'}
        )""" %name
    query_result = s.query(
            {'query_type': 'general_search', 'query': query})
    out = query_result['result_data']['query_result']
    ssdistnr = out['candidate']['ssdistnr']
    distnr = out['candidate']['distpsnr1']
    sgscore = out['candidate']['sgscore1']
    return ssdistnr, distnr, sgscore


if __name__=="__main__":
    # The very first observation of this field: 58345.44143477606
    # The last observation of this field before the break: 58386.4729638
    search_start = 2458345.94143478
    search_end = 2458386.9729638

    s = logon()
    #q = get_query(search_start, search_end)
    #r = s.query(query=q)
    #print(r['result']['query_result'])
    names = []
    dat = pandas.read_json("../docs/later_query.json")
    for item in dat['query_result']:
        name = item['objectId']
        names.append(name)
        #names.append(item['objectId'])
    names = np.unique(np.array(names))
    print(len(names))

    keep = []
    for ii,name in enumerate(names):
        print(ii)
        ssdistnr, distnr, sgscore = get_details(name)
        notrock = ssdistnr<0
        star = np.logical_and(distnr<1,sgscore>0.8)
        if np.logical_and(notrock, ~star):
            keep.append(name)

    np.savetxt("keep_2.txt", keep, dtype=str)
