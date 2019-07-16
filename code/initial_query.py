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
import glob
from kowalski_login import logon


if __name__=="__main__":
    s = logon()
    allfields = np.loadtxt("fieldlist.txt")

    for field in allfields:
        print("running field %s" %field)
        outf = "recent_field%s_filter1.txt" %field
        if len(glob.glob(outf))==0:
            q = {"query_type": "find",
                 "query": {
                     "catalog": "ZTF_alerts",
                     "filter": {
                             'candidate.jd': {'$gt': 2458585, '$lt': 2458682},
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

            np.savetxt(outf, names, fmt='%s')
