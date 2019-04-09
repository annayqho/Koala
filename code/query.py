""" Query for fast-evolving transients

Criteria:
    - more than 10 detections at a given position with RB score > 0.5
    - rising or fading at a rate >= 1 mag/day
    - jdendhist - jdstarthst < 1 year
"""

from kowalski_login import logon


if __name__=="__main__":
    s = logon()
    q = {"query_type": "general_search", 
         "query": "db['ZTF_alerts'].aggregate([{'$group':{'_id':'$objectId', 
                                                          'count':{'$sum':{'$cond':[{'$gt':['$candidate.rb',0.5]}, 1, 0]} } } }, 
                                               {'$match':{'count':{'$gt':20}}}, 
                                               {'$project':{'objectId':'$_id','_id':0}}], allowDiskUse=True)" }
    r = s.query(query=q)
    print(r['result']['query_result'])
