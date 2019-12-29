import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time


class Position:
    """ Get a celestial position in various formats """
    def __init__(self, ra, dec, rastr, decstr):
        """ RA and Dec given in decimal degrees """
        self.coord = SkyCoord(ra,dec,unit='deg')
        self.ra = ra
        self.dec = dec
        self.rastr = rastr
        self.decstr = decstr
        # In Galactic coordinates as well (to see latitude)
        self.l = self.coord.galactic.l.deg
        self.b = self.coord.galactic.b.deg



def get_lc(s, name):
    """ Retrieve LC for object """
    det_alerts = get_dets(s, name)
    
    # If the detection resulted in an alert
    isalert = []
        
    # Detections
    jd = []
    mag = []
    emag = []
    filt = []
    program = []

    # Non-detections
    limjds = []
    limmags = []
    limfilts =[]
    limprogram = []
            
    for det in det_alerts:
        cand = det['candidate']

        current_jd = cand['jd']
        current_mag = cand['magpsf']
        current_emag = cand['sigmapsf']
        current_filter = cand['fid'] 
        current_prog = cand['programid']
                
        if current_jd not in jd:
            isalert.append(True)
            jd.append(current_jd)
            mag.append(current_mag)
            emag.append(current_emag)
            filt.append(current_filter)
            program.append(current_prog)

    det_prv = get_prv_dets(s, name)
    for prv_cand in det_prv:
        if 'magpsf' in prv_cand.keys():
            if prv_cand['jd'] not in jd:
                isalert.append(False)
                jd.append(prv_cand['jd'])
                mag.append(prv_cand['magpsf'])
                emag.append(prv_cand['sigmapsf'])
                filt.append(prv_cand['fid'])
                program.append(prv_cand['pid'])
        else:
            limjds.append(prv_cand['jd'])
            limmags.append(prv_cand['diffmaglim'])
            limfilts.append(prv_cand['fid'])
            limprogram.append(prv_cand['pid'])

    isalert = np.array(isalert)
    jd = np.array(jd)
    mag = np.array(mag)
    emag = np.array(emag)
    filt = np.array(filt)
    program = np.array(program)
    limjds = np.array(limjds)
    limmags = np.array(limmags)
    limfilts = np.array(limfilts)
    limprogram = np.array(limprogram)

    # Sort in order of jd
    order = np.argsort(jd)
    isalert = isalert[order]
    jd = jd[order]
    mag = mag[order]
    emag = emag[order]
    filt = filt[order]
    program = program[order]

    order = np.argsort(limjds)
    limjds = limjds[order]
    limmags = limmags[order]
    limfilts = limfilts[order]
    limprogram = limprogram[order]

    return isalert,jd,mag,emag,filt,program,limjds,limmags,limfilts,limprogram


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

