""" Count the number of unique sources """
import glob
import numpy as np

def calc(search):
    allf = glob.glob(search)

    cands = np.loadtxt(allf[0], dtype=str)
    for f in allf:
        #print(f)
        newcands = np.loadtxt(f, dtype=str)
        #print(newcands)
        cands = np.append(cands, newcands)

    ucands = np.unique(cands)
    np.savetxt("passed_lc_check.txt", ucands, fmt='%s')
    return len(allf), len(ucands)


if __name__=="__main__":
    #print(calc("field*filter1.txt"))
    #print(calc("field*nostars.txt"))
    print(calc("field*goodlc.txt"))
