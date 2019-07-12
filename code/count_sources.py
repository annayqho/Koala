""" Count the number of unique sources """
import glob
import numpy as np

def calc(search):
    allf = glob.glob(search)

    cands = np.loadtxt(allf[0], dtype=str)
    for f in allf:
        cands = np.append(cands, np.loadtxt(f, dtype=str))

    ucands = np.unique(cands)
    return len(allf), len(ucands)


if __name__=="__main__":
    print(calc("field*filter1.txt"))
    print(calc("field*nostars.txt"))
