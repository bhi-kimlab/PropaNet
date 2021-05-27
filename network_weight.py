import pandas as pd
import argparse
import numpy as np
from scipy.stats import pearsonr
import multiprocessing


def isNum(x):
    try:
        float(x)
        return True
    except:
        return False


def myCorr(x):
    g1, g2 = sorted(x)
    if g1 == g2:
        val = 0.0
    else:
        val, pval = pearsonr(lst_exps[g1], lst_exps[g2])
    return (g1, g2, val)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="python %(prog)s -nwk nwkFile -exp expFile -o out"
    )
    parser.add_argument("-nwk", required=True, help="Network file")
    parser.add_argument("-exp", required=True, help="gene expression File")
    parser.add_argument("-p", type=int, metavar="N", default="1", help="total process")
    parser.add_argument("-o", required=True, help="out File")
    args = parser.parse_args()

    lst_exps = dict()
    with open(args.exp) as f:
        lines = f.readlines()
    for line in lines:
        s = line.strip().split("\t")
        if not isNum(s[1]):
            continue
        else:
            gene, exps = s[0], list(map(float, s[1:]))
            lst_exps[gene] = exps
    lst_pairs = list()
    with open(args.nwk) as f2:
        for line in f2.readlines():
            s = line.strip().split()
            if s[0] not in lst_exps or s[1] not in lst_exps:
                continue
            lst_pairs.append([s[0], s[1]])
    pool = multiprocessing.Pool(args.p)
    res = pool.imap_unordered(myCorr, lst_pairs)
    with open(args.o, "w") as f3:
        for g1, g2, val in res:
            if g1 == g2:
                continue
            f3.write("\t".join([g1, g2, str(val)]) + "\n")

    print("network done")
