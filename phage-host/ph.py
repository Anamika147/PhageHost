
import argparse
import math
#import numpy

import multiprocessing
from pathlib import Path
from typing import List, Set, Union, Dict

import numpy as np
import pandas as pd



def get_arguments():


    desc = f'Ph predicts hosts from phage (meta)genomic data'
    p = argparse.ArgumentParser(description=desc)
    p.add_argument('virus_dir', metavar='virus_dir')
    p.add_argument('host_dir', metavar='host_dir')
    p.add_argument('out', metavar='output_file', 
    type=argparse.FileType('w'))

    p.add_argument('--p', dest='p', type=float,
     default=0.75, metavar=None)

    p.add_argument('--k', dest='k', type=int,
                   default=30)
                   
    p.add_argument('--t', dest='num_threads', type=int,
                   default=multiprocessing.cpu_count())

    args = p.parse_args()
    return args


def rbo_out(hi, vi, p=0.75):
    hname = hnames[hi]
    vname = vnames[vi]

    return hi, vi, rbo(hd[hname], vd[vname], p=p)

def rbo(l1: List[Set[Union[str, int]]],
        l2: List[Set[Union[str, int]]],
        p: float = 0.75) -> float:
   

    # Finding short (S) and long (L) list
    sl, ll = sorted([(len(l1), l1), (len(l2), l2)])
    #s=len(sl)
    #l=len(ll)
    s, S = sl
    l, L = ll


    # ss= items in shortest list
    #ls=items in larget list
    ss = set()
    ls = set()

    x_d = {0: 0} 
    # no of shared items between S and L 
    # no of items in S
    len_ss = {0: 0}
    #no of items in L
    len_ls = {0: 0}

    sum1 = 0
    
    for i in range(l):
        d = i + 1

        xset = L[i]
        yset = S[i] if i < s else None

        ls.update(xset)
        if yset is not None:
            ss.update(yset)

        x_d[d] = x_d[d-1]   
        len_ss[d] = len_ss[d-1]
        len_ls[d] = len_ls[d-1]


        for x in xset:
            x_d[d] += 1 if x in ss else 0
            len_ls[d] += 1
        
        
        if yset is not None:
            for y in yset:
                x_d[d] += 0 if y not in ls or y in xset else 1
                len_ss[d] += 1
        

        agreement = (2 * x_d[d] / (len_ss[d] + len_ls[d]))
        sum1 += p ** d * agreement

    
    
    #rbo score calculation
    agreement = (2 * x_d[s] / (len_ss[s] + len_ls[s]))
    x_s = agreement * s
    
    sum2 = sum([p ** d * x_s * (d - s) / s / d for d in range(s+1, l+1)])
    
    term1 = (1 - p) / p * (sum1 + sum2)

    agreement = (2 * x_d[l] / (len_ss[l] + len_ls[l]))
    x_l = agreement * s
    
    term2 = p ** l * ((x_l - x_s) / l + x_s / s)
    
    return term1 + term2


"""def weight(d: int, p: float) -> float:
    
    p1 = 1 - pow(p, d-1)
    p2 = ((1-p)/p)*d
    p3 = math.log(1.0/(1-p))
    p4 = sum([pow(p, i)/i for i in range(1, d)])
 
    weight = p1 + p2 * (p3-p4)
 
    return weight"""


def read_list(filename: Union[str, Path], k: int = 0) -> List[Set[str]]:
    
      #read a list from a file
    fh = open(filename)
    lst = [set(l.strip().split(',')) for l in fh]
    fh.close()

    return lst[:k] if k else lst


def get_lists(directory: str, k: int = 0) -> Dict[str, List[Set[str]]]:
    #read all lists from a directory
    path = Path(directory)
    d = {}
    for f in path.iterdir():
        d[f.stem] = read_list(f, k)
    return d


if __name__ == '__main__':
    args = get_arguments()

    vd = get_lists(args.virus_dir, args.k)
    hd = get_lists(args.host_dir, args.k)

    vnames = sorted(vd)
    hnames = sorted(hd)

    
    vds = []
    for i, vname in enumerate(vnames):
        vs = set()
        for rset in vd[vname]:
            for item in rset:
                vs.add(item)
        vds.append(vs)

    q = []
    for hi, hname in enumerate(hnames):
        hs = set([item for rset in hd[hname] for item in rset])
        for vi, vs in enumerate(vds):
            if vs.intersection(hs):
                q.append((hi, vi, args.p))

      #calculate rbo scores between every phage and every host
    data = np.zeros((len(hnames), len(vnames)))
    with multiprocessing.Pool(args.num_threads) as pool:
        res = pool.starmap(rbo_out, q)
        for hi, vi, score in res:
            data[hi, vi] = score

     
    df = pd.DataFrame(data=data, index=hnames, columns=vnames)
    df.to_csv(f'{args.out.name}.matrix.csv')


    #for each phage find host with the highest score
    args.out.write(f'phage,host,score\n')
    for vname in vnames:
        max_score = df[vname].max()
        host_names = df.loc[df[vname] == max_score].index

        for hname in host_names:
            args.out.write(f'{vname},{hname},{max_score}\n')

    args.out.close()