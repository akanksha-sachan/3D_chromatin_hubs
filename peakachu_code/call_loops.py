#!/usr/env/bin python

import gc
import sys
import numpy as np
import pandas as pd
from collections import defaultdict
import argparse
import peakacluster

def main(args):

    res = args.resolution
    
    clusters, score_pool = peakacluster.parse_peakachu(args.infile, args.threshold, res)
    with open(args.outfile, 'w') as out:
        for c in clusters:
            for p in clusters[c]:
                if p in score_pool[c]:
                    s1 = str(p[0]*res)
                    e1 = str(p[0]*res+res)
                    s2 = str(p[1]*res)
                    e2 = str(p[1]*res+res)
                    prob = str(score_pool[c][p][0])
                    raw_signal = str(score_pool[c][p][1])
                    line = [c, s1, e1, c, s2, e2, prob, raw_signal]
                    out.write('\t'.join(line)+'\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--threshold', help='Path to the trained model', type=float, default=0.9)
    parser.add_argument('--infile', help='Path to the input Hi-C file')
    parser.add_argument('--outfile', help='Output filtered loops BEDfile path')
    parser.add_argument('--resolution', help='Resolution of the Hi-C file', type=int, default=10000)
    args = parser.parse_args()
    main(args)