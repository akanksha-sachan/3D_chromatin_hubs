#!/usr/bin/env python

import joblib, os
import numpy as np
import scoreUtils, utils
import cooler
import argparse

def main(args):
    
    np.seterr(divide='ignore', invalid='ignore')

    if os.path.exists(args.output):
        os.remove(args.output)

    model = joblib.load(args.model)

    # decide which normalization method to use
    if args.clr_weight_name.lower() == 'raw':
        correct = False
    else:
        correct = args.clr_weight_name

    # deduce the width parameter used during the training
    width = int((np.sqrt(model.feature_importances_.size) - 1) / 2)

    # more robust to check if a file is .hic
    hic_info = utils.read_hic_header(args.path)
    if hic_info is None:
        hic = False
        Lib = cooler.Cooler(args.path)
        chromosomes = Lib.chromnames[:]
        #nam = args.path.split('.cool')[0]
    else:
        hic = True
        chromosomes = utils.get_hic_chromosomes(args.path, args.resolution)
        #nam = args.path.split('.hic')[0]
    #nam = nam.split('/')[-1]

    queue = []
    for key in chromosomes:

        chromlabel = key.lstrip('chr')
        if (not args.chroms) or (chromlabel.isdigit() and '#' in args.chroms) or (chromlabel in args.chroms):
            queue.append(key)
    
    for key in queue:

        if key.startswith('chr'):
            cname = key
        else:
            cname = 'chr'+key
            
        if not hic:
            if correct:
                M = utils.tocsr(Lib.matrix(balance=correct, sparse=True).fetch(key))
                raw_M = utils.tocsr(Lib.matrix(balance=False, sparse=True).fetch(key))
                weights = Lib.bins().fetch(key)[correct].values
                X = scoreUtils.Chromosome(M, model=model, raw_M=raw_M, weights=weights,
                                          cname=cname, lower=args.lower,
                                          upper=args.upper, res=args.resolution,
                                          width=width)
            else:
                M = utils.tocsr(Lib.matrix(balance=False, sparse=True).fetch(key))
                X = scoreUtils.Chromosome(M, model=model, raw_M=M, weights=None,
                                          cname=cname, lower=args.lower,
                                          upper=args.upper, res=args.resolution,
                                          width=width)
        else:
            if correct:
                M = utils.csr_contact_matrix('KR', args.path, key, key, 'BP', args.resolution)
                raw_M = utils.csr_contact_matrix('NONE', args.path, key, key, 'BP', args.resolution)
                X = scoreUtils.Chromosome(M, model=model, raw_M=raw_M, weights=None,
                                          cname=cname, lower=args.lower,
                                          upper=args.upper, res=args.resolution,
                                          width=width)
            else:
                M = utils.csr_contact_matrix('NONE', args.path, key, key, 'BP', args.resolution)
                X = scoreUtils.Chromosome(M, model=model, raw_M=M, weights=None,
                                          cname=cname, lower=args.lower,
                                          upper=args.upper, res=args.resolution,
                                          width=width)

        result, R = X.score(thre=args.minimum_prob)
        X.writeBed(args.output, result, R)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--model', type=str, help='path to the trained model')
    parser.add_argument('--path', type=str, help='path to the Hi-C data')
    parser.add_argument('--output', type=str, help='output file')
    parser.add_argument('--resolution', type=int, default=10000, help='resolution of the Hi-C data')
    parser.add_argument('--clr_weight_name', type=str, default='KR', help='normalization method')
    parser.add_argument('--lower', type=int, default=6, help='lower bound of the diagonal')
    parser.add_argument('--upper', type=int, default=300, help='upper bound of the diagonal')
    parser.add_argument('--chroms', type=str, default=None, help='chromosome names')
    parser.add_argument('--minimum_prob', type=float, default=0.5, help='minimum probability of a loop')
    args = parser.parse_args()
    main(args)