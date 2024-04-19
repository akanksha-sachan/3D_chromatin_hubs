#!/usr/bin/env python

import os, joblib
import numpy as np
import argparse
import scoreUtils, utils
import cooler

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
    else:
        hic = True

    if not hic:
        Lib = cooler.Cooler(args.path)

    # ccname is consistent with chromosome labels in .hic / .cool
    ccname = args.chrom
    cikada = 'chr' + ccname.lstrip('chr')  # cikada always has prefix "chr"

    if not hic:
        if correct:
            M = utils.tocsr(Lib.matrix(balance=correct, sparse=True).fetch(ccname))
            raw_M = utils.tocsr(Lib.matrix(balance=False, sparse=True).fetch(ccname))
            weights = Lib.bins().fetch(ccname)[correct].values
            X = scoreUtils.Chromosome(M, model=model, raw_M=raw_M, weights=weights,
                                      cname=cikada, lower=args.lower,
                                      upper=args.upper, res=args.resolution,
                                      width=width)
        else:
            M = utils.tocsr(Lib.matrix(balance=False, sparse=True).fetch(ccname))
            X = scoreUtils.Chromosome(M, model=model, raw_M=M, weights=None,
                                      cname=cikada, lower=args.lower,
                                      upper=args.upper, res=args.resolution,
                                      width=width)
    else:
        if correct:
            M = utils.csr_contact_matrix('KR', args.path, ccname, ccname, 'BP', args.resolution)
            raw_M = utils.csr_contact_matrix('NONE', args.path, ccname, ccname, 'BP', args.resolution)
            X = scoreUtils.Chromosome(M, model=model, raw_M=raw_M, weights=None,
                                      cname=cikada, lower=args.lower,
                                      upper=args.upper, res=args.resolution,
                                      width=width)
        else:
            M = utils.csr_contact_matrix('NONE', args.path, ccname, ccname, 'BP', args.resolution)
            X = scoreUtils.Chromosome(M, model=model, raw_M=M, weights=None,
                                      cname=cikada, lower=args.lower,
                                      upper=args.upper, res=args.resolution,
                                      width=width)
    
    result, R = X.score(thre=args.minimum_prob)
    X.writeBed(args.output, result, R)

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--model', help='Path to the trained model')
    parser.add_argument('--path', help='Path to the input Hi-C file')
    parser.add_argument('--chrom', help='Chromosome to score')
    parser.add_argument('--output', help='Output BEDfile path')
    parser.add_argument('--clr_weight_name', help='Normalization method', type=str, default='raw')
    parser.add_argument('--resolution', help='Resolution of the Hi-C file', type=int, default=10000)
    parser.add_argument('--lower', help='Lower bound of the anchor distance', type=int, default=6)
    parser.add_argument('--upper', help='Upper bound of the anchor distance', type=int, default=300)
    parser.add_argument('--minimum_prob', help='Minimum probability to consider', type=float, default=0.5)
    args = parser.parse_args()
    main(args)