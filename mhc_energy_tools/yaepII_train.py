#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Trains Yet Another Epitope Predictor (YAEP), class II
Further description in yaepII.py # TODO propagate

@author: Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu 
"""

import argparse, csv, os, shutil
from epilib.yaepII import BindingData, Cores, YAEPIITrainer
from epilib.netmhcII import NetMHCII

def setup_parser():
    """Creates an ArgumentParser and sets up all its arguments.
    returns: ArgumentParser"""

    parser = argparse.ArgumentParser(description='Trains Yet Another Epitope Predictor (YAEP), class II', 
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog="""Additional notes:
    - architecture for Andreatta 2015: 10,15,40,60 at 10 reps
    - specified allele name(s) must match the training data
    - predefined allele sets follow NetMHCII naming
    - can use fixed cores for training (faster during basic development), rather than iterative predicting using current model and then updating weights according to best core. # TODO: package up scripts to make this usable, if seems of general utility, else drop?""")
    # nn parameters
    parser.add_argument('--units', help='comma-separated list of different numbers of hidden units (default: %(default))', default='10,15,40,60')
    parser.add_argument('--reps', help='number of repetitions in ensemble (default: %(default))', type=int, default=10)
    parser.add_argument('--max_epochs', help='max number of epochs to train (default: %(default))', type=int, default=250)
    parser.add_argument('--validation_split', help='portion of training data to set aside for early stopping (default: %(default)); use 0 for no early stopping', type=float, default=0.15)
    # TODO: specs for encoding
    # TODO: other specs for training
    
    # which allele(s)
    alleles = parser.add_mutually_exclusive_group()
    alleles.add_argument('--allele_set', help='name of predefined set of alleles', choices=NetMHCII.allele_sets.keys())
    alleles.add_argument('--alleles', help='comma-separated list of allele names')
    
    # spec for training files
    training = parser.add_mutually_exclusive_group(required=True)
    training.add_argument('--netmhcii_folds_dir', help='holds the training data split into folds, in files test1.txt ... test5.txt, from http://www.cbs.dtu.dk/suppl/immunology/NetMHCIIpan-3.2/ (or other data in same format)')
    
    # spec for cores
    cores = parser.add_mutually_exclusive_group()
    cores.add_argument('--variable_cores', action='store_true') # TODO: this is default, but not clear how to specify that
    cores.add_argument('--netmhcii_predicted_cores_dir', help="use fixed pre-predicted cores for training (doesn't much impact accuracy, just speed of training); see notes for details")
        
    # outputs
    parser.add_argument('--out_base', help='base directory for models. etc. (nested within) (default: %(default))', default='models')
    parser.add_argument('--overwrite', help='must set this in order to overwrite existing model', action='store_true')

    return parser

def train_one_allele(args, allele):
    """Sets up an epitope predictor training run based on the parsed args.
    args: ArgumentParser""" 

    # output directory for this allele   
    base = args.out_base if args.out_base is not None else 'yaepII'
    allele_dir = base+'/'+allele+'/'
    if os.path.exists(allele_dir):
        if args.overwrite:
            print('overwriting',allele)
            shutil.rmtree(allele_dir)
        else:
            raise Exception('output directory "'+allele_dir+'" already exists; specify --overwrite if you want to overwrite it')
    elif not os.path.exists(base):
        os.mkdir(base)

    # architecture
    # TODO: error checking
    units = [int(u) for u in args.units.split(',')]

    yaep = YAEPIITrainer(units=units, reps=args.reps, max_epochs=args.max_epochs, validation_split=args.validation_split)

    # training data
    if args.netmhcii_folds_dir is not None:
        binding_data = BindingData.load_nm_training([args.netmhcii_folds_dir+'/test'+str(fold+1)+'.txt' for fold in range(5)], allele)
        if any(len(fold)==0 for fold in binding_data.folds):
            raise Exception('empty fold for allele '+allele)
        print('loaded folds of size', [len(fold) for fold in binding_data.folds])
    fixed_cores = None
    if args.netmhcii_predicted_cores_dir is not None:
        nm_cores = Cores.load_nm_predictions([args.netmhcii_predicted_cores_dir+'/'+allele+'-fold'+str(fold+1)+'.out' for fold in range(5)], allele)
        fixed_cores = nm_cores.cores
        print('loaded',len(fixed_cores),'cores')
           
    # train it
    yaep.xval_train(binding_data, fixed_cores=fixed_cores)

    # outputs
    # model
    yaep.save_model(allele_dir) # note that this creates the directory, which must not exist previously
    # frozen model
    cmd = ' '.join(['freeze_graph', '--input_saved_model_dir', allele_dir, '--output_node_names','prediction,main_init', '--initializer_nodes','main_init', '--output_graph',allele_dir+'/frozen.pb'])
    print('>', cmd)
    os.system(cmd)
    # TODO: errors?
    
    # predictions during training
    flat_predicted = list(yaep.predicted.values())
    flat_actual = [binding_data.values[p] for p in yaep.predicted]
    with open(allele_dir+'summary.txt','w') as outfile:
        print(allele, 'overall xval: auc %f; f1 %f; pcc %f' % yaep.evaluate_prediction(flat_predicted, flat_actual), file=outfile)

    with open(allele_dir+'xval-pred.csv','w') as outfile:
        outcsv = csv.writer(outfile)
        outcsv.writerow(['fold','peptide','actual','predicted'])
        for (fold, fold_peptides) in enumerate(binding_data.folds):
            for peptide in fold_peptides:
                outcsv.writerow([fold,peptide,binding_data.values[peptide],yaep.predicted[peptide]])  
    
def main(args):
    if args.allele_set is not None:
        # note that the allowed choices were limited to the dictionary keys, so should be legit
        alleles = NetMHCII.allele_sets[args.allele_set]
    elif args.alleles is not None:
        alleles = args.alleles.split(',')
    if len(alleles) == 0:
        print('!!! no alleles specified -- this will be easy :)')
        
    for allele in alleles:
        print('=================',allele)
        train_one_allele(args, allele)
    
# python train_yaepII.py --netmhcii_folds_dir training-data/netmhcIIpan-3.2/folds/ --alleles DRB1_0101
# python train_yaepII.py --netmhcii_folds_dir training-data/netmhcIIpan-3.2/folds/ --netmhcii_predicted_cores_dir training-data/netmhcIIpan-3.2/pred/ --alleles DRB1_0101
if __name__ == '__main__':
    parser = setup_parser()
    args = parser.parse_args()
    main(args)
