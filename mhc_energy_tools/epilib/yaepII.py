#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Yet Another Epitope Predictor (YAEP), class II

Tensorflow based implementation of a standard neural network based MHC class II epitope predictor.
Currently just a pretty straight-up encoding of the architecture and peptide encoding described in
  Improved methods for predicting peptide binding affinity to MHC class II molecules. 
  Jensen KK, Andreatta M, Marcatili P, Buus S, Greenbaum JA, Yan Z, Sette A, Peters B, and Nielsen M. 
  PMID: 29315598
But general enough to readily enable variations in all aspects.

TODO: this is all tested using tf 1.13, and is generating deprecation warnings for tf 2
but since tf 2 is still in beta, not yet addressing that

@author: Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu 
"""

import csv, math, os
import tensorflow as tf
import numpy as np, scipy.stats, sklearn.metrics

from epilib.epitope_predictor import EpitopePredictor, EpitopeScore

# ============================================================================
# training models

AAs = 'ACDEFGHIKLMNPQRSTVWY'
AAs21 = 'ACDEFGHIKLMNPQRSTVWYX' # X shows up in some of the training data
AApad = '-' # a character not actually used in a peptide, to pad the N- and C-terminal flanks to a uniform length of 3
Flank = AApad+AApad+AApad # an empty flank

# ------------------------------------------------------------------------

class BindingData (object):
    """Holds MHC:peptide binding data, potentially split into folds.
       - name: string, arbitrary name for the dataset
       - values: { string : float }, mapping from peptides to binding probabilities
       - folds: [[ string ]], partitioning the peptides into different folds (if not given, all in one fold)"""

    nM_base = 50000 # for converting IC50 in nM into probability, following NetMHCII convention
    
    @staticmethod
    def nM_to_prob(nM):
        """Converts an IC50 in nM to a probability by way of logging to the nM_base"""
        return max(0,1-math.log(nM,BindingData.nM_base))
    
    @staticmethod
    def prob_to_nM(prob):
        """Converts a probability to an IC50 in nM by way of exponentiating the nM_base"""
        return math.pow(BindingData.nM_base,1-prob)

    def __init__(self, values, folds=None):
        self.values = values
        self.folds = folds if folds is not None else [ list(values.keys()) ]
        
    @classmethod
    def load_iedb_benchmark(cls, filenames):
        """Loads benchmark data in the IEDB data from the specified files
        #TODO: url
        filenames: list of the names, one per fold
        returns: new instance of the class"""
        values = {}
        folds = []
        for (fold, filename) in enumerate(filenames):
            fold_peptides = []
            folds.append(fold_peptides)
            with open(filename, 'r') as infile:
                for row in csv.DictReader(infile, fieldnames=['species','allele','len','tbd','sequence','inequality','meas'], delimiter='\t'):
                    # TODO: ignoring inequality; could change model to leverage it (one-sided loss)
                    peptide = row['sequence']
                    if len(peptide)<9:
                        print('!!! skipping short',peptide)
                        continue
                    value = float(row['meas'])
                    if peptide in values:
                        print('!!! ignoring duplicate entry for',peptide,':',values[peptide],'vs.',value)
                    else:
                        values[peptide] = BindingData.nM_to_prob(value)
                        fold_peptides.append(peptide)
        return cls(values, folds)

    @classmethod
    def load_nm_training(cls, filenames, allele):
        """Loads training data in the NetMHCIIPan format for the specified allele from the specified files
        Ex: http://www.cbs.dtu.dk/suppl/immunology/NetMHCIIpan-3.2/train[1-5]
        filenames: list of the names, one per fold
        allele: name of the allele (since NetMHCIIPan includes all alleles in a single file)
        returns: new instance of the class"""
        values = {}
        folds = []
        for (fold, filename) in enumerate(filenames):
            fold_peptides = []
            folds.append(fold_peptides)
            for row in open(filename, 'r'):
                (peptide,value,allele2) = row.split()
                if allele2 != allele: continue
                if len(peptide)<9:
                    print('!!! skipping short',peptide)
                elif peptide in values:
                    print('!!! ignoring duplicate entry for',peptide,':',values[peptide],'vs.',value)
                else:
                    values[peptide] = float(value)
                    fold_peptides.append(peptide)
        return cls(values, folds)
        
    @classmethod
    def load_iedb_predictions(cls, filename):
        """Loads the predicted values from the specified IEDB benchmark file (note: single file contains predictions for all folds)
        #TODO: url
        returns: new instance of the class"""
        values = {}
        # TODO: allow copying folds from training data instance?
        with open(filename, 'r') as infile:
            for row in csv.DictReader(infile, delimiter='\t'):
                peptide = row['peptide']
                if len(peptide)<9:
                    print('!!! skipping short',peptide)
                    continue
                value = float(row['IC50_p'])
                if peptide in values:
                    print('!!! ignoring duplicate entry for',peptide,':',values[peptide],'vs.',value)
                else:
                    values[peptide] = value
        return cls(values)
    
    @classmethod
    def load_nm_predictions(cls, filenames, allele):
        """Loads predictions in the NetMHCIIPan format for the specified allele from the specified files
        #TODO: how generated
        filenames: list of the names, one per fold
        allele: name of the allele (in case the file contains predicitons for multiple alleles)
        returns: new instance of the class"""
        values = {}
        folds = []
        for (fold, filename) in enumerate(filenames):
            fold_peptides = []
            folds.append(fold_peptides)
            for row in open(filename,'r'):
                if len(row)==0 or row[0]=='#': continue
                cols = row.split()
                if len(cols)==0 or cols[0] != allele: continue
                peptide = cols[2]
                value = float(cols[5])
                if peptide in values:
                    print('!!! ignoring duplicate entry for',peptide,':',values[peptide],'vs.',value)
                    continue
                values[peptide] = value
                fold_peptides.append(peptide)
        return cls(values, folds)

# --------------------

class Cores (object):
    """Holds a specification of binding cores for peptides, to expedite training during development"""

    def __init__(self, cores):
        """cores: { string : int }, the core index for each peptide (i.e., the core 9mer starts at that position in the peptide)"""
        self.cores = cores

    @classmethod
    def load_iedb_predictions(cls, filename):
        """Loads the predicted cores from the specified IEDB benchmark file"""
        cores = {}
        with open(filename, 'r') as infile:
            for row in csv.DictReader(infile, delimiter='\t'):
                peptide = row['peptide']
                if len(peptide)<9:
                    print('!!! skipping short',peptide)
                    continue
                if peptide in cores:
                    print('!!! ignoring duplicate entry for',peptide)
                else:
                    core = int(row['core'])
                    if core<0:
                        print('!!! negative core for',peptide,'; setting to 0')
                        core = 0
                    elif core+9>len(peptide):
                        print('!!! core at',core,'for',peptide,'does not leave a 9-mer; stepping back to', len(peptide)-9)
                        core = len(peptide)-9
                    cores[peptide] = [core]

        return cls(cores)
    
    @classmethod
    def load_nm_predictions(cls, filenames, allele):
        """Loads the predicted cores for the specified allele from the specified NetMCHIIPan predictions files"""
        cores = {}
        for filename in filenames:
            for row in open(filename,'r'):
                if len(row)==0 or row[0]=='#': continue
                cols = row.split()
                if len(cols)==0 or cols[0] != allele: continue
                peptide = cols[2]
                if peptide in cores:
                    print('!!! ignoring duplicate entry for',peptide)
                    continue
                core = int(cols[4])
                if core<0:
                    print('!!! negative core for',peptide,'; setting to 0')
                    core = 0
                elif core+9>len(peptide):
                    print('!!! core at',core,'for',peptide,'does not leave a 9-mer; stepping back to', len(peptide)-9)
                    core = len(peptide)-9
                cores[peptide] = [core]
        return cls(cores)
    
# ------------------------------------------------------------------------

class Subst (object):
    """A substitution matrix
    aa_vectors: { char : [ int ] }, for each AA, a vector of substitution scores
    col_aas : [ char ], the AAs corresponding to those substitution scores
    """

    def __init__(self, aa_vectors, col_aas):
        self.aa_vectors = aa_vectors
        self.col_aas = col_aas
        self.width = len(self.col_aas)

    def encoding(self, aa):
        """The vector for the AA"""
        return self.aa_vectors[aa]      
    
    @classmethod
    def parse_blast(cls, raw, row_aas=AAs21, col_aas=AAs):
        """Creates a substitution matrix from the raw string format as used in BLAST:
        first row labels all the AAs in the columns
        remaining rows give an AA and its substitution scores"""
        rows = [row.strip().split() for row in raw.strip().split('\n')]
        header_indices = dict(zip(rows[0], range(len(rows[0]))))
        col_indices = [1+header_indices[aa] for aa in col_aas] # header starts at 0, rest start at 1
        aa_vectors = dict((row[0], np.asarray([float(row[i]) for i in col_indices])) for row in rows[1:] if row[0] in row_aas)
        return cls(aa_vectors, col_aas)

    @classmethod
    def make_hot(cls, row_aas=AAs21, col_aas=AAs, on=1, off=0, ambig=0.05):
        """A one hot encoding, potentially relaxed by not being binary for "on" and "off",
        and allowing a different value for the "ambiguous" AA 'X'"""
        aa_vectors = dict((aa, np.asarray([on if aa2==aa else off for aa2 in col_aas])) for aa in row_aas)
        if ambig is not None:
            aa_vectors['X'] = np.asarray([ambig for aa in col_aas])
        return cls(aa_vectors, col_aas)
    
# TODO: just stuffed this in here instead of reading from a file
Subst.blosum50 = Subst.parse_blast("""
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *
A  5 -2 -1 -2 -1 -1 -1  0 -2 -1 -2 -1 -1 -3 -1  1  0 -3 -2  0 -2 -2 -1 -1 -5
R -2  7 -1 -2 -4  1  0 -3  0 -4 -3  3 -2 -3 -3 -1 -1 -3 -1 -3 -1 -3  0 -1 -5
N -1 -1  7  2 -2  0  0  0  1 -3 -4  0 -2 -4 -2  1  0 -4 -2 -3  5 -4  0 -1 -5
D -2 -2  2  8 -4  0  2 -1 -1 -4 -4 -1 -4 -5 -1  0 -1 -5 -3 -4  6 -4  1 -1 -5
C -1 -4 -2 -4 13 -3 -3 -3 -3 -2 -2 -3 -2 -2 -4 -1 -1 -5 -3 -1 -3 -2 -3 -1 -5
Q -1  1  0  0 -3  7  2 -2  1 -3 -2  2  0 -4 -1  0 -1 -1 -1 -3  0 -3  4 -1 -5
E -1  0  0  2 -3  2  6 -3  0 -4 -3  1 -2 -3 -1 -1 -1 -3 -2 -3  1 -3  5 -1 -5
G  0 -3  0 -1 -3 -2 -3  8 -2 -4 -4 -2 -3 -4 -2  0 -2 -3 -3 -4 -1 -4 -2 -1 -5
H -2  0  1 -1 -3  1  0 -2 10 -4 -3  0 -1 -1 -2 -1 -2 -3  2 -4  0 -3  0 -1 -5
I -1 -4 -3 -4 -2 -3 -4 -4 -4  5  2 -3  2  0 -3 -3 -1 -3 -1  4 -4  4 -3 -1 -5
L -2 -3 -4 -4 -2 -2 -3 -4 -3  2  5 -3  3  1 -4 -3 -1 -2 -1  1 -4  4 -3 -1 -5
K -1  3  0 -1 -3  2  1 -2  0 -3 -3  6 -2 -4 -1  0 -1 -3 -2 -3  0 -3  1 -1 -5
M -1 -2 -2 -4 -2  0 -2 -3 -1  2  3 -2  7  0 -3 -2 -1 -1  0  1 -3  2 -1 -1 -5
F -3 -3 -4 -5 -2 -4 -3 -4 -1  0  1 -4  0  8 -4 -3 -2  1  4 -1 -4  1 -4 -1 -5
P -1 -3 -2 -1 -4 -1 -1 -2 -2 -3 -4 -1 -3 -4 10 -1 -1 -4 -3 -3 -2 -3 -1 -1 -5
S  1 -1  1  0 -1  0 -1  0 -1 -3 -3  0 -2 -3 -1  5  2 -4 -2 -2  0 -3  0 -1 -5
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  2  5 -3 -2  0  0 -1 -1 -1 -5
W -3 -3 -4 -5 -5 -1 -3 -3 -3 -3 -2 -3 -1  1 -4 -4 -3 15  2 -3 -5 -2 -2 -1 -5
Y -2 -1 -2 -3 -3 -1 -2 -3  2 -1 -1 -2  0  4 -3 -2 -2  2  8 -1 -3 -1 -2 -1 -5
V  0 -3 -3 -4 -1 -3 -3 -4 -4  4  1 -3  1 -1 -3 -2  0 -3 -1  5 -3  2 -3 -1 -5
B -2 -1  5  6 -3  0  1 -1  0 -4 -4  0 -3 -4 -2  0  0 -5 -3 -3  6 -4  1 -1 -5
J -2 -3 -4 -4 -2 -3 -3 -4 -3  4  4 -3  2  1 -3 -3 -1 -2 -1  2 -4  4 -3 -1 -5
Z -1  0  0  1 -3  4  5 -2  0 -3 -3  1 -1 -4 -1  0 -1 -2 -2 -3  1 -3  5 -1 -5
X -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -5
* -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5  1
""")
Subst.one_hot = Subst.make_hot()
Subst.relaxed_one_hot = Subst.make_hot(on=0.95, off=0.05)

# ------------------------------------------------------------------------

class PeptideEncoder (object):
    """Encodes a string peptide into a tensor of features to be input to a model"""

    def __init__(self, peptide, expt_len, encoding):
        self.peptide = peptide
        self.expt_len = expt_len
        self.encoding = encoding
        self.width = encoding.shape[1].value # TODO: always right?

    @classmethod
    def andreatti15(cls, raw_subst=Subst.blosum50, pad_char=AApad):
        """Encoding from Andreatti 2015 Immunogenetics"""

        # make an internal encoding of the substitution matrix
        # add the padding character to the internal representation of the substitution matrix
        pad = tf.constant(pad_char)
        aas = list(raw_subst.aa_vectors.keys()) + [pad_char]
        subst = np.array(list(raw_subst.aa_vectors.values()) + [[0]*len(raw_subst.col_aas)])
        # first from AA string tensor to an integer index
        table = tf.contrib.lookup.index_table_from_tensor(tf.constant(aas))
        # then from integer index to tensor from substitution matrix
        encoder = tf.constant(subst, dtype=tf.float32)

        # input is a 15-character string; can be part of a larger peptide in expertimental data for training
        pep15 = tf.placeholder(tf.string, shape=(None,), name='peptide')
        expt_len = tf.placeholder_with_default(15.0*tf.ones_like(pep15, dtype=tf.float32), shape=(None,))
        # split it into 15 1-character strings
        aas = tf.sparse.to_dense(tf.string_split(pep15, delimiter=''), default_value=pad_char)
        # extract the flanks and core
        nflank3 = aas[:,:3]
        core9 = aas[:,3:12]
        cflank3 = aas[:,12:15]

        # build a list of various types of features; will flatten later
        feats = []
        # substitution-matrix encoding of core 9mer
        core9_enc = tf.gather(encoder, table.lookup(core9))
        feats.append(tf.reshape(core9_enc, shape=(-1, 9*encoder.shape[1])))
        # encoding of flanks is the average of the substitution matrices
        # but ignore the pad character
        nflank_naa = tf.reduce_sum(tf.cast(tf.not_equal(nflank3, pad), tf.float32), 1)
        cflank_naa = tf.reduce_sum(tf.cast(tf.not_equal(cflank3, pad), tf.float32), 1)
        feats.append(tf.reduce_sum(tf.gather(encoder, table.lookup(nflank3)), 1) / tf.expand_dims(tf.math.maximum(1.0,nflank_naa), 1))
        feats.append(tf.reduce_sum(tf.gather(encoder, table.lookup(cflank3)), 1) / tf.expand_dims(tf.math.maximum(1.0,cflank_naa), 1))
        # length features
        expt_len_enc = 1/(1+tf.exp(expt_len-15)/2)
        feats.append(tf.stack([nflank_naa/(1+nflank_naa), 1-nflank_naa/(1+nflank_naa), cflank_naa/(1+cflank_naa), 1-cflank_naa/(1+cflank_naa), expt_len_enc, 1-expt_len_enc], axis=1))

        # now put it all together
        encoding = tf.concat(feats, axis=1, name='encoding')

        return cls(pep15, expt_len, encoding)

    def encode_peptides(self, peptides, sess, possible_cores=None):
        """Pre-encodes the peptides, broken down into their possible cores, for feeding into a model.
        peptides: [ string ], the variable length (>= 9) peptides
        possible_cores: { string : [ int ] }, the various indices into the peptide that should be considered as possible cores
           if None, then all possibilities are considered
        returns: { string: [ tensor ] }, for each peptide, the encodings according to different cores"""
        encoded = {}
        for peptide in peptides:
            # 15mers (padded as needed) for each possible core
            cores = possible_cores[peptide] if possible_cores is not None else range(len(peptide)-8)
            # by adding --- to the N-terminus, the core index is now at the start of the N-flank
            formatted = [ (Flank+peptide+Flank)[core : (core+15)] for core in cores ]
            encoded[peptide] = sess.run(self.encoding, feed_dict={ self.peptide:formatted, self.expt_len: [len(peptide)]*len(cores) })
        return encoded

# ------------------------------------------------------------------------

class YAEPIITrainer (object):
    """Handles the training of an ensemble of models
    Various parameters specify the encoding, the NN architecture, and the training setup
    The main method is train_xval, given a BindingData set split into folds for cross-validation
    The models for the different folds comprise the core ensemble.
    The ensemble is further magnified by different architectures (numbers of hidden units) 
    and multiple replicates for each architecture (different random initializations)"""
    
    def __init__(self,
                 peptide_encoder=None, 
                 units=None, reps=10, activation='tanh', 
                 optimizer='adam', patience=10, max_epochs=250, batch_size=25, validation_split=0.15, verbose=0,
                 binder_thresh=None):

        self.graph = tf.Graph()
        self.sess = tf.Session(graph=self.graph)
        tf.keras.backend.set_session(self.sess) # note: tf.keras actually has its own distinct session otherwise! caused much pain to figure that out

        with self.graph.as_default(): # TODO: feels a bit clunky, but don't see a cleaner way if localizing graph to this trainer
            self.peptide_encoder = peptide_encoder if peptide_encoder is not None else PeptideEncoder.andreatti15()

        self.units = units if units is not None else [40]
        self.reps = reps
        self.activation = activation
        self.optimizer = optimizer
        self.patience = patience
        self.max_epochs = max_epochs
        self.batch_size = batch_size
        self.validation_split = validation_split
        self.verbose = verbose
        self.binder_thresh = binder_thresh if binder_thresh is not None else BindingData.nM_to_prob(500)

    def build_ensemble(self, folds, base_name='yaepII'):
        """Creates a new ensemble according to the architecture and replicate setup"""
        self.ensemble = []
        self.fold_ensembles = []
        for fold in range(folds):
            fold_ensemble = []
            for size in self.units:
                for rep in range(self.reps):
                    model = tf.keras.models.Sequential(name=base_name+'-%d-%d-%d' % (fold, size, rep))
                    model.add(tf.keras.layers.Dense(size, activation=self.activation, input_shape=(self.peptide_encoder.width,)))
                    model.add(tf.keras.layers.Dense(1, activation='sigmoid'))
                    model.compile(optimizer=self.optimizer, loss='mse')
                    fold_ensemble.append(model)
            self.fold_ensembles.append(fold_ensemble)
            self.ensemble.extend(fold_ensemble)
        
    def decode_predictions(self, peptides, encodings, predictions):
        """Multiple cores can be tested for each peptide, but in this approach only one matters.
        So this function determines which, yielding one predicted score per peptide
        peptides: [ string ], the variable length (>= 9) peptides
        encodings: { string: [ tensor ] }, from encode_peptides -- for each peptide, the encodings according to different cores
        predictions: [ float ], the predicted scores corresponding to cores
        returns: [ float ], [ int ], the scores corresponding to peptides, the selected core indices"""
        values = [0]*len(peptides); cores = [None]*len(peptides)
        pred_idx = 0
        for (pep_idx,peptide) in enumerate(peptides):
            best_core = None; best_value = None
            for core_idx in range(len(encodings[peptide])):
                if best_core is None or predictions[pred_idx] > best_value:
                    best_core = core_idx; best_value = predictions[pred_idx]
                pred_idx += 1
            values[pep_idx] = best_value
            cores[pep_idx] = best_core
        return (values, cores)   

    def train_fixed_cores(self, models, peptides, encodings, values):
        """Trains the models to predict values from peptides, in the situation where there is exactly one fied core to consider for each peptide
        peptides: [ string ], the variable length (>= 9) peptides
        encodings: { string: [ tensor ] }, from encode_peptides -- for each peptide, the encoding according to the fixed cores (as a list for consistency)
        values: [ float ], the scores corresponding to peptides
        """
        if self.validation_split==0:
            callbacks = []
        else:
            callbacks = [tf.keras.callbacks.EarlyStopping(monitor='val_loss', min_delta=0, patience=self.patience)]
        for model in models:
            print(model.name)
            # preshuffle to get new validation data
            # TODO: similarity reduce it?
            shuffle = np.random.permutation(peptides)
            inputs = np.stack([encoded_core for peptide in shuffle for encoded_core in encodings[peptide]])
            outputs = np.array([values[peptide] for peptide in shuffle])
            history = model.fit(inputs, outputs, epochs=self.max_epochs, batch_size=self.batch_size, validation_split=self.validation_split, callbacks=callbacks, verbose=self.verbose)
            if self.validation_split>0: 
                print('fit', len(history.history['val_loss']), history.history['val_loss'][-1])

    def train_unknown_cores(self, models, peptides, encodings, values):
        """Trains the models to predict values from peptides, considering all possible cores for each peptide and using the strongest score.
        peptides: [ string ], the variable length (>= 9) peptides
        encodings: { string: [ tensor ] }, from encode_peptides -- for each peptide, the encodings according to different cores
        values: [ float ], the scores corresponding to peptides
        """
        for model in models:
            print(model.name)
            # preshuffle to get new validation data
            # TODO: similarity reduce it?
            shuffle = np.random.permutation(peptides)
            inputs_all_cores = np.stack([encoded_core for peptide in shuffle for encoded_core in encodings[peptide]])
            outputs_single_core = np.array([values[peptide] for peptide in shuffle]) # stays the same throughout epochs, even as core changes
            epoch = 0
            best_epoch = None; best_score = None
            while epoch < self.max_epochs:
                # find strongest cores
                (_predicted_values, best_cores) = self.decode_predictions(shuffle, encodings, model.predict(inputs_all_cores))
                # now fix them and train
                inputs_single_core = np.stack([encodings[peptide][best_cores[i]] for (i,peptide) in enumerate(shuffle)])
                history = model.fit(inputs_single_core, outputs_single_core, epochs=1, batch_size=self.batch_size, validation_split=self.validation_split, verbose=self.verbose)
                if self.validation_split==0:
                    print(epoch)
                else:
                    score = history.history['val_loss'][-1]
                    if best_epoch is None or score < best_score:
                        print(epoch, score, 'best in', epoch-(best_epoch if best_epoch is not None else 0))
                        best_epoch = epoch
                        best_score = score
                    elif best_epoch+self.patience == epoch:
                        # TODO revert model?
                        print('converged', epoch, score)
                        break
                    else:
                        print(epoch, score)
                epoch += 1

    def class_of(self, v):
        """1 for binder (proability above the thresh), 0 for non-binder"""
        if v >= self.binder_thresh: return 1
        return 0

    def evaluate_prediction(self, predicted_values, actual_values):
        """Computes some metrics evaluating how well the predicted and actual values match
        returns (auc, f1, pcc)"""
        actual_classes = [self.class_of(v) for v in actual_values]
        predicted_classes = [self.class_of(v) for v in predicted_values]
        auc = sklearn.metrics.roc_auc_score(actual_classes, predicted_values)
        f1 = sklearn.metrics.f1_score(actual_classes, predicted_classes)
        pcc = scipy.stats.pearsonr(actual_values, predicted_values)
        return (auc, f1, pcc[0])

    def evaluate_fold(self, ensemble, peptides, encodings, actual_values):
        """Given an ensemble trained on a fold's training data, evaluate its performance on the testing data
        Augments self.predicted with the predicted peptide -> value entries
        ensemble: [ model ]
        peptides: [ string ], the testing peptides
        encodings: { string: [ tensor ] }, from encode_peptides -- for each peptide, the encodings according to different cores
        values: [ float ], the actual scores corresponding to peptides
        """
        flattened_cores = np.stack([c for p in peptides for c in encodings[p]])
        ensemble_values = []
        nans = 0
        # Note: different models could have different favorite cores, so treat them separately
        # I.e., for each model, find its best value (for whatever core) and then average,
        # rather than averaging predictions over models for each core and then finding the 
        # best value of the averages.
        for model in ensemble:
            (values,cores) = self.decode_predictions(peptides, encodings, model.predict(flattened_cores))
            for (i,v) in enumerate(values):
                if math.isnan(v): 
                    values[i] = 0
                    nans += 1
            ensemble_values.append(values)
        if nans>0: print('!',nans,'nans') # sign of bad training
        for (i,peptide) in enumerate(peptides):
            self.predicted[peptide] = np.average([model_values[i] for model_values in ensemble_values])
        fold_predicted = [self.predicted[peptide] for peptide in peptides]
        fold_actual = [actual_values[peptide] for peptide in peptides]
        print('fold auc %f; f1 %f; pcc %f' % self.evaluate_prediction(fold_predicted, fold_actual))

    def xval_train(self, binding_data, fixed_cores=None):
        """Creates and trains an ensemble of models trained on the cross-validation splits in binding_data.
        Also stores the cross-validated predictions in self.predicted
        binding_data: BindingData
        fixed_cores: { string : int }, the specific core allowed for each peptide, if fixed; else None to consider all"""
        
        with self.graph.as_default():        
            self.build_ensemble(len(binding_data.folds))
        
            # now that encoder and ensemble are in place, make sure all variables are initialized
            self.sess.run(tf.tables_initializer())
            self.sess.run(tf.global_variables_initializer())
            self.sess.run(tf.local_variables_initializer())
        
            encodings = self.peptide_encoder.encode_peptides([peptide for fold in binding_data.folds for peptide in fold], self.sess)
            if fixed_cores is not None:
                fixed_encodings = self.peptide_encoder.encode_peptides([peptide for fold in binding_data.folds for peptide in fold], self.sess, fixed_cores)
            self.predicted = {}
            for (fold, fold_ensemble) in enumerate(self.fold_ensembles):
                training_peptides = [peptide for f in range(len(binding_data.folds)) if f != fold for peptide in binding_data.folds[f]]
                print('*** fold', fold, len(training_peptides), len(binding_data.folds[fold]))
                if fixed_cores is not None:
                    self.train_fixed_cores(fold_ensemble, training_peptides, fixed_encodings, binding_data.values)
                else:
                    self.train_unknown_cores(fold_ensemble, training_peptides, encodings, binding_data.values)
                testing_peptides = binding_data.folds[fold]
                self.evaluate_fold(fold_ensemble, testing_peptides, encodings, binding_data.values)
            flat_predicted = list(self.predicted.values())
            flat_actual = [binding_data.values[p] for p in self.predicted]
            print('***','overall xval: auc %f; f1 %f; pcc %f' % self.evaluate_prediction(flat_predicted, flat_actual))
        
            # package it all up, from peptide to prediction, averaging over the ensemble 
            self.peptide = self.peptide_encoder.peptide
            pep_in = tf.keras.Input(tensor=self.peptide_encoder.encoding)
            cat = tf.keras.layers.concatenate([model(pep_in) for model in self.ensemble])
            self.prediction = tf.reduce_mean(cat, axis=1, name='prediction')

    def save_model(self, base):
        # https://www.tensorflow.org/tfx/serving/serving_basic
        # TODO: this saves it for serving, which maybe isn't the standard case now, so could simplify
        # but then again it doesn't much hurt to have this available (tried playing with it -- pretty cool to have an epitope server)

        # TODO: note, currently not saving out expt_len as an input, since expected use case is to feed a 15mer
        
        with self.graph.as_default(): 
            tf.keras.backend.set_learning_phase(0)  # Ignore dropout at inference
    
            exporter = tf.saved_model.builder.SavedModelBuilder(base)
            sig = tf.saved_model.signature_def_utils.build_signature_def(
                inputs={'peptide':tf.saved_model.utils.build_tensor_info(self.peptide)},
                outputs={'prediction': tf.saved_model.utils.build_tensor_info(self.prediction)},
                method_name=tf.saved_model.signature_constants.PREDICT_METHOD_NAME
                )
            exporter.add_meta_graph_and_variables(
                self.sess,[tf.saved_model.tag_constants.SERVING],
                signature_def_map={tf.saved_model.signature_constants.DEFAULT_SERVING_SIGNATURE_DEF_KEY: sig},
                main_op=tf.tables_initializer(name='main_init')
                )
            exporter.save()

# ============================================================================
# using models
            
# TODO: background distributions; percentile ranks

class YAEPIIAlleleModel (object):
    """Model for a single allele; allows predicting binding affinity for a 15-mer peptide"""
    
    def __init__(self, allele, graph, sess, peptide, prediction):
        self.allele = allele
        self.graph = graph
        self.sess = sess
        self.peptide = peptide
        self.prediction = prediction

    @classmethod
    def load_saved(cls, allele, models_dir):
        """The original model saved by YAEPIITrainer.save_model"""
        graph = tf.Graph()
        sess = tf.Session(graph=graph)
        with graph.as_default():
            model = tf.saved_model.loader.load(sess, [tf.saved_model.tag_constants.SERVING], models_dir+'/'+allele)
            peptide = graph.get_tensor_by_name(model.signature_def['serving_default'].inputs['peptide'].name)
            prediction = graph.get_tensor_by_name(model.signature_def['serving_default'].outputs['prediction'].name)
        return cls(allele, graph, sess, peptide, prediction)
    
    @classmethod
    def load_frozen(cls, allele, models_dir):
        """The model frozen from the saved model."""
        from tensorflow.python.platform import gfile 
        graph = tf.Graph()
        sess = tf.Session(graph=graph)
        with graph.as_default():
            raw = gfile.FastGFile(models_dir+'/'+allele+'/frozen.pb','rb').read()   
            graph_def = tf.GraphDef()
            graph_def.ParseFromString(raw)
            tf.import_graph_def(graph_def, name='yaep')
            init = graph.get_operation_by_name('yaep/main_init')
            sess.run(init)
            peptide = graph.get_tensor_by_name('yaep/peptide:0')
            prediction = graph.get_operation_by_name('yaep/prediction').outputs[0]
            return cls(allele, graph, sess, peptide, prediction)

    def score_peptide(self, peptide):
        res = self.sess.run(self.prediction, feed_dict={self.peptide:[peptide]})
        return res[0]

class YAEPII (EpitopePredictor):
    """Model(s) for one or more allele;s allows predicting binding affinity for a 15-mer peptide and thresholding as hit/not for each allele"""

    allele_sets = { # TODO: read from a spec file in model_base? for now, just copied from NetMHCII
        'test':
            ['DRB1_0101'],
        
        # Greenbaum J, Sidney J, Chung J, Brander C, Peters B, Sette A.
        # Functional classification of class II human leukocyte antigen (HLA) molecules reveals seven different supertypes and a surprising degree of repertoire sharing across supertypes
        # Immunogenetics. 2011 Jun;63(6):325-35
        # https://www.ncbi.nlm.nih.gov/pubmed/21305276
        'greenbaum11':
            ['DRB1_0101', 'DRB1_0301', 'DRB1_0401', 'DRB1_0405', 'DRB1_0701', 'DRB1_0802', 'DRB1_0901', 'DRB1_1101', 'DRB1_1201', 'DRB1_1302', 'DRB1_1501',
             'DRB3_0101', 'DRB3_0202', 'DRB4_0101', 'DRB5_0101',
             'HLA-DPA10103-DPB10201', 'HLA-DPA10103-DPB10401', 'HLA-DPA10201-DPB10101', 'HLA-DPA10201-DPB10501', 'HLA-DPA10201-DPB11401', 'HLA-DPA10301-DPB10402',
             'HLA-DQA10101-DQB10501', 'HLA-DQA10102-DQB10602', 'HLA-DQA10301-DQB10302', 'HLA-DQA10401-DQB10402', 'HLA-DQA10501-DQB10201', 'HLA-DQA10501-DQB10301'],

        # Paul S, Lindestam Arlehamn CS, Scriba TJ, Dillon MB, Oseroff C, Hinz D, McKinney DM, Carrasco Pro S, Sidney J, Peters B, Sette A.
        # Development and validation of a broad scheme for prediction of HLA class II restricted T cell epitopes
        # J Immunol Methods. 2015 Jul;422:28-34
        # https://www.ncbi.nlm.nih.gov/pubmed/25862607
        'paul15':
            ['DRB1_0301', 'DRB1_0701', 'DRB1_1501', 'DRB3_0101', 'DRB3_0202', 'DRB4_0101', 'DRB5_0101'],
        }

    @staticmethod
    def find_models_dir():
        # TODO: look in rosetta database, etc.
         if os.getenv('YAEPII') is not None: return os.getenv('YAEPII')
         return 'models'
        
    def __init__(self, models, models_dir=None, prob_thresh=None):
        super().__init__('yaepii', alleles=[m.allele for m in models], peptide_length=15, overhang=3)
        self.models = models
        self.models_dir = models_dir if models_dir is not None else YAEPII.find_models_dir()
        self.prob_thresh = prob_thresh if prob_thresh is not None else BindingData.nM_to_prob(500)

    def set_alleles(self, alleles):
        # TODO: make sure supported
        self.models = [YAEPIIAlleleModel.load_frozen(allele, self.models_dir) for allele in alleles]
        self.alleles = [m.allele for m in self.models]

    @classmethod
    def load_saved(cls, alleles, models_dir=None, prob_thresh=None):
        if models_dir == None: models_dir = YAEPII.find_models_dir()
        models = [YAEPIIAlleleModel.load_saved(allele, models_dir) for allele in alleles]
        return cls(models, models_dir=models_dir, prob_thresh=prob_thresh)

    @classmethod
    def load_frozen(cls, alleles, models_dir=None, prob_thresh=None):
        if models_dir == None: models_dir = YAEPII.find_models_dir()
        models = [YAEPIIAlleleModel.load_frozen(allele, models_dir) for allele in alleles]
        return cls(models, models_dir=models_dir, prob_thresh=prob_thresh)

    def score_peptide(self, pep):
        details = [m.score_peptide(pep) for m in self.models]
        return EpitopeScore(sum(1 for s in details if s>self.prob_thresh), details)

    # TODO: score_peptides in batch -- feed whole list
