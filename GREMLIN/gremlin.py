#!/usr/bin/env python3

# ------------------------------------------------------------------
# "THE BEERWARE LICENSE" (Revision 42)
# ------------------------------------------------------------------
# <so@g.harvard.edu> and <pkk382@g.harvard.edu> wrote this code.
# As long as you retain this notice, you can do whatever you want
# with this stuff. If we meet someday, and you think this stuff
# is worth it, you can buy us a beer in return.
# --Sergey Ovchinnikov and Peter Koo
# ------------------------------------------------------------------
# The original MATLAB code for GREMLIN was written by Hetu Kamisetty
# ------------------------------------------------------------------

# This is a modified version of GREMLIN_TF v2.1.BETA.1 (81eaa6f)
# Modifications were done in the get_mtx function to return the raw tensor:
# Authors: Moritz Ertelt (moritz.ertelt@stud.uni-regensburg.de), Samuel Schmitz (Samuel.Schmitz@vanderbilt.edu)

import numpy as np
import tensorflow as tf
import matplotlib.pylab as plt
from scipy import stats
from scipy.spatial.distance import pdist,squareform
import pandas as pd
import pdb

######## GREMLIN CODE ###############
#####################################

# note: if you are modifying the alphabet
# make sure last character is "-" (gap)
################
alphabet = "ARNDCQEGHILKMFPSTWYV-"
states = len(alphabet)

# map amino acids to integers (A->0, R->1, etc)
a2n = dict((a,n) for n,a in enumerate(alphabet))
aa2int = lambda x: a2n.get(x,a2n['-'])

# from fasta
def parse_fasta(filename):
  '''function to parse fasta'''
  header = []
  sequence = []
  lines = open(filename, "r")
  for line in lines:
    line = line.rstrip()
    if line[0] == ">":
      header.append(line[1:])
      sequence.append([])
    else:
      sequence[-1].append(line)
  lines.close()
  sequence = [''.join(seq) for seq in sequence]
  return np.array(header), np.array(sequence)

def filt_gaps(msa, gap_cutoff=0.5):
  '''filters alignment to remove gappy positions'''
  frac_gaps = np.mean((msa == states-1).astype(np.float),0)
  non_gaps = np.where(frac_gaps < gap_cutoff)[0]
  return msa[:,non_gaps], non_gaps

def get_eff(msa, eff_cutoff=0.8):
  '''compute effective weight for each sequence'''  
  msa_sm = 1.0 - squareform(pdist(msa,"hamming"))
  msa_w = (msa_sm >= eff_cutoff).astype(np.float)
  msa_w = 1.0/np.sum(msa_w,-1)
  return msa_w

def str2int(x):
  '''convert a list of strings into list of integers'''
  # Example: ["ACD","EFG"] -> [[0,4,3], [6,13,7]]
  if x.dtype.type is np.str_:
    if x.ndim == 0: return np.array([aa2int(aa) for aa in x])
    else: return np.array([[aa2int(aa) for aa in seq] for seq in x])
  else:
    return x
  
def split_train_test(seqs, frac_test=0.1):
  # shuffle data
  x = np.copy(seqs)
  np.random.shuffle(x[1:])

  # fraction of data used for testing
  split = int(len(x) * (1.0-frac_test))

  # split training/test datasets
  return x[:split], x[split:]

def mk_msa(seqs, gap_cutoff=0.5, eff_cutoff=0.8):
  '''converts list of sequences to MSA (Multiple Sequence Alignment)'''
  # =============================================================================
  # The function takes a list of sequences (strings) and returns a (dict)ionary
  # containing the following:
  # =============================================================================
  # BEFORE GAP REMOVAL
  # -----------------------------------------------------------------------------
  # msa_ori   msa
  # ncol_ori  number of columns
  # -----------------------------------------------------------------------------
  # AFTER GAP REMOVAL
  # little complicated, as we need to keep track of which positions were removed.
  # -----------------------------------------------------------------------------
  # msa       msa
  # ncol      number of columns
  # v_idx     index of positions kept
  # -----------------------------------------------------------------------------
  # weights   weight for each sequence (based on sequence identity)
  # nrow      number of rows (sequences)
  # neff      number of effective sequences sum(weights)
  # =============================================================================
  
  msa_ori = str2int(seqs)

  # remove positions with more than > 50% gaps
  msa, v_idx = filt_gaps(msa_ori, gap_cutoff)
  
  # compute effective weight for each sequence
  msa_weights = get_eff(msa, eff_cutoff)
    
  return {"msa_ori":msa_ori,
          "msa":msa,
          "weights":msa_weights,
          "neff":np.sum(msa_weights),
          "v_idx":v_idx,
          "nrow":msa.shape[0],
          "ncol":msa.shape[1],
          "ncol_ori":msa_ori.shape[1]}

# optimizer
def opt_adam(loss, name, var_list=None, lr=1.0, b1=0.9, b2=0.999, b_fix=False):
  # adam optimizer
  # Note: this is a modified version of adam optimizer. More specifically, we replace "vt"
  # with sum(g*g) instead of (g*g). Furthmore, we find that disabling the bias correction
  # (b_fix=False) speeds up convergence for our case.
  
  if var_list is None: var_list = tf.trainable_variables() 
  gradients = tf.gradients(loss,var_list)
  if b_fix: t = tf.Variable(0.0,"t")
  opt = []
  for n,(x,g) in enumerate(zip(var_list,gradients)):
    if g is not None:
      ini = dict(initializer=tf.zeros_initializer,trainable=False)
      mt = tf.get_variable(name+"_mt_"+str(n),shape=list(x.shape), **ini)
      vt = tf.get_variable(name+"_vt_"+str(n),shape=[], **ini)
      
      mt_tmp = b1*mt+(1-b1)*g
      vt_tmp = b2*vt+(1-b2)*tf.reduce_sum(tf.square(g))
      lr_tmp = lr/(tf.sqrt(vt_tmp) + 1e-8)

      if b_fix: lr_tmp = lr_tmp * tf.sqrt(1-tf.pow(b2,t))/(1-tf.pow(b1,t))

      opt.append(x.assign_add(-lr_tmp * mt_tmp))
      opt.append(vt.assign(vt_tmp))
      opt.append(mt.assign(mt_tmp))
        
  if b_fix: opt.append(t.assign_add(1.0))
  return(tf.group(opt))


def GREMLIN(msa,
            opt_iter=100,
            opt_rate=1.0,
            batch_size=None,
            lam_v=0.01,
            lam_w=0.01,
            scale_lam_w=True,
            v=None,
            w=None,
            ignore_gap=True):
  
  '''fit params of MRF (Markov Random Field) given MSA (multiple sequence alignment)'''
  # ==========================================================================
  # this function takes a MSA (dict)ionary, from mk_msa() and returns a MRF
  # (dict)ionary containing the following:
  # ==========================================================================
  # len       full length
  # v_idx     index of positions (mapping back to full length)
  # v         2-body term
  # w         2-body term
  # ==========================================================================
  # WARNING: The mrf is over the msa after gap removal. "v_idx" and "len" are
  # important for mapping the MRF back to the original MSA.
  # ==========================================================================
  
  ########################################
  # SETUP COMPUTE GRAPH
  ########################################
  # reset tensorflow graph
  tf.reset_default_graph()
  
  # length of sequence
  ncol = msa["ncol"] 
  
  # input msa (multiple sequence alignment) 
  MSA = tf.placeholder(tf.int32,shape=(None,ncol),name="msa")
  
  # input msa weights
  MSA_weights = tf.placeholder(tf.float32, shape=(None,), name="msa_weights")
  
  # one-hot encode msa
  OH_MSA = tf.one_hot(MSA,states)
  
  if ignore_gap:
    ncat = states - 1
    NO_GAP = 1.0 - OH_MSA[...,-1] 
    OH_MSA = OH_MSA[...,:ncat]
    
  else:
    ncat = states
  
  ########################################
  # V: 1-body-term of the MRF
  ########################################
  V = tf.get_variable(name="V",
                          shape=[ncol,ncat],
                          initializer=tf.zeros_initializer)
  
  ########################################
  # W: 2-body-term of the MRF
  ########################################
  W_tmp = tf.get_variable(name="W",
                          shape=[ncol,ncat,ncol,ncat],
                          initializer=tf.zeros_initializer)  
  
  # symmetrize W
  W = W_tmp + tf.transpose(W_tmp,[2,3,0,1])
  
  # set diagonal to zero
  W = W * (1-np.eye(ncol))[:,None,:,None]

  ########################################
  # Pseudo-Log-Likelihood
  ########################################
  # V + W
  VW = V + tf.tensordot(OH_MSA,W,2)
  
  # hamiltonian
  H = tf.reduce_sum(OH_MSA*VW,-1)
  
  # local Z (parition function)
  Z = tf.reduce_logsumexp(VW,-1)

  PLL = H - Z
  if ignore_gap:
    PLL = PLL * NO_GAP  

  PLL = tf.reduce_sum(PLL,-1)  
  PLL = tf.reduce_sum(MSA_weights * PLL)/tf.reduce_sum(MSA_weights)

  ########################################
  # Regularization
  ########################################
  L2 = lambda x: tf.reduce_sum(tf.square(x))
  L2_V = lam_v * L2(V)
  L2_W = lam_w * L2(W) * 0.5
  
  if scale_lam_w:
    L2_W = L2_W * (ncol-1) * (states-1)
  
  ########################################
  # Loss Function
  ########################################
  # loss function to minimize
  loss = -PLL + (L2_V + L2_W) / msa["neff"]
  
  # optimizer
  opt = opt_adam(loss,"adam",lr=opt_rate)
  
  ########################################
  # Input Generator
  ########################################
  all_idx = np.arange(msa["nrow"])
  def feed(feed_all=False):
    if batch_size is None or feed_all:
      return {MSA:msa["msa"], MSA_weights:msa["weights"]}
    else:
      batch_idx = np.random.choice(all_idx,size=batch_size)
      return {MSA:msa["msa"][batch_idx], MSA_weights:msa["weights"][batch_idx]}
  
  ########################################
  # OPTIMIZE
  ########################################
  with tf.Session() as sess:
    
    # initialize variables V and W
    sess.run(tf.global_variables_initializer())

    # initialize V
    if v is None:
      oh_msa = np.eye(states)[msa["msa"]]
      if ignore_gap: oh_msa = oh_msa[...,:-1]
      
      pseudo_count = 0.01 * np.log(msa["neff"])
      f_v = np.einsum("nla,n->la",oh_msa,msa["weights"])
      V_ini = np.log(f_v + pseudo_count)
      if lam_v > 0:
        V_ini = V_ini - np.mean(V_ini,axis=-1,keepdims=True)
      sess.run(V.assign(V_ini))
      
    else:
      sess.run(V.assign(v))

    # initialize W
    if w is not None:
      sess.run(W_tmp.assign(w * 0.5))
      
    # compute loss across all data
    get_loss = lambda: np.round(sess.run(loss,feed(True)) * msa["neff"],2)

    print("starting",get_loss())      
    for i in range(opt_iter):
      sess.run(opt,feed())  
      if (i+1) % int(opt_iter/10) == 0:
        print("iter",(i+1),get_loss())
    
    # save the V and W parameters of the MRF
    V_ = sess.run(V)
    W_ = sess.run(W)
    
  ########################################
  # return MRF
  ########################################
  no_gap_states = states - 1
  mrf = {"v": V_[:,:no_gap_states],
         "w": W_[:,:no_gap_states,:,:no_gap_states],
         "v_idx": msa["v_idx"],
         "len": msa["ncol_ori"]}
  
  return mrf

########################## INPUT MSA #########################

###############
## FUNCTIONS
###############
def normalize(x):
  x = stats.boxcox(x - np.amin(x) + 1.0)[0]
  x_mean = np.mean(x)
  x_std = np.std(x)
  return((x-x_mean)/x_std)

###############
## FUNCTIONS
###############
def get_mtx(mrf):
    '''convert MRF (Markov Random Field) to MTX (Matrix or Contact-map)
      new version: uses APC on the whole tensor to better capture the difference in score
      between different amino acid combinations'''
  #--------------------------- Original Gremlin Function ---------------------------------
  # raw (l2norm of each 20x20 matrix)
    raw_sq = np.sqrt(np.sum(np.square(mrf["w"]),(1,3)))
    raw = squareform(raw_sq, checks=False)
    
  
  # apc (average product correction)
    ap_sq = np.sum(raw_sq,0,keepdims=True) * np.sum(raw_sq,1,keepdims=True)/np.sum(raw_sq)
    apc = squareform(raw_sq - ap_sq, checks=False)
    i, j = np.triu_indices_from(raw_sq,k=1)
    mtx = {
         "i": mrf["v_idx"][i],
         "j": mrf["v_idx"][j],
         "raw": raw,
         "apc": apc,
         "zscore": normalize(apc),
         "len": mrf["len"]}  
 #----------------------------Custom Function to get AA combinations--------------------
    # Just accessing the raw tensor 
    raw_all = np.zeros((20,20,len(raw)))

    for i in range(20):
        for j in range(20):
            raw_all_singleAA = squareform(mrf['w'][:,i,:,j], checks=False)
            raw_all[i,j,:] = raw_all_singleAA
   

   # Different implementations of Average Product Correction (APC)
    raw_sq_all = np.absolute(mrf["w"])

    # calculating the APC correction term for every Value in every 20x20 matrix
    ap_sq_all = np.sum(raw_sq_all, (2,3), keepdims=True) * np.sum(raw_sq_all, (0,1), keepdims=True)/np.sum(raw_sq_all)        

    # Correcting the values in our tensor
    resLength = mrf['w'].shape[0]
    tensorAPC = np.zeros((resLength, 20, resLength, 20))
    for i in range(20):
        for j in range(20):
            apcTensor = raw_sq_all[:,i,:,j] - ap_sq_all[:,i,:,j]
            tensorAPC[:,i,:,j] = apcTensor

    flatTensorAPC = tensorAPC.flatten()
    
    # setting values that are negative after correction to 0
    for index, value in enumerate(flatTensorAPC):
        if value < 0:
            flatTensorAPC[index] = 0

    # we lost the sign with our L2 normalization so here we reassign it
    for index, value in enumerate(mrf['w'].flatten()):   
        if value < 0:
            flatTensorAPC[index] = flatTensorAPC[index] * -1

    correctedTensorAPC = flatTensorAPC.reshape(resLength,20,resLength,20)
    
    # converting from dim[i,aa1,j,aa2] to [i,j,aa12]
    arrayAPC = np.zeros((20, 20, len(raw)))
    for i in range(20):
        for j in range(20):
            tensorAPCsingleAA = squareform(correctedTensorAPC[:,i,:,j], checks=False)
            arrayAPC[i,j,:] = tensorAPCsingleAA

   # Different APC correction where we set negative values to 0 before APC correction

    mrf_flat = mrf['w'].flatten()
    for index, value in enumerate(mrf_flat):
        if value < 0:
            mrf_flat[index] = 0
    
    raw_sq_all0 = np.absolute(mrf_flat.reshape(mrf['w'].shape))
    
    ap_sq_all0 = np.sum(raw_sq_all0, (2,3), keepdims=True) * np.sum(raw_sq_all0, (0,1), keepdims=True)/np.sum(raw_sq_all0)        

    # Correcting the values in our tensor
    resLength = mrf['w'].shape[0]
    tensorAPC0 = np.zeros((resLength, 20, resLength, 20))
    for i in range(20):
        for j in range(20):
            apcTensor0 = raw_sq_all0[:,i,:,j] - ap_sq_all0[:,i,:,j]
            tensorAPC0[:,i,:,j] = apcTensor0

    flatTensorAPC0 = tensorAPC0.flatten()
    
    # setting values that are negative after correction to 0
    for index, value in enumerate(flatTensorAPC0):
        if value < 0:
            flatTensorAPC0[index] = 0

    correctedTensorAPC0 = flatTensorAPC0.reshape(resLength,20,resLength,20)
    
    # converting from dim[i,aa1,j,aa2] to [i,j,aa12]
    arrayAPC0 = np.zeros((20, 20, len(raw)))
    for i in range(20):
        for j in range(20):
            tensorAPCsingleAA0 = squareform(correctedTensorAPC0[:,i,:,j], checks=False)
            arrayAPC0[i,j,:] = tensorAPCsingleAA0

    i, j = np.triu_indices_from(raw_sq,k=1)
    
    mtx_all = {
         "i": np.tile((mrf["v_idx"][i]+1), 400),
         "j": np.tile((mrf["v_idx"][j]+1), 400),
         "raw": raw_all.flatten(),
         "apc": arrayAPC.flatten(),
         "apc0": arrayAPC0.flatten(),
         "len": mrf["len"]} 
    
    return mtx, mtx_all, correctedTensorAPC, correctedTensorAPC0

def plot_mtx(mtx):
  '''plot the mtx'''
  plt.figure(figsize=(15,5))
  for n, key in enumerate(["raw","apc","zscore"]):
    
    # create empty mtx
    m = np.ones((mtx["len"],mtx["len"])) * np.nan
    
    # populate
    m[mtx["i"],mtx["j"]] = mtx[key]
    m[mtx["j"],mtx["i"]] = m[mtx["i"],mtx["j"]]
    
    #plot
    plt.subplot(1,3,n+1)
    plt.title(key)
    if key == "zscore": plt.imshow(m, cmap='Blues', vmin=1, vmax=3)
    else: plt.imshow(m, cmap='Blues')
    plt.grid(False)
    
  plt.savefig("ContactMap.png", dpi=1000)
