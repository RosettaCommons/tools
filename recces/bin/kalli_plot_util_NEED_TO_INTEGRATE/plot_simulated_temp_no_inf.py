#!/Users/kkappel/anaconda/bin/python

import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib import rc
import data_k as data


def main(args):
	native = data.SingleHistSimulation('./non_native_a20_run1/', uucg_torsion_volume(20,180))
	#non_native = data.SingleHistSimulation('./non_native_a15_run1/', uucg_torsion_volume(15,180))
	colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', '0.75']
	temps = ['0.8', '1', '1.4', '1.8', '3', '7', '30']
	scores = native.hist_scores
	native_hists = native.histograms
	#non_native_hists = non_native.histograms
	native_kTs = native.kt_list
	#non_native_kTs = non_native.kt_list
	native_samples = []
	#non_native_samples = []
	native_total_samples = 0
	#non_native_total_samples = 0
	for i in xrange(0,len(native_hists)): 
		native_total_samples+=native_hists[i].sum()
	#for i in xrange(0,len(non_native_hists)): 
	#	non_native_total_samples+=non_native_hists[i].sum()
	for i in [0,1,2,3,5,6,7]:
	#for i in range(0,len(temps)):
		native_normhist = native_hists[i]/float(native_hists[i].sum())
		native_samples.append(native_hists[i].sum()/float(native_total_samples))
		plt.bar(scores,native_normhist,color=colors[i], width=0.2, linewidth=0, alpha=0.75)
	plt.axis([-50,120,0,0.05])
	plt.xlabel('Score')
	plt.ylabel('Probability')
	plt.legend(["T %.1f, %1.2f" %(native_kTs[0], native_samples[0]), "T %.1f, %1.2f" %(native_kTs[1], native_samples[1]), 
		"T %.1f, %1.2f" %(native_kTs[2], native_samples[2]), "T %.1f, %1.2f" %(native_kTs[3], native_samples[3]), 
		"T %.1f, %1.2f" %(native_kTs[5], native_samples[4]), "T %.1f, %1.2f" %(native_kTs[6], native_samples[5]), 
		"T %.1f, %1.2f" %(native_kTs[7], native_samples[6])], fontsize=8, loc=2)
	plt.savefig("non_native_a20_st_bar.pdf")
	plt.show()

def uucg_torsion_volume(bb_range,chi_range):
        return (np.radians(bb_range)*2)**25*(np.radians(chi_range)*2)**4

if __name__ == '__main__':
	parser=argparse.ArgumentParser(description="Plot histograms from thermal sampler")
	parser.add_argument('-o', '--outfile', help='Name of output file')
	args = parser.parse_args()
	main(args)
