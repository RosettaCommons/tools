#!/Users/kkappel/anaconda/bin/python

import numpy as np
import argparse
import matplotlib.pyplot as plt
import data_k as data


def main(args):
	native = data.SingleHistSimulation('%s' %(args.datadir), tandemga_torsion_volume(float(args.bb_ang), float(args.chi_ang))) 
	#non_native = data.SingleHistSimulation('./non_native_a15_run1/', uucg_torsion_volume(15,180))
	colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', '0.75', 'b', 'g', 'r', 'c', 'm', 'y']
	#colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', '0.75']
	temps = ['0.8', '1', '1.4', '1.8', '3', '7', '30', '75', '150', '300', '600', '1200', 'Inf']
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
	for i in range(0,len(native_hists)):
		native_normhist = native_hists[i]/float(native_hists[i].sum())
		native_samples.append(native_hists[i].sum()/float(native_total_samples))
		plt.bar(scores,native_normhist,color=colors[i], width=10, linewidth=0, alpha=0.75)
	plt.axis([-100,10000,0,0.007])
	#plt.axis([-100,10000,0,0.025])
	plt.xlabel('Score')
	plt.ylabel('Probability')
	#plt.legend(["T %.1f, %1.2f" %(native_kTs[0], native_samples[0]), "T %.1f, %1.2f" %(native_kTs[1], native_samples[1]), 
	plt.savefig("non_native_infT_a%s_c%s_st_bar.pdf" %(args.bb_ang, args.chi_ang))
	plt.show()

def uucg_torsion_volume(bb_range,chi_range):
        return (np.radians(bb_range)*2)**25*(np.radians(chi_range)*2)**4

def tandemga_torsion_volume(bb_range,chi_range):
        return (np.radians(bb_range)*2)**50*(np.radians(chi_range)*2)**8

if __name__ == '__main__':
	parser=argparse.ArgumentParser(description="Plot histograms from thermal sampler")
	parser.add_argument('-d', '--datadir', required=True, help='Directory that contains the data')
	parser.add_argument('-b', '--bb_ang', required=True, help='Backbone angle range')
	parser.add_argument('-c', '--chi_ang', required=True, help='Chi angle range')
	parser.add_argument('-o', '--outfile', help='Name of output file')
	args = parser.parse_args()
	main(args)
