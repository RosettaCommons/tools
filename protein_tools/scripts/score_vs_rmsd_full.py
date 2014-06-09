#!/usr/bin/env python2.7
import sys
import amino_acids
from Bio.Align.Applications import ClustalwCommandline as cw
import warnings
from multiprocessing import Pool
import argparse
import os
from Bio import AlignIO
import rosettaScore_beta as rsb
from Bio.PDB import PDBExceptions
from rosettautil.protein import util, pdbStat
from rosettautil import resfile
import glob

clustalw_exe = r"/sb/apps/clustalw/Linux2/i686/v1.83/clustalw"

warnings.simplefilter('ignore', PDBExceptions.PDBConstructionWarning)

if sys.version_info < (2, 7):
    raise Exception("You must use python2.7 to run this")


# usage
usage = "%prog [options] -n native.pdb *pdbs"
parser = argparse.ArgumentParser(prog="score_vs_rmsd_full.py", formatter_class=argparse.RawTextHelpFormatter,
                                 description="\n\nA rosetta score vs rmsd file encompassing a full analysis.\n\n\
                                 -This recapitulates the most accurate representation of score vs rmsd by doing the following:\n\
                                 1. First does an alignment of sequences to find out amino acid correspondance\n\
                                 2. Iterates through atoms attempting to match them up. If they are not identical types, \n\
                                 \tthen they will be rejected in the calculation\n\
                                 3. (optional) An RMSD will be calculated all after all the atoms are aligned that were considered from steps 1 and 2.\n\
                                 4. An optional resfile can be passed that will only calculate the RMSD of those residues defined.\n\
                                 5. Scores will be calculated for each model for all available scoring terms.\n\
                                 6. (optional) If a resfile is passed, scores for just that part of the model will be calculated for all scoring terms\n\
                                 \n\tTo determine which atoms are rejected please use debug mode")

parser.add_argument("pdbs", metavar="*pdb", nargs=argparse.REMAINDER, help="The pdb files that you want superimposed")
parser.add_argument(
    "--multi", "-m", dest="multi", help="run using multiple processors", action='store_true', default=False)
parser.add_argument("--native", "-n", dest="native", help="native structure")
parser.add_argument("--out", "-o", dest="table", help="the name of the prefix of the  output file", default="score_vs_rmsd")
parser.add_argument("--res", "-r", dest="residues",
                    help="the res file to align by too. If not specified, it will align the whole molecule", default="")
parser.add_argument(
    "--tag", "-g", dest="tag", help="put a tag in your output file. This is useful if you want to combine a bunch of these files and use them in mysql or plot function in R", default="")
parser.add_argument("--debug", "-d", dest="debug", action="store_true", default=False,
                    help="will output a file called debug_atoms.txt to tell you exatctly what was calculated in the RMSD")
args = parser.parse_args()

if len(args.pdbs) < 1:
    parser.error("specify at least 1 protein to compare to native")

def get_fasta_from_pdb(pdb):
    holder = []
    for i in pdb.get_residues():
        try:
            holder.append((i,amino_acids.longer_names[i.get_resname()]))
        except KeyError:
            if args.debug:
                print "Cant find short name for {0}, putting \"X\" instead".format(i)
                holder.append((i,"X"))
    return holder

def get_pairwise_alignment(native,target,target_name):
    alignment_dict = {}
    fasta1 = get_fasta_from_pdb(native)
    fasta2 = get_fasta_from_pdb(target)
    file = "clustal_tmp_in_{0}.fasta".format(target_name)
    with open(file,'w') as f:
        f.write(">native\n")
        for letter1 in fasta1:
            f.write(letter1[1])
        f.write("\n>target\n")
        for letter2 in fasta2:
            f.write(letter2[1])
    cline = cw(clustalw_exe,infile=file)
    assert os.path.isfile(clustalw_exe), "Clustal W binary is missing"
    stdout,stderr = cline()
    alignment = AlignIO.read(file.split('.fasta')[0]+'.aln','clustal')
    alignment_dict['native'] = str(alignment[0,:].seq)
    alignment_dict['target'] = str(alignment[1,:].seq)
    fil = glob.glob("clustal_tmp*")
    #for f in fil:
    #    os.remove(f)
    return alignment_dict



def main(args):
    if args.multi:
        try:
            p = Pool()
            tag_list = p.map(find_rmsd, args.pdbs)
        except KeyboardInterrupt:
            p.terminate()
            quit()
    else:
        tag_list = []
        for i in args.pdbs:
        # stats is a dictionary
            stats = find_rmsd(i)
            tag_list.append(stats)
    make_table(tag_list, args.table, args.tag)
    temp_files = glob.glob('clustal_tmp*')
    for i in temp_files:
        os.remove(i)



def make_table(name_score_rmsd, out_name, tag_label):
    additional_headers = name_score_rmsd[0].values()[0]["pose_score"].keys()
    header = ["MODEL", "CA_RMSD", "BB_RMSD", "ALL_ATOM_RMSD"] + additional_headers + ["\n"]
    with open(out_name + "_align_all_model.tsv", 'w') as f:
        f.write("\t".join(header))
        for entity in name_score_rmsd:
            for model_entity in entity:
                model = model_entity
                all_rmsds = entity[model]["rms_all"]
                pose_scores = entity[model]["pose_score"]
                f.write(model + "\t")
                f.write("{0:.2f}\t".format(all_rmsds["ca_rmsd"]))
                f.write("{0:.2f}\t".format(all_rmsds["bb_rmsd"]))
                f.write("{0:.2f}\t".format(all_rmsds["all_rmsd"]))
                for scores in pose_scores:
                    f.write("{0:.2f}\t".format(pose_scores[scores]))
                f.write("\n")
    if args.residues:
        additional_headers = name_score_rmsd[0].values()[0]["residue_scores"].keys()
        header = ["MODEL", "CA_RMSD", "BB_RMSD", "ALL_ATOM_RMSD"] + additional_headers + ["\n"]
        with open(out_name + "_align_by_residue.tsv", 'w') as f:
            f.write("\t".join(header))
            for entity in name_score_rmsd:
                for model_entity in entity:
                    model = model_entity
                    residue_rmsds = entity[model]["rms_res"]
                    residue_scores = entity[model]["residue_scores"]
                    f.write(model + "\t")
                    f.write("{0:.2f}\t".format(residue_rmsds["ca_rmsd"]))
                    f.write("{0:.2f}\t".format(residue_rmsds["bb_rmsd"]))
                    f.write("{0:.2f}\t".format(residue_rmsds["all_rmsd"]))
                    for scores in residue_scores:
                        f.write("{0:.2f}\t".format(residue_scores[scores]))
                    f.write("\n")


def find_rmsd(model):
    if args.debug:
        print "processing %s" %(model)
    name_score_rmsd = {}
    native = args.native
    rms_all = []
    rms_res = []
    residue_scores = {}
    pwa = get_pairwise_alignment(util.load_pdb(args.native), util.load_pdb(model),os.path.basename(model))
    # will return three rmsds aligned by everything in decoy
    # give it an empty residue file and it will autmatically calculate all rmsds
    rms_all = pdbStat.calculate_all_superpositions(util.load_pdb(native), util.load_pdb(model), [], args.debug,pwa)
    score_table = rsb.ScoreTable(model)
    # dictionary with all the pose scores
    pose_score = score_table.get_pose_all_scores()
    if args.residues:
        # dictionary with all the residue scores
        residues = resfile.Residue_File(args.residues).get_designed_entities()
        rms_res = pdbStat.calculate_all_superpositions(
            util.load_pdb(native), util.load_pdb(model), residues, args.debug,pwa)
        for i in residues:
            score_table_per_residue = score_table.get_all_score_terms(chain=i[0], pdbres=int(i[1]),)
            for i in score_table_per_residue.iterkeys():
                try:
                    residue_scores[i] += score_table_per_residue[i]
                except KeyError:
                    residue_scores[i] = score_table_per_residue[i]
    name_score_rmsd[model] = {"rms_all": rms_all, "rms_res": rms_res,
                              "pose_score": pose_score, "residue_scores": residue_scores}
    return name_score_rmsd


if __name__ == "__main__":
    main(args)
