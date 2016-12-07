#!/usr/bin/env python
import os
import stat
import sys
import tempfile
import argparse
import shutil
import glob
import math

def is_clustalw_line(line):
    nstring = len(line.strip().split())
    if  nstring < 2 or nstring > 3: return False
    if line.strip().split()[1].find("*") != -1: return False
    if line.strip().split()[1].find(".") != -1: return False
    if line.strip().split()[1].find(":") != -1: return False
    if line.strip().split()[1].find("X") != -1: return False
    return True

class SingleAlignment:
    def __init__(self):
        self.query_sequence = ""

        self.query_tag = ""
        self.target_tag = ""
        self.target_pdb_code = ""
        self.target_pdb_chain = ""

        self.query_start = -1
        self.query_aln_seq = ""

        self.target_start = -1
        self.target_aln_seq = ""

        self.score_line = ""

    def read_grishin(self, lines):
        aln_block = False
        for line in lines:
            if line[:2] == "##":
                aln_block = True
                is_query = True
                buff = line.strip().split()
                self.query_tag = buff[1]
                self.target_tag = buff[2]
                self.target_pdb_code = self.target_tag[:4]
                self.target_pdb_chain = self.target_tag[4]
                continue

            if line[:2] == "--":
                if aln_block:
                    return self
                else:
                    continue

            if line.strip() == "#": continue
            if line.find("scores_from_program") == 0:
                self.score_line = line.strip()
            if not line[0].isdigit(): continue

            if aln_block:
                if is_query:
                    query_start_local, self.query_aln_seq = line.strip().split()
                    self.query_start = int(query_start_local) + 1
                    is_query = False
                else:
                    target_start_local,self.target_aln_seq = line.strip().split()
                    self.target_start = int(target_start_local) + 1
        return self

    def read_fasta(self, lines):
        aln_line_numbers = [iline for iline in range(len(lines)) if lines[iline][0] == '>']
        assert (len(aln_line_numbers) == 2), "Wrong format for FASTA alignment file"

        tag = lines[aln_line_numbers[0]][1:].strip()
        self.target_pdb_code = tag[:4]
        self.target_pdb_chain = tag[4]
        self.target_tag = "%s%s_thread"%(self.target_pdb_code,self.target_pdb_chain)
        self.target_start = 1
        for line in lines[aln_line_numbers[0]+1:aln_line_numbers[1]]:
            self.target_aln_seq += line.strip()

        self.query_tag = lines[aln_line_numbers[1]][1:].strip()
        self.query_start = 1
        for line in lines[aln_line_numbers[1]+1:]:
            self.query_aln_seq += line.strip()

    def read_hhsearch(self, lines):
        aln_block = False
        for line in lines:
            if line[:3] == "No ":
                if not aln_block:
                    aln_block = True
                    rank = int(line.strip().split()[1])
                    continue
                else:
                    aln_block = False
                    return self

            if aln_block:
                if line.strip() == "": continue
                if line[0] == ">":
                    target_pdb = line[1:].strip().split()[0]
                    self.target_pdb_code = target_pdb.split('_')[0]
                    self.target_pdb_chain = target_pdb.split('_')[1]
                    self.target_tag = "%s%s_%d"%(self.target_pdb_code, self.target_pdb_chain, rank+200)
                if line[:2] == "Q " or line[:2] == "T ":
                    if line.find("ss_pred") != -1: continue
                    if line.find("Consensus") != -1: continue
                    if line.find("ss_dssp") != -1: continue
                    buff = line.strip().split()
                    if line[:2] == "Q ":
                        self.query_tag = buff[1]
                        if self.query_start == -1:
                            self.query_start = int(buff[2])
                        self.query_aln_seq += buff[3]
                    if line[:2] == "T ":
                        if self.target_start == -1:
                            self.target_start = int(buff[2])
                        self.target_aln_seq += buff[3]
        return self

    def read_clustalw(self, lines):
        block = []
        for line in lines:
            if line.find("CLUSTAL ") == 0: continue
            if line.strip() != "\n":
                if is_clustalw_line(line):
                    block.append(line)
            else:
                if len(block) != 0:
                    self.query_start = 1
                    self.target_start = 1

                    buff = block[0].strip().split()
                    self.target_pdb_code = buff[0][:4]
                    self.target_pdb_chain = buff[0][4]
                    self.target_tag = "%s%s_thread"%(self.target_pdb_code,self.target_pdb_chain)
                    self.target_aln_seq += buff[1]

                    buff = block[-1].strip().split()
                    self.query_tag = buff[0]
                    self.query_aln_seq += buff[1]
                    block = []

        if len(block) != 0:
            self.query_start = 1
            self.target_start = 1
	    target_line = block[0].strip().split()
	    self.target_pdb_code = target_line[0][:4]
	    self.target_pdb_chain = target_line[0][4]
	    self.target_tag = "%s%s_thread"%(self.target_pdb_code,self.target_pdb_chain)
	    query_line = block[1].strip().split()
	    self.query_tag = query_line[0]

	    # Sequences can span muliple lines
	    for i in range (0,len(block) - 1,2):
		    buff = block[i].strip().split()
		    self.target_aln_seq += buff[1]

		    buff = block[i+1].strip().split()
		    self.query_aln_seq += buff[1]
	    block = []
	    return self

    def grishin_lines(self):
        outlines = []
        outlines.append("## %s %s\n"%(self.query_tag, self.target_tag))
        #outfile.write("## %s %s%s_%d\n"%(self.target_tags, self.target_pdb_and_chain[i_aln][0], self.target_pdb_and_chain[i_aln][1], i_aln+thread))
        outlines.append("#  \n")
        if self.score_line == "":
            outlines.append("scores_from_program: 0\n")
        else:
            outlines.append("%s\n"%self.score_line)
        outlines.append("%d %s\n"%(self.query_start-1, self.query_aln_seq))
        outlines.append("%d %s\n"%(self.target_start-1, self.target_aln_seq))
        outlines.append("--\n")
        return outlines

    def seqence_identity(self):
        assert(len(self.query_aln_seq) == len(self.target_aln_seq))
        seq_id = 0
        query_lenth = 0

        for i in range(len(self.query_aln_seq)):
            if self.query_aln_seq[i] != '-':
                query_lenth += 1
                if self.query_aln_seq[i] == self.target_aln_seq[i]:
                    seq_id += 1
        return float(seq_id)/float(query_lenth)

class Alignment:
    def __init__(self):
        self.alignments = []

    def read_grishin(self, lines):
        aln_block = False
        for line in lines:
            if line[:2] == "##":
                single_aln_lines = []
                aln_block = True
                single_aln_lines.append(line)
                continue

            if line[:2] == "--":
                single_aln_lines.append(line)
                aln = SingleAlignment()
                aln.read_grishin(single_aln_lines)
                self.alignments.append(aln)
                aln_block = False
                continue

            single_aln_lines.append(line)
        return

    def read_hhsearch(self, lines):
        aln_block = False
        for line in lines:
            if line[:3] == "No ":
                if aln_block:
                    aln = SingleAlignment()
                    aln.read_hhsearch(single_aln_lines)
                    self.alignments.append(aln)

                single_aln_lines = []
                single_aln_lines.append(line)
                aln_block = True
                continue

            if aln_block:
                single_aln_lines.append(line)

        aln = SingleAlignment()
        aln.read_hhsearch(single_aln_lines)
        self.alignments.append(aln)
        return

    def read_clustalw(self, lines):
        n_alignments = -1
        block = []
        for line in lines:
            if line.find("CLUSTAL ") == 0: continue
            if is_clustalw_line(line):
                block.append(line)
            else:
                if len(block) != 0:
                    break
        n_alignments = len(block)-1

        for i_aln in range(n_alignments):
            single_aln_lines = []
            block_line_count = 0
            for line in lines:
                if line.find("CLUSTAL ") == 0:
                    single_aln_lines.append(line)
                if line == "\n":
                    single_aln_lines.append(line)
                    block_line_count = 0
                else:
                    block_line_count += 1
                    if block_line_count == i_aln + 1 or block_line_count > n_alignments:
                        single_aln_lines.append(line)
            aln = SingleAlignment()
            aln.read_clustalw(single_aln_lines)
            self.alignments.append(aln)
        return

    def read_fasta(self, lines):
        aln_line_numbers = [iline for iline in range(len(lines)) if lines[iline][0] == '>']
        aln_line_numbers.append(len(lines))
        assert (len(aln_line_numbers) >= 3), "Wrong format for FASTA alignment file"
        n_alignments = len(aln_line_numbers) - 1
        for i_aln in range(1,n_alignments):
            single_aln_lines = []
            # add target
            for iline in range(aln_line_numbers[i_aln], aln_line_numbers[i_aln+1]):
                single_aln_lines.append(lines[iline])
            # add query
            for iline in range(aln_line_numbers[0], aln_line_numbers[1]):
                single_aln_lines.append(lines[iline])
            aln = SingleAlignment()
            aln.read_fasta(single_aln_lines)
            self.alignments.append(aln)
        return

    def read_vienna(self, lines):
        aln_block = False
        is_query = True

        for line in lines:
            if line[0] == ">":
                aln_block = True
                tag = line[1:].strip()
                if is_query:
                    query_tag = tag
                else:
                    target_tag = tag.replace(".pdb","")
                    target_tag += "_%d"%(thread+len(self.query_aln_seq))
                    print target_tag
                continue
            else:
                 # sequence line
                if len(line.strip()) == 0: continue
                alignment_line = line.strip()
                if is_query:
                    query_line = alignment_line
                else:
                    single_aln = SingleAlignment()

                    single_aln.query_tag = query_tag

                    single_aln.target_tag = target_tag
                    single_aln.target_aln_seq = alignment_line
                    single_aln.target_start = 1

                    pdbtag = target_tag
                    single_aln.target_pdb_code = pdbtag[:4]
                    single_aln.target_pdb_chain = pdbtag[4]

                    single_aln.query_aln_seq = query_line
                    single_aln.query_start = 1

                    self.add_alignment(single_aln)
                is_query=False
        return

    def read_modeller(self, lines):
        aln_block = False
        aln_info_defined = False
        is_query = True

        for line in lines:
            # parse first line
            if line[:4] == ">P1;":
                aln_block = True
                tag = line[4:].strip()
                continue

            # parse second line
            if aln_block and not aln_info_defined:
                buff = line.strip().split(':')
                if buff[0] == "sequence":
                    is_query = True
                    query_tag = tag
                else:
                    is_query = False
                    tgt_pdb_code = buff[1]
                    tgt_chain = buff[3]
                    tgt_start_seq = int(buff[2])

                alignment_line = ""
                aln_info_defined= True
                continue

            # parse aligned sequence
            if aln_block and aln_info_defined:
                s = line.strip()
                if s[-1] == "*":
                    alignment_line += s[:-1]
                    if is_query:
                        query_line = alignment_line
                    else:
                        self.target_pdb_and_chain.append((tgt_pdb_code,tgt_chain))
                        self.target_aln_seq.append(alignment_line)
                        self.target_start.append(tgt_start_seq)
                        self.target_tags.append(tag)

                    # end of alignment process
                    aln_block = False
                    aln_info_defined = False
                    continue
                else:
                    alignment_line += s
                    continue

        for i_aln in range(len(self.target_aln_seq)):
            self.query_aln_seq.append(query_line)
            self.query_start.append(1)
        return

    def write_grishin(self, out_fn):
        assert(len(self.alignments) > 0), "Input alignment empty!!"
	outfile = open(out_fn,'w')
        for aln in self.alignments:
            for line in aln.grishin_lines():
                outfile.write(line)
        outfile.close()

    def filter_by_seqence_identity(self,thr):
        for i in range(len(self.alignments)-1, -1, -1):
            if self.alignments[i].seqence_identity() >= thr:
                self.alignments.pop(i)

    def convert(self, aln_in, aln_in_fmt, aln_out, aln_out_fmt):
        lines_in = open(aln_in).readlines()

        if aln_in_fmt == "modeller":
            self.read_modeller(lines_in)
        elif aln_in_fmt == "vie":
            self.read_vienna(lines_in)
        elif aln_in_fmt == "hhsearch":
            self.read_hhsearch(lines_in)
        elif aln_in_fmt == "clustalw":
            self.read_clustalw(lines_in)
        elif aln_in_fmt == "fasta":
            self.read_fasta(lines_in)
        else:
            print "Do not understand input alignment format: %s"%aln_in_fmt
            sys.exit(-1)

        if aln_out_fmt == "grishin":
            self.write_grishin(aln_out)
        else:
            print "Do not understand output alignment format: %s"%aln_out_fmt
            sys.exit(-1)
        return


def write_flags(flag_fn, fasta_fn, xml_fn, silent_fn):
    flag_file=open(flag_fn,'w')
    flag_file.write("# i/o\n")
    flag_file.write("-in:file:fasta %s\n"%fasta_fn)
    flag_file.write("-nstruct 20\n")
    flag_file.write("-parser:protocol %s\n\n"%xml_fn)

    flag_file.write("# relax options\n")
    flag_file.write("-relax:dualspace\n")
    flag_file.write("-relax:jump_move true\n")
    flag_file.write("-default_max_cycles 200\n")
    flag_file.write("-beta_cart\n")
    flag_file.write("-hybridize:stage1_probability 1.0\n")

def write_xml(fn, template_filenames):
    xml_file=open(fn,'w')
    xml_file.write("<ROSETTASCRIPTS>\n")
    xml_file.write("    <TASKOPERATIONS>\n")
    xml_file.write("    </TASKOPERATIONS>\n")
    xml_file.write("    <SCOREFXNS>\n")
    xml_file.write("        <ScoreFunction name=\"stage1\" weights=\"score3\" symmetric=\"0\">\n")
    xml_file.write("            <Reweight scoretype=\"atom_pair_constraint\" weight=\"0.1\"/>\n")
    xml_file.write("        </ScoreFunction>\n")
    xml_file.write("        <ScoreFunction name=\"stage2\" weights=\"score4_smooth_cart\" symmetric=\"0\">\n")
    xml_file.write("            <Reweight scoretype=\"atom_pair_constraint\" weight=\"0.1\"/>\n")
    xml_file.write("        </ScoreFunction>\n")
    xml_file.write("        <ScoreFunction name=\"fullatom\" weights=\"beta_cart\" symmetric=\"0\">\n")
    xml_file.write("            <Reweight scoretype=\"atom_pair_constraint\" weight=\"0.1\"/>\n")
    xml_file.write("        </ScoreFunction>\n")
    xml_file.write("    </SCOREFXNS>\n")
    xml_file.write("    <FILTERS>\n")
    xml_file.write("    </FILTERS>\n")
    xml_file.write("    <MOVERS>\n")
    xml_file.write("        <Hybridize name=\"hybridize\" stage1_scorefxn=\"stage1\" stage2_scorefxn=\"stage2\" fa_scorefxn=\"fullatom\" batch=\"1\" stage1_increase_cycles=\"1.0\" stage2_increase_cycles=\"1.0\">\n")
    for tmpl in template_filenames:
        xml_file.write("            <Template pdb=\"%s\" cst_file=\"AUTO\" weight=\"1.0\" />\n"%(tmpl))
    xml_file.write("        </Hybridize>\n")
    xml_file.write("    </MOVERS>\n")
    xml_file.write("    <APPLY_TO_POSE>\n")
    xml_file.write("    </APPLY_TO_POSE>\n")
    xml_file.write("    <PROTOCOLS>\n")
    xml_file.write("        <Add mover=\"hybridize\"/>\n")
    xml_file.write("    </PROTOCOLS>\n")
    xml_file.write("</ROSETTASCRIPTS>\n")
    xml_file.close()

def download_pdb(pdb_id, dest_dir):
    #print "downloading %s" % ( pdb_id )
    url      = 'http://www.rcsb.org/pdb/files/%s.pdb.gz' % ( pdb_id.upper() )
    dest     = '%s.pdb.gz' % ( pdb_id.lower() )
    wget_cmd = 'curl --silent -L %s -o %s' % ( url, dest )
    if args.verbose: print "Running %s"%wget_cmd
    log_lines = os.popen( wget_cmd ).readlines()
    os.system("gunzip -f %s"%dest)

    if ( os.path.exists(dest[:-3]) ):
        return dest[:-3]
    else:
        print "Error: cannot download PDB %s!"%pdb_id

def read_zipped_pdb_file(fn):
    lines = os.popen( 'zcat %s'%(fn),'r').readlines()
    return lines

def keep_chain(pdblines, chainID):
    outlines = []
    for line in pdblines:
        if line[:6] != "ATOM  ": continue
        if line[21] in chainID:
            outlines.append(line)
    return outlines

def save_pdb(lines, out_fn):
    out = open(out_fn, 'w')
    for line in lines:
        out.write(line)
    out.close()
    return out_fn

def download_templates(aln_file):
    lines = open(aln_file).readlines()
    tmpl_ids = []
    tmpl_fullnames = []
    for line in lines:
        if line[:2] == "##":
            buff = line.strip().split()
            tmpl_ids.append(buff[2][:5])
            tmpl_fullnames.append("%s.pdb"%buff[2])
    uniq_tmpl_ids = []
    for x in tmpl_ids:
        if x not in uniq_tmpl_ids: uniq_tmpl_ids.append(x)

    saved_tmpl_fn=[]
    for tmpl in uniq_tmpl_ids:
        pdbcode=tmpl[:4]
        chainid=tmpl[4]
        dl_dir = pdbcode[1:3]
        #if (not os.path.exists(dl_dir)): os.system("mkdir %s"%dl_dir)
        dl_pdb = download_pdb(pdbcode, dl_dir)
        pdblines=open(dl_pdb).readlines()
        pdblines=keep_chain(pdblines,chainid)
        saved_pdb=save_pdb(pdblines, "%s.pdb"%tmpl)
        assert (os.path.exists(saved_pdb)), "File %s doesn't exist"%saved_pdb
        saved_tmpl_fn.append(saved_pdb)
    return saved_tmpl_fn, tmpl_fullnames

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Set up and run RosettaCM')
    parser.add_argument('--fasta', action="store", help='input fasta file',required=True)
    parser.add_argument('--alignment', action="store", default="", help='input alignment file')
    parser.add_argument('--alignment_format', action="store", default="grishin", help='input alignment file format (grishin, modeller, hhsearch, clustalw, or fasta)')
    parser.add_argument('--templates', nargs='*', default=[], help="input target templates")
    parser.add_argument('--rosetta_bin', action="store", default="~/Rosetta/main/source/bin", help="rosetta path")
    parser.add_argument('--build', action="store", default="default", help="build name in rosetta executable")
    parser.add_argument('--platform', action="store", default="linux", help="os platform name in rosetta executable")
    parser.add_argument('--compiler', action="store", default="gcc", help="compiler name in rosetta executable")
    parser.add_argument('--compiling_mode', action="store", default="release", help="compiling mode name in rosetta executable")
    parser.add_argument('--setup_script', action="store", default="", help="setup.pl path")
    parser.add_argument('--run', action="store_true", default=False, help="run the modeling process")
    parser.add_argument('-j', action="store",type=int, default=1, help="number of processors (requires GNU parallel in path!)")
    parser.add_argument('--keep_files', action="store_true", default=False, help="keep intermediate files")
    parser.add_argument('--run_dir', action="store", default="rosetta_cm", help="running dir")
    parser.add_argument('--use_dna', action="store_true", default=False, help="Add DNA from template to target")
    parser.add_argument('--verbose', action="store_true", default=False, help="print extra info")
    args = parser.parse_args()
    rosetta_bin = os.path.expanduser(args.rosetta_bin)
    rosetta_db = os.path.abspath("%s/../../database"%rosetta_bin)

	# fpd use parallel
    #if args.j > 1:
        #if args.build != 'mpi':
            #print "Add \"--build mpi\" for multiple processors"
            #sys.exit(-1)

    rosetta_exe_ext = "%s.%s%s%s"%(args.build, args.platform, args.compiler, args.compiling_mode)

    run_dir = os.path.abspath(args.run_dir)

    if args.verbose: print "Using rosetta from",rosetta_bin
    assert(os.path.exists(rosetta_bin)), "%s doesn't exist"%rosetta_bin
    assert(os.path.exists(rosetta_db)), "%s doesn't exist"%rosetta_db

    fasta_fn = os.path.abspath(args.fasta)
    if args.alignment != "":
        alignment_fn = os.path.abspath(args.alignment)
        assert(os.path.exists(alignment_fn)), "%s doesn't exist"%alignment_fn

    if len(args.templates) != 0:
        input_templates = [os.path.abspath(x) for x in args.templates]
        for x in input_templates:
            assert(os.path.exists(x)), "%s doesn't exist"%x


    curr_dir = os.getcwd()
    if (not args.keep_files):
        tempdir = tempfile.mkdtemp()
        if args.verbose: print "Switching to directory: %s"%tempdir
        os.chdir(tempdir)

    # get alignment
    converted_aln = os.path.abspath("converted_alignment.aln")
    if args.alignment != "":
        if args.alignment_format not in ["grishin", "modeller", "hhsearch", "clustalw", "fasta"]:
            print "Do not understand input alignment format: %s"%args.alignment_format
            sys.exit(-1)

        if args.alignment_format == "grishin":
            os.system("cp %s %s"%(alignment_fn, converted_aln))
            grishin_lines = open(converted_aln).readlines()
            alignment = Alignment()
            alignment.read_grishin(grishin_lines)
        elif args.alignment_format in ["modeller", "hhsearch", "clustalw", "fasta"]:
            alignment = Alignment()
            alignment.convert(alignment_fn, args.alignment_format, converted_aln, "grishin")

    if (not os.path.exists(run_dir)): os.makedirs(run_dir)

    assert (os.path.exists(converted_aln)), "File %s doesn't exist"%converted_aln
    os.system( "cp " + converted_aln + " " + run_dir )
    os.chdir(run_dir)

    if args.alignment_format == "clustalw":
        aln_file = open(converted_aln,'r')
        # Check that target sequence is the first one in the aln file
        for line in aln_file:
            buff = line.strip().split()
            if buff[0] == "##":
                print buff
                assert (buff[1] == fasta_fn.split('/')[-1].split('.')[0]), "\
		The first sequence ID in a clustalw alignment file must match the fasta file name. \n \
					   -> You may have to rename or swap the order of sequences in your alignment file."


    # special case, using a different script to set up RosettaCM
    if len(args.setup_script.strip()) > 0:
        # get template files
        templates, template_fullnames = download_templates(converted_aln)

        if os.path.exists(args.setup_script):
            command = "%s %s %s hybrid --autofrags"%(args.setup_script, args.fasta, converted_aln)
            os.system(command)
            # run RosettaCM
            exe_name="%s/rosetta_scripts.%s"%(args.rosetta_bin, rosetta_exe_ext)
            flags=glob.glob("hybrid/flags[0-9]_*")
            for flag in flags:
                command="%s @flags_common @%s -database %s -nstruct 10"%(exe_name, os.path.basename(flag), rosetta_db)
                if not args.run:
                    print "Run command in the \"hybrid\" directory: %s"%command
                    sys.exit(-1)

    # set up RosettaCM
    if (not os.path.exists("hybrid/hybridize0_C1.xml")):
        # get template files
        if len(args.templates) == 0:
            templates, template_fullnames = download_templates(converted_aln)
            template_filenames = " %s"*len(templates)%(tuple(templates))
        else:
            template_filenames = " %s"*len(input_templates)%(tuple(input_templates))
            template_fullnames = ["%s/%s.pdb"%(run_dir, os.path.basename(x)) for x in input_templates]

        thread_fullnames = ["%s/%s.pdb"%(run_dir, x.target_tag) for x in alignment.alignments ];

        # get partial thread
        exe_name="%s/partial_thread.%s.%s%s%s"%(args.rosetta_bin, args.build, args.platform, args.compiler, args.compiling_mode)
        command = "%s -database %s -mute all -in:file:fasta %s -in:file:alignment %s -in:file:template_pdb %s -ignore_unrecognized_res"%(exe_name, rosetta_db, fasta_fn, converted_aln, template_filenames)
        if args.verbose: print "Running %s"%command
        loglines = os.popen( command ).readlines()
        for name in thread_fullnames: assert (os.path.exists(name)), "File %s doesn't exist"%(name)

        flag_fn = "flags"
        xml_fn = "rosetta_cm.xml"
        silent_out = "rosetta_cm.silent"
        write_flags(flag_fn, fasta_fn, xml_fn, silent_out)
        write_xml(xml_fn, thread_fullnames)

        # run hybridization
        exe_name="%s/rosetta_scripts.%s"%(os.path.expanduser(args.rosetta_bin), rosetta_exe_ext)
        if args.j == 1:
            command="%s @%s -database %s -nstruct 20"%(exe_name, flag_fn, rosetta_db)
        else:
            command="parallel -j %d %s @%s -database %s -nstruct %d -out::suffix _{} ::: {1..%d}"%(args.j, exe_name, flag_fn, rosetta_db, math.ceil(20/args.j), args.j)

        if not args.run:
            print "Run command: %s"%command
            sys.exit(-1)
        else:
            assert (os.path.exists(exe_name)), "Executable %s doesn't exist"%exe_name

            if args.verbose: print "Running %s"%command
            loglines = os.popen( command ).readlines()

            # selection
            #score_lines= os.popen('grep SCORE: %s'%silent_fn)
            #score_columns =[]
            #for line in score_lines:
            #    if line.find(" score ") != -1:
            #        score_columns=line.strip().split()
            #        break

    os.chdir(curr_dir)
    if not args.keep_files:
        if args.verbose: print "Removing %s"%tempdir
        shutil.rmtree(tempdir)

