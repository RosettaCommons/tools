#!/usr/bin/ruby
require "optparse"
require "set"
require "pathname"
def get_fasta(fn)
  fasta_fl = File.open(fn)
  line = fasta_fl.gets()
  fasta = ""
  while(line = fasta_fl.gets())
    fasta += line.strip
  end
  return(fasta)
end

def generate_profile(in_aln_fn,fasta_fn,out_profile_fn)
  in_aln = File.open(in_aln_fn,"r")
  out_profile = File.open(out_profile_fn,"w")
  aa_order="ACDEFGHIKLMNPQRSTVWY"
  native = get_fasta(fasta_fn)
  cts = Array.new
  for ii in 0..native.size-1 do
    cts.push([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
  end
  #collect all counts
  while(line=in_aln.gets)
    seq=in_aln.gets.strip
    for ii in 0..seq.size-1 do
      aa = seq[ii]
      pos = aa_order.index(aa)
      if(pos != nil)
        cts[ii][pos]+=1
      end
    end
  end
  #print out counts
  out_profile << ("aa\t")
  for jj in 0..aa_order.size
    out_profile << ("#{aa_order[jj]}    \t")
  end
  out_profile << ("\n")
  for ii in 0..native.size-1 do
    out_profile << ("#{native[ii]}\t")
    total_cts = 0
    for jj in 0..aa_order.size-1
      total_cts+=cts[ii][jj]
  end
    tmp = 0
    for jj in 0..aa_order.size-1
      score = -Math.log((cts[ii][jj]+1).to_f/(total_cts+20).to_f)
      out_profile << sprintf("%.2f\t",score)
    end
   out_profile << ("\n")
  end
end

def get_fasta(fn)
  fasta_fl = File.open(fn)
  line = fasta_fl.gets()
  fasta = ""
  while(line = fasta_fl.gets())
    fasta += line.strip
  end
  return(fasta)
end

def generate_cstfile(fasta,cstFileOut_fn)
  cstFileOut = File.open(cstFileOut_fn,"w")
  fasta = get_fasta(fasta)
  for ii in 0..(fasta.size-1)
    cstFileOut << "SequenceProfile #{ii+1} profile\n"
  end 
end

def generate_blueprint(pdb_fn,fasta_fn,dssp_script)
    dssp_fn = "start.ss_dssp"
    system("#{dssp_script} #{pdb_fn} > #{dssp_fn}")
    fasta = get_fasta(fasta_fn)
    dssp = get_fasta(dssp_fn)
    blueprint_fn = "start.blueprint"
    blueprint = File.open(blueprint_fn,"w")
    for ii in 0..fasta.size-1 
        blueprint<< "#{ii+1} #{fasta[ii].chr} #{dssp[ii].chr} .\n"
    end
end

options = {}
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#This line should automatically detect the location of the script
structProf_path = File.expand_path(File.dirname(__FILE__))
pdb2vall_path = Pathname(structProf_path).parent
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

optparse = OptionParser.new do |opts|
    options[:pdb] = ""
    opts.on('-p','--pdb FILE',String,"input pdb:REQ") do |f|
        options[:pdb] = f
    end
    opts.on( '-h','--help', 'Display this screen') do
        puts opts
        exit
    end
end
optparse.parse!
pdb_path = options[:pdb]
pdb = File.basename(pdb_path,".pdb")
Dir.mkdir(pdb)
system("cp #{pdb_path} #{pdb}/start.pdb")
Dir.chdir(pdb)
system("#{pdb2vall_path}/pdb_scripts/pdb2fasta.py start.pdb > start.fasta")
system("#{pdb2vall_path}/structure_profile_scripts/DEPTH-CLONE-2.8.7/DEPTH -thread 20 -i start.pdb -o start.depth")
system("#{pdb2vall_path}/structure_profile_scripts/make_sequence_fragments.pl -n_frags 100 -n_candidates 1000 -frag_sizes 9 -depth start.depth-residue.depth -native start.pdb start.fasta")
system("#{pdb2vall_path}/structure_profile_scripts/make_alignment_from_fragfile_rmsCutoff.pl start.fasta")
generate_profile("start.100.9mers.ali.fasta","start.fasta","profile")
generate_cstfile("start.fasta","MSAcst")
dssp_script = "#{pdb2vall_path}/structure_profile_scripts/dssp2threestateSS.pl"
generate_blueprint("start.pdb","start.fasta",dssp_script)
Dir.chdir("..")
