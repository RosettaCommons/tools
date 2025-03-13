#!/usr/bin/env python3.9
# ./easy_fragments.py --secfile 2jsv.ss 2jsv.ss --seqfile 2jsv.fasta --vall vall --bin ~/src.git/Rosetta/main/source/bin/
#
# ./easy_fragments.py  --sequence MQYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE --secondary CEEEEEECCCCCCEEEEEECCHHHHHHHHHHHHHHHHCCCCEEEEECCCCEEEEEC --vall ~/src.git/Rosetta/tools/fragment_tools/vall.jul19.2011.gz --bin ~/src.git/Rosetta/main/source/bin/

import argparse
import os
import subprocess

def run_command_with_flagfile(command: str, flags: dict):
    with open("flagfile", "w") as f:
        for flag, value in flags.items():
            if isinstance(value, list):
                for v in value:
                    print(flag+"\t\t"+v, file=f)
            else:
                print(flag+"\t\t"+value, file=f)
    
    try:
        result = subprocess.run([command, "@flagfile"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        print("Output:", result.stdout)
    except subprocess.CalledProcessError as e:
        print("Error:", e.stderr)


def read_fasta(filename):
    """Read a FASTA file and return the sequence as a string."""
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    sequence = "".join(line.strip() for line in lines if not line.startswith('>'))
    return sequence


def convert_to_ss2(sequence, structure):
    """Convert amino acid sequence and secondary structure to ss2 format."""
    header = "# PSIPRED VFORMAT (PSIPRED V2.6 by David Jones)\n"
    output_lines = [header]
    
    for i, (aa, ss) in enumerate(zip(sequence, structure), start=1):
        probabilities = [0.0, 0.0, 0.0]
        if ss == 'C':
            probabilities[0] = 1.0
        elif ss == 'H':
            probabilities[1] = 1.0
        elif ss == 'E':
            probabilities[2] = 1.0
        
        output_lines.append(f"{i:4d} {aa} {ss}   {probabilities[0]:.3f}  {probabilities[1]:.3f}  {probabilities[2]:.3f}")
    
    return "\n".join(output_lines)


def write_quota_config(ss2_files, pdb_file=None):
    n_pools = len(ss2_files)
    with open("quota.def", "w") as f:
        fraction = 1.0 / n_pools
        print("#pool_id	pool_name	fraction",file=f)
        for i in range(n_pools):
            print(f"{i}\tpool{i}\t{fraction}",file=f)

    with open("quota-protocol.wghts", "w") as f:
        print("# score name          priority	wght   min_allowed  extras",file=f)
        for i in range(len(ss2_files)):
            print(f"SecondarySimilarity	350	100.0	-	pool{i}",file=f)
            print(f"RamaScore		150	100.0	-	pool{i}",file=f)
        if pdb_file:
            print(f"FragmentCrmsd		50	0.0	{pdb_file}",file=f)



def run_fragment_picker(fasta_file, ss2_files, vall_pathfile, rosetta_bin_path, **kwargs):

    extension = kwargs.get("extension", "linuxgccrelease").strip('.')
    pdb_file = kwargs.get("pdb_file", None)
    lengths = kwargs.get("lengths", [3,9])
    command = os.path.join(rosetta_bin_path, "fragment_picker."+extension)
    ss2_inputs = " ".join([f"{fname} pool{i}" for (i, fname) in enumerate(ss2_files)])
    flags = {"-frags::n_candidates": "1000",
             "-frags::n_frags": "200",
             "-frags::scoring::config": "quota-protocol.wghts",
             "-frags::ss_pred": ss2_inputs,
             "-in::file::vall": vall_pathfile,
             "-frags::frag_sizes": " ".join([str(l) for l in lengths]),
             "-in:file:fasta": fasta_file,
	     "-frags::picking::quota_config_file": "quota.def",
	     "-frags::scoring::config": "quota-protocol.wghts",
	     "-out::level": "1000",
	     "-out::mute": ["protocols.frag_picker.VallProvider",
	    	    "protocols.frag_picker.quota.QuotaCollector",
	    	    "protocols.frag_picker.quota.SecondaryStructurePool",
	    	    "protocols.frag_picker.scores.SecondarySimilarity",
                    "protocols.frag_picker.quota.QuotaSelector",
	    	    "core",
	    	    "basic",
                ]
            }
    if pdb_file:
        flags["-in::file::s"] = pdb_file
        flags["-frags::describe_fragments"] = "frags.fsc"
    run_command_with_flagfile(command, flags)


def main():
    parser = argparse.ArgumentParser(description="Convert amino acid sequence and secondary structure to ss2 format.")
    parser.add_argument("--sequence", type=str, help="Amino acid sequence as a string.")
    parser.add_argument("--secondary", type=str, help="Secondary structure sequence as a string.")
    parser.add_argument("--seqfile", type=str, help="FASTA file containing the amino acid sequence.")
    parser.add_argument("--secfile", type=str, nargs='+', help="FASTA file containing the secondary structure (HEC sequence).")
    parser.add_argument("--native", type=str, help="If the native protein structure is provided, a score table will be printed, including crmsd of fragments")
    parser.add_argument("--vall", type=str, help="Database file with protein structures to extract fragments from.")
    parser.add_argument("--bin", type=str, help="Path to Rosetta binaries")
    parser.add_argument("-n", type=int, nargs='+', default=[3,9], help="Size(s) of fragments to pick")
    args = parser.parse_args()
    
    fasta_file = None
    if args.seqfile:
        sequence = read_fasta(args.seqfile)
        fasta_file = args.seqfile
    elif args.sequence:
        sequence = args.sequence
        fasta_file = "input.fasta"
        with open(fasta_file,"w") as f:
            print(f">input\n{sequence}\n", file=f)
    else:
        print("Error: No sequence provided.")
        return
    
    secondary = []
    if args.secfile:
        for ss2 in args.secfile:
            secondary.append(read_fasta(ss2))
    elif args.secondary:
        secondary.append(args.secondary)
    else:
        print("Error: No secondary structure provided.")
        return
    if len(secondary) == 1:		# duplicate secondary if just one was given: the quota protocol need at least 2 pools
        secondary.append(secondary[0])

    if len(sequence) != len(secondary[0]):
        print("Error: Sequence and secondary must be of the same length.")
        return
    
    if args.vall:
        vall = args.vall
    else:
        print("Error: No vall database provided.")
        return

    if args.bin:
        bin_path = args.bin
    else:
        print("Error: No binary path provided.")
        return

    ss2_fnames = []
    for (i,si) in enumerate(secondary):
        ss2_output = convert_to_ss2(sequence, si)
        fname = f"secondary{i}.ss2"
        ss2_fnames.append(fname)
        with open(fname,"w") as f:
            print(ss2_output, file=f)

    pdb_file = None if not args.native else args.native

    write_quota_config(ss2_fnames, pdb_file)
    extra_args = {"pdb_file": pdb_file, "lengths": args.n}
    run_fragment_picker(fasta_file, ss2_fnames, vall, bin_path, **extra_args)


if __name__ == "__main__":
    main()
