#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   antibody.py
## @brief  Pre-processing script for antibody protocol
## @author Sergey Lyskov
## @author Modified by Daisuke Kuroda

import os, sys, re, json, commands, shutil

from optparse import OptionParser, IndentedHelpFormatter

_script_path_ = os.path.dirname( os.path.realpath(__file__) )

_framework_names_ = ['FRL', 'FRH', 'light', 'heavy', 'L1', 'L2', 'L3', 'H1', 'H2', 'H3', 'light_heavy']


'''
_alignment_legend_to_pretty_legend = {
    'subject-id': 'Subject id',
    'resolution':'Resolution',
    '%-identity': '% identity',
    'alignment-length': 'alignment length',
    'mismatches': 'mismatches',
    'gap-opens': 'gap-opens',
    'q.start': 'q.start',
    'q.end': 'q.end',
    's.start': 's.start',
    's.end': 's.end',
    'evalue': 'evalue',
    'bit-score': 'bit-score'
}'''


def main(args):
    ''' Script for preparing detecting antibodys and preparing info for Rosetta protocol.
    '''
    parser = OptionParser(usage="usage: %prog [OPTIONS] [TESTS]")
    parser.set_description(main.__doc__)

    parser.add_option('-L','--light-chain',
      action="store",
      help="Specify the light chain.",
    )

    parser.add_option('-H','--heavy-chain',
      action="store",
      help="Specify the heavy chain.",
    )

    parser.add_option('--prefix',
      action="store", default='output/',
      help="Prefix for output files. Should be dir name. Default is ./ string.",
    )

    parser.add_option('--blast',
      action="store", default='blastp',
      help="Specify path+name for 'blastall' executable. Default is blastp [blast+].",
    )

    parser.add_option('--profit',
      action="store", default='profit',
      help="Specify path+name for 'ProFIt' executable. Default is profit.",
    )

    parser.add_option('--blast-database',
      action="store", default=None,
      help="Specify path of blast database dir.",
    )

    parser.add_option('--antibody-database',
      action="store", default=None,
      help="Specify path of antibody database dir.",
    )

    parser.add_option('--rosetta-database',
      action="store", default=None,
      help="Specify path of rosetta database dir.",
    )

    parser.add_option('--homologue_exclusion',
      default=200, type="int",
      help="Specify the cut-off for homologue exclusion during template selections.",
    )

    parser.add_option('--homologue_exclusion_cdr',
      default=200, type="int",
      help="Specify the cut-off for homologue exclusion during ***CDR*** template selections.",
    )

    parser.add_option('--homologue_exclusion_fr',
      default=200, type="int",
      help="Specify the cut-off for homologue exclusion during ***FR*** template selections.",
    )

    parser.add_option('--rosetta-bin',
      action="store", default=None,
      help="Specify path to 'rosetta/source/bin' dir where antibody_graft', idealize and relax executable expected to be found. Default is '$ROSETTA/main/source/bin, then <script location>/bin' (plasce symlink there) and if not found corresponding steps will be skipped.",
    )

    parser.add_option('--rosetta-platform',
      action="store", default=None,
      help="Specify full extra+compier+build type for rosetta biniaries found in --rosetta-bin. For example use static.linuxgccrelease for static build on Linux. Default is dynamic release build of current OS",
    )

    parser.add_option("--idealize",
      default=0, type="int",
      help="Specify if idealize protocol should be running on final model [0/1]. (default: 0, which mean do not run idealize protocol)",
    )

    parser.add_option("--relax",
      default=1, type="int",
      help="Specify if relax protocol should be running on final model [0/1]. (default: 1, which mean run relax protocol)",
    )

    parser.add_option("--quick","-q",
      default=0, type="int",
      help="Specify fast run (structure will have clashes).  Prevents stem optimization and turns off relax, idealize.",
    )

    #parser.add_option('--rosetta-options',
    #  action="store", default='',
    #  help="Specify extra options for antibody_graft run.",
    #)

    for name in _framework_names_:
        parser.add_option('--' + name,
          action="store", default='',
          help="Specify path or PDB code for %s template. If specified this will overwrite blast selection." % name,
        )

    parser.add_option('--self-test',
      action="store_true",
      help="Perform self test by using data in test/ dir and exit.",
    )

    parser.add_option('--self-test-dir',
      action="store", default='self-test/',
      help="Specify path for self test dir [default:self-test/].",
    )

    parser.add_option('-v', "--verbose",
      action="store_true", default=False,
      help="Generate verbose output.",
    )

    # Filter list of 'filter-function:default state' pairs here, and  extend it to add more filters
    global Filters;  Filters = { filter_by_sequence_length:True, filter_by_alignment_length:False, filter_by_template_resolution:True,
                                 filter_by_outlier:True, filter_by_template_bfactor:True, filter_by_sequence_homologue:True }

    for f in Filters: parser.add_option('--' + f.func_name.replace('_', '-'), type="int", default=int(Filters[f]),
                                        help="Boolen option [0/1] that control filetering results with %s function." % f.func_name)

    (options, args) = parser.parse_args(args=args[1:])
    global Options;  Options = options

    global sid_cutoff
    global sid_cutoff_cdr
    global sid_cutoff_fr

    sid_cutoff     = Options.homologue_exclusion
    sid_cutoff_cdr = Options.homologue_exclusion_cdr
    sid_cutoff_fr  = Options.homologue_exclusion_fr

    global frlh_info
    global frl_info
    global frh_info

    frlh_info, legend = {}, ''
    for l in file( _script_path_ + '/info/frlh_info' ):
        if l.startswith('# '): legend = l[2:].split()
        elif len(l)>8: frlh_info[l.split()[0]] =  dict( zip(legend, l.split() ) )

    frl_info, legend = {}, ''
    for l in file( _script_path_ + '/info/frl_info' ):
        if l.startswith('# '): legend = l[2:].split()
        elif len(l)>8: frl_info[l.split()[0]] =  dict( zip(legend, l.split() ) )

    frh_info, legend = {}, ''
    for l in file( _script_path_ + '/info/frh_info' ):
        if l.startswith('# '): legend = l[2:].split()
        elif len(l)>8: frh_info[l.split()[0]] =  dict( zip(legend, l.split() ) )

    #os.getcwd() #os.chdir()
    script_dir = os.path.dirname(__file__)

    if Options.prefix and Options.prefix[-1] != '/': Options.prefix += '/'

    if not Options.blast_database:    Options.blast_database    = script_dir + '/blast_database'
    if not Options.antibody_database: Options.antibody_database = script_dir + '/antibody_database'
    if not Options.rosetta_bin:
        if 'ROSETTA' in os.environ:
            Options.rosetta_bin = os.path.abspath(os.environ['ROSETTA']) + '/main/source/bin'
        else: Options.rosetta_bin = script_dir + '/bin'
    if not Options.rosetta_database:
        if os.path.isdir(script_dir + '/rosetta_database'):
            Options.rosetta_database = os.path.abspath( script_dir + '/rosetta_database' )
        elif 'ROSETTA3_DB' in os.environ:
            Options.rosetta_database = os.path.abspath(os.environ['ROSETTA3_DB'])
        elif os.path.isdir(os.environ['HOME'] + '/rosetta_database'):
            Options.rosetta_database = os.path.abspath(os.environ['HOME'] + '/rosetta_database')

    Options.blast_database    = os.path.abspath( Options.blast_database )
    Options.antibody_database = os.path.abspath( Options.antibody_database )
    Options.rosetta_bin       = os.path.abspath( Options.rosetta_bin )

    if not Options.rosetta_platform:
        if sys.platform.startswith("linux"): Options.rosetta_platform = 'linuxgccrelease'
        elif sys.platform == "darwin" : Options.rosetta_platform = 'macosgccrelease'
        else: Options.rosetta_platform = '_unknown_'

    if Options.self_test: self_test();  return

    if not(options.light_chain and options.heavy_chain):
        print 'Script for preparing detecting antibodys and preparing info for Rosetta protocol.'
        print 'At miminum you need to specify options --light-chain and --heavy-chain.'
        print 'For full list of options run "antibody.py --help"\nERROR: No input chains was specifiede... exiting...'
        sys.exit(1)

    light_chain = read_fasta_file(options.light_chain)
    heavy_chain = read_fasta_file(options.heavy_chain)

    print 'Rosetta Antibody script [Python, version 2.0]. Starting...'

    prefix_details = Options.prefix + 'details/'
    if not os.path.isdir(Options.prefix): print 'Could not find %s... creating it...' % Options.prefix;  os.makedirs(Options.prefix)
    if not os.path.isdir(prefix_details): print 'Could not find %s... creating it...' % prefix_details;  os.makedirs(prefix_details)

    print 'Prefix:', Options.prefix
    print 'Blast database:', Options.blast_database
    print 'Antibody database:', Options.antibody_database
    print 'rosetta_bin:', Options.rosetta_bin, ' [platform: %s]' % Options.rosetta_platform
    print 'rosetta_database:', Options.rosetta_database


    antibody_database_files = os.listdir(Options.antibody_database)
    bz2ipped_files  = [f for f in antibody_database_files if f.endswith('.bz2') and f[:-4] not in antibody_database_files]
    if bz2ipped_files:
        print 'Unpacking rosetta_database files (this need to be done only once)...'
        commandline = 'cd %s && bunzip2 -k %s' % (Options.antibody_database, ' '.join(bz2ipped_files))
        res, output = commands.getstatusoutput(commandline)
        if Options.verbose and not res: print commandline+'\n', output
        if res: print commandline+'\n', 'ERROR: Unpacking antibody database files failed with code %s and message: %s' % (res, output);  sys.exit(1)


    if   sid_cutoff > 100 and sid_cutoff_cdr > 100 and sid_cutoff_fr > 100:
        print 'Using full antibody database (no homolog exclusion)'
    elif sid_cutoff_cdr <= 100 or sid_cutoff_fr <= 100:
        if sid_cutoff_cdr <= 100:
            print '\n!!! Homologues will be excluded with %s SID cut-off during ***CDR*** template selections !!!' % sid_cutoff_cdr
        else:
            print '\n!!! Homologues will not be excluded during ***CDR*** template selection (default)    !!!'

        if sid_cutoff_fr <= 100:
            print '\n!!! Homologues will be excluded with %s SID cut-off during ***FR*** template selections !!!' % sid_cutoff_fr
        else:
            print '\n!!! Homologues will not be excluded during ***FR*** template selection (default)    !!!'
    elif sid_cutoff <= 100:
        print     '\n!!! Homologues will be excluded with %s SID cut-off during template selections !!!' % sid_cutoff
        sid_cutoff_cdr = sid_cutoff
        sid_cutoff_fr = sid_cutoff

    if Options.quick:
        Options.relax = 0
        Options.idealize = 0

    print 'Idealize:', bool(Options.idealize)
    print 'Relax:', bool(Options.relax)

    for name in _framework_names_:
        if getattr(Options, name): print 'Custom %s template:' % name, getattr(Options, name)
    print
    print "Light chain: %s" % light_chain
    print "Heavy chain: %s" % heavy_chain

    CDRs = IdentifyCDRs(light_chain, heavy_chain)
    CDRs.update( Extract_FR_CDR_Sequences(**CDRs) )
    CDRs['light_heavy'] = CDRs['light'] + CDRs['heavy']
    write_results(CDRs, prefix_details)

    if Options.verbose: print 'CDR:', json.dumps(CDRs, sort_keys=True, indent=2)
    else:
        c = dict(CDRs);  c.pop('numbering_L');  c.pop('numbering_H')
        #print 'CDR:', json.dumps(c, sort_keys=True, indent=2)

    alignment, legend = run_blast(CDRs, prefix=prefix_details, blast=Options.blast, blast_database=Options.blast_database, verbose=Options.verbose)

    create_virtual_template_pdbs(prefix=prefix_details)
    thread_template_pdbs(CDRs, prefix=prefix_details)
    superimpose_templates(CDRs, prefix=prefix_details)

    if Options.rosetta_database:
        run_rosetta(CDRs, prefix=Options.prefix, rosetta_bin=Options.rosetta_bin, rosetta_platform=Options.rosetta_platform, rosetta_database=Options.rosetta_database)
    else:
        print 'Rosetta database was not found... skipping rosetta run...'

    results = { 'cdr':{}, 'numbering':{}, 'alignment':alignment, 'alignment.legend':legend }
    for k in [i for i in CDRs if not i.startswith('numbering_')]: results['cdr'][k] = CDRs[k]
    for n in [i for i in CDRs if i.startswith('numbering_')]: results['numbering'][n] = CDRs[n]

    with file(Options.prefix+'results.json', 'w') as f: json.dump(results, f, sort_keys=True, indent=2)
    print 'Done!'


def read_fasta_file(file_name):
    return ''.join( [l.rstrip() for l in file(file_name) if not l.startswith('>') ] ) . replace(' ', '')


def write_fasta_file(file_name, data, prefix):
    with file(prefix+file_name+'.fasta', 'w') as f: f.write('> %s\n' % file_name);  f.write(data); f.write('\n')


def write_results(CDRs, prefix):
    #with file(prefix+'cdr.json', 'w') as f: json.dump(CDRs, f, sort_keys=True, indent=2)
    for k in [i for i in CDRs if not i.startswith('numbering_')]: write_fasta_file(k, CDRs[k], prefix)
    for n in [i for i in CDRs if i.startswith('numbering_')]:
        with file(prefix + n +'.txt', 'w') as f:
            f.write('\n'.join( [ '%s %s' % (CDRs[n][k], k) for k in sorted(CDRs[n].keys(), key=lambda x: (int_(x), x) ) ]) + '\n')


def safelen(seq):
    return 0 if not seq else len(seq)

def IdentifyCDRs(light_chain, heavy_chain):
    ''' Identift CDR region and return them as dict with keys: 'FR_H1', 'FR_H2', 'FR_H3', 'FR_H4', 'FR_L1', 'FR_L2', 'FR_L3', 'FR_L4', 'H1', 'H2', 'H3', 'L1', 'L2', 'L3'
    '''
    light_first, light_second = (light_chain[:60], light_chain[50:50+60]) if len(light_chain) > 120 else (light_chain[:60], light_chain[50:])
    heavy_first, heavy_second = (heavy_chain[:60], heavy_chain[50:50+90]) if len(heavy_chain) > 140 else (heavy_chain[:60], heavy_chain[50:])

    # L1
    res = re.search( r'C[A-Z]{1,17}(WYL|WLQ|WFQ|WYQ|WYH|WVQ|WVR|WWQ|WVK|WYR|WLL|WFL|WVF|WIQ|WYR|WNQ|WHL|WHQ|WYM|WYY)', light_first)
    L1 = res.group()[1:-3] if res else False
    print "L1 detected: ", L1, " (",safelen(L1),"residues )"

    # L3
    res = re.search( r'C[A-Z]{1,15}(F|V|S)G[A-Z](G|Y)', light_second)
    L3 = res.group()[1:-4] if res else False
    print "L3 detected: ", L3, " (",safelen(L3),"residues )"

    if L1 and L3:
        L1_start = light_chain.index(L1)
        L1_end = L1_start + len(L1) - 1

        L2_start = L1_end + 16
        L2_end = L2_start + 7 - 1

        L3_start = light_chain.index(L3)
        L3_end = L3_start + len(L3) - 1

        L2 = light_chain[L2_start:L2_start+7]  # L2 is identified here. Current implementation can deal with only 7-resiue L2
        print "L2 detected: ", L2, " (",safelen(L2),"residues )"

        FR_L1 = light_chain[:L1_start]
        FR_L2 = light_chain[L1_end + 1 : L1_end + 1+ 15]
        FR_L3 = light_chain[L2_end+1 : L2_end+1 +L3_start - L2_end - 1 ]
        FR_L4 = light_chain[L3_end + 1 : L3_end + 1 + 12]

        print "FR_L1: ", FR_L1
        print "FR_L2: ", FR_L2
        print "FR_L3: ", FR_L3
        print "FR_L4: ", FR_L4
        print "L segments: ",FR_L1,L1,FR_L2,L2,FR_L3,L3,FR_L4

        # Light chain sub-type classification. This is useful in the future. But currently this is not used.
        # ... skipped, see Google doc for details
        # FR classification by AHo. This might be useful in the future. But currently this is not used.
        # ... skipped, see Google doc for details

    # H1
    res = re.search( r'C[A-Z]{1,16}(W)(I|V|F|Y|A|M|L|N|G)(R|K|Q|V|N|C|G)(Q|K|H|E|L|R)', heavy_first) # jeff's mod for ATHM set
    #res = re.search( r'C[A-Z]{1,16}(W)(I|V|F|Y|A|M|L|N|G)(R|K|Q|V|N|C)(Q|K|H|E|L|R)', heavy_first)
    H1 = res.group()[4:-4] if res else False
    print "H1 detected: ", H1, " (",safelen(H1),"residues )"

    # H3
    res = re.search( r'C[A-Z]{1,33}(W)(G|A|C)[A-Z](Q|S|G|R)', heavy_second)
    H3 = res.group()[3:-4] if res else False  #H3_and_stem = res.group()[0:-4] if res else False
    print "H3 detected: ", H3, " (",safelen(H3),"residues )"

    if H1 and H3:
        H1_start = heavy_chain.index(H1)
        H1_end = H1_start + len(H1) - 1

        H3_start = heavy_chain.index(H3)
        H3_end = H3_start + len(H3) - 1

        H2_start = H1_end + 15
        H2_end = H3_start - 33

        H2 = heavy_chain[H2_start:H2_start + H2_end-H2_start+1]
        print "H2 detected: ", H2, " (",len(H2),"residues )"

        FR_H1 = heavy_chain[:H1_start]
        FR_H2 = heavy_chain[H1_end + 1: H1_end + 1 + H2_start - H1_end - 1]
        FR_H3 = heavy_chain[H2_end + 1: H2_end + 1 + H3_start - H2_end - 1]
        FR_H4 = heavy_chain[H3_end + 1: H3_end + 1 + 12]

        print "FR_H1: ", FR_H1
        print "FR_H2: ", FR_H2
        print "FR_H3: ", FR_H3
        print "FR_H4: ", FR_H4
        print "H segments: ",FR_H1,H1,FR_H2,H2,FR_H3,H3,FR_H4

    if not (L1 and L3 and H1 and H3):
        if not L1: print 'ERROR: CDR L1 cannot be recognized !!!  L1 pattern: C[A-Z]{1,17}(WYL|WLQ|WFQ|WYQ|WYH|WVQ|WVR|WWQ|WVK|WYR|WLL|WFL|WVF|WIQ|WYR|WNQ|WHL|WHQ|WYM|WYY)'
        if not L3: print 'ERROR: CDR L3 cannot be recognized !!!  L3 pattern: C[A-Z]{1,15}(F|V|S)G[A-Z](G|Y)'
        if not H1: print 'ERROR: CDR H1 cannot be recognized !!!  H1 pattern: C[A-Z]{1,16}(W)(I|V|F|Y|A|M|L|N|G)(R|K|Q|V|N|C)(Q|K|H|E|L|R)'
        if not H3: print 'ERROR: CDR H3 cannot be recognized !!!  H3 pattern: C[A-Z]{1,33}(W)(G|A|C)[A-Z](S|G|R)'
        sys.exit(1)

    res = dict(L1=L1, L2=L2, L3=L3, H1=H1, H2=H2, H3=H3,  FR_L1=FR_L1, FR_L2=FR_L2, FR_L3=FR_L3, FR_L4=FR_L4,  FR_H1=FR_H1, FR_H2=FR_H2, FR_H3=FR_H3, FR_H4=FR_H4)
    #if Options.verbose: print 'L1: %(L1)s\nL2: %(L2)s\nL3: %(L3)s\nH1: %(H1)s\nH2: %(H2)s\nH3: %(H3)s' % res

    return res


def int_(s): return int( re.sub('[A-Z]', '', s) )  #v = int( re.sub('[A-Z]', '', new_number_FR_L1) )  # $new_number_FR_L1[$i] =~ s/[A-Z]//     #new_number_FR_L1[i] = string.translate(new_number_FR_L1[i], None, string.ascii_letters)

def Extract_FR_CDR_Sequences(L1='', L2='', L3='', H1='', H2='', H3='', FR_L1='', FR_L2='', FR_L3='', FR_L4='', FR_H1='', FR_H2='', FR_H3='', FR_H4=''):
    # LIGHT CHAIN
    # FR_L1	How can we handle missing residue in C/N-terminals?
    if re.search( r'[A-Z][QE][A-Z]{9}[A-Z][A-Z]{4}[LVIMF][A-Z]C', FR_L1): # Change G to [A-Z] (3G04)
        if   len(FR_L1) == 19: new_number_FR_L1="5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
        elif len(FR_L1) == 20: new_number_FR_L1="4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
        elif len(FR_L1) == 21: new_number_FR_L1="3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
        elif len(FR_L1) == 22: new_number_FR_L1="2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
        elif len(FR_L1) == 23: new_number_FR_L1="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
        elif len(FR_L1) == 24:
            #new_number_FR_L1="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
            FR_L1 = FR_L1[1:]  # Remove 0th residue 10/24/2012
            new_number_FR_L1="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
        else: print "ERROR: FR_L1 matches [A-Z][QE][A-Z]{9}[A-Z][A-Z]{4}[LVIMF][A-Z]C but length",len(FR_L1),"is not between 19 and 24"

    elif re.search( r'[A-Z][QE][A-Z]{8}[A-Z][A-Z]{4}[LVIMF][A-Z]C', FR_L1):  # Change G to [A-Z] (3G04)
        if   len(FR_L1) == 19: new_number_FR_L1="4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
        elif len(FR_L1) == 20: new_number_FR_L1="3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
        elif len(FR_L1) == 21: new_number_FR_L1="2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
        elif len(FR_L1) == 22: new_number_FR_L1="1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
        elif len(FR_L1) == 23: new_number_FR_L1="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
        elif len(FR_L1) == 24:
            #new_number_FR_L1="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
            FR_L1 = FR_L1[1:]  # Remove 0th residue 10/24/2012
            new_number_FR_L1="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23"
        else: print "ERROR: FR_L1 matches [A-Z][QE][A-Z]{8}[A-Z][A-Z]{4}[LVIMF][A-Z]C but length",len(FR_L1),"is not between 19 and 24"
    else:
        print 'ERROR: Current code could not assign Chothia numbering of FR_L1 in the query sequence!!! Exiting...'

    # L1
    if   len(L1) ==  8: new_number_L1="24,25,26,27,28,29,30,34"
    elif len(L1) ==  9: new_number_L1="24,25,26,27,28,29,30,33,34"
    elif len(L1) == 10: new_number_L1="24,25,26,27,28,29,30,32,33,34"
    elif len(L1) == 11: new_number_L1="24,25,26,27,28,29,30,31,32,33,34"
    elif len(L1) == 12: new_number_L1="24,25,26,27,28,29,30,30A,31,32,33,34"
    elif len(L1) == 13: new_number_L1="24,25,26,27,28,29,30,30A,30B,31,32,33,34"
    elif len(L1) == 14: new_number_L1="24,25,26,27,28,29,30,30A,30B,30C,31,32,33,34"
    elif len(L1) == 15: new_number_L1="24,25,26,27,28,29,30,30A,30B,30C,30D,31,32,33,34"
    elif len(L1) == 16: new_number_L1="24,25,26,27,28,29,30,30A,30B,30C,30D,30E,31,32,33,34"
    elif len(L1) == 17: new_number_L1="24,25,26,27,28,29,30,30A,30B,30C,30D,30E,30F,31,32,33,34"
    else: print "ERROR: L1 length",len(L1),"is not between 8 and 17"

    # FR_L2
    if len(FR_L2) == 15: new_number_FR_L2="35,36,37,38,39,40,41,42,43,44,45,46,47,48,49"
    else: print "ERROR: FR_L2 length",len(FR_L2),"is not 15"

    # L2
    if   len(L2) ==  7: new_number_L2="50,51,52,53,54,55,56"
    elif len(L2) ==  8: new_number_L2="50,51,52,53,54,54A,55,56"
    elif len(L2) ==  9: new_number_L2="50,51,52,53,54,54A,54B,55,56"
    elif len(L2) == 10: new_number_L2="50,51,52,53,54,54A,54B,54C,55,56"
    elif len(L2) == 11: new_number_L2="50,51,52,53,54,54A,54B,54C,54D,55,56"
    else: print "ERROR: L2 length",len(L2),"is not between 7 and 11"

    # FR_L3
    if   len(FR_L3) == 32: new_number_FR_L3="57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88"
    elif len(FR_L3) == 33: new_number_FR_L3="57,58,59,60,61,62,63,64,65,66,66A,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88"
    elif len(FR_L3) == 34: new_number_FR_L3="57,58,59,60,61,62,63,64,65,66,66A,66B,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88"
    else: print "ERROR: FR_L3 length",len(FR_L3),"is not between 32 and 34"

    # L3
    if   len(L3) ==  5: new_number_L3="89,90,91,92,97"
    elif len(L3) ==  6: new_number_L3="89,90,91,92,93,97"
    elif len(L3) ==  7: new_number_L3="89,90,91,92,93,94,97"
    elif len(L3) ==  8: new_number_L3="89,90,91,92,93,94,95,97"
    elif len(L3) ==  9: new_number_L3="89,90,91,92,93,94,95,96,97"
    elif len(L3) == 10: new_number_L3="89,90,91,92,93,94,95,95A,96,97"
    elif len(L3) == 11: new_number_L3="89,90,91,92,93,94,95,95A,95B,96,97"
    elif len(L3) == 12: new_number_L3="89,90,91,92,93,94,95,95A,95B,95C,96,97"
    elif len(L3) == 13: new_number_L3="89,90,91,92,93,94,95,95A,95B,95C,95D,96,97"
    elif len(L3) == 14: new_number_L3="89,90,91,92,93,94,95,95A,95B,95C,95D,95E,96,97"
    elif len(L3) == 15: new_number_L3="89,90,91,92,93,94,95,95A,95B,95C,95D,95E,95F,96,97"
    else: print "ERROR: L3 length",len(L3),"is not between 5 and 15"

    # FR_L4
    new_number_FR_L4="98,99,100,101,102,103,104,105,106,107,108,109"

    # HEAVY CHAIN
    # FR_H1	How can we handle missing residue in C/N-terminals?
    if   len(FR_H1) == 16: new_number_FR_H1="10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
    elif len(FR_H1) == 17: new_number_FR_H1="9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
    elif len(FR_H1) == 18: new_number_FR_H1="8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
    elif len(FR_H1) == 19: new_number_FR_H1="7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
    elif len(FR_H1) == 20: new_number_FR_H1="6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
    elif len(FR_H1) == 21: new_number_FR_H1="5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
    elif len(FR_H1) == 22: new_number_FR_H1="4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
    elif len(FR_H1) == 23: new_number_FR_H1="3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
    elif len(FR_H1) == 24: new_number_FR_H1="2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
    elif len(FR_H1) == 25: new_number_FR_H1="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
    elif len(FR_H1) == 26:
        #new_number_FR_H1="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
        new_number_FR_H1="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25"
        FR_H1 = FR_H1[1:]  # Remove 0th residue 10/24/2012
    else: print "ERROR: FR_H1 length",len(FR_H1),"is not between 16 and 26"

    # H1
    if   len(H1) ==  6: new_number_H1="26,27,32,33,34,35"
    elif len(H1) ==  7: new_number_H1="26,27,28,32,33,34,35"
    elif len(H1) ==  8: new_number_H1="26,27,28,29,32,33,34,35"
    elif len(H1) ==  9: new_number_H1="26,27,28,29,30,32,33,34,35"
    elif len(H1) == 10: new_number_H1="26,27,28,29,30,31,32,33,34,35"
    elif len(H1) == 11: new_number_H1="26,27,28,29,30,31,31A,32,33,34,35"
    elif len(H1) == 12: new_number_H1="26,27,28,29,30,31,31A,31B,32,33,34,35"
    elif len(H1) == 13: new_number_H1="26,27,28,29,30,31,31A,31B,31C,32,33,34,35"
    else: print "ERROR: H1 length",len(H1),"is not between 6 and 13"

    # FR_H2
    if len(FR_H2) == 14: new_number_FR_H2="36,37,38,39,40,41,42,43,44,45,46,47,48,49"
    else: print "ERROR: FR_H2 length",len(FR_H2),"is not 14"

    # H2
    if   len(H2) == 12: new_number_H2="50,51,52,57,58,59,60,61,62,63,64,65"
    elif len(H2) == 13: new_number_H2="50,51,52,56,57,58,59,60,61,62,63,64,65"
    elif len(H2) == 14: new_number_H2="50,51,52,55,56,57,58,59,60,61,62,63,64,65"
    elif len(H2) == 15: new_number_H2="50,51,52,54,55,56,57,58,59,60,61,62,63,64,65"
    elif len(H2) == 16: new_number_H2="50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65"
    elif len(H2) == 17: new_number_H2="50,51,52,52A,53,54,55,56,57,58,59,60,61,62,63,64,65"
    elif len(H2) == 18: new_number_H2="50,51,52,52A,52B,53,54,55,56,57,58,59,60,61,62,63,64,65"
    elif len(H2) == 19: new_number_H2="50,51,52,52A,52B,52C,53,54,55,56,57,58,59,60,61,62,63,64,65"
    elif len(H2) == 20: new_number_H2="50,51,52,52A,52B,52C,52D,53,54,55,56,57,58,59,60,61,62,63,64,65"
    elif len(H2) == 21: new_number_H2="50,51,52,52A,52B,52C,52D,52E,53,54,55,56,57,58,59,60,61,62,63,64,65"
    elif len(H2) == 22: new_number_H2="50,51,52,52A,52B,52C,52D,52E,52F,53,54,55,56,57,58,59,60,61,62,63,64,65"
    else: print "ERROR: H2 length",len(H2),"is not between 12 and 22"

    # FR_H3
    if   len(FR_H3) == 30: new_number_FR_H3="66,67,68,69,70,71,72,73,76,77,78,79,80,81,82,82A,82B,82C,83,84,85,86,87,88,89,90,91,92,93,94"
    elif len(FR_H3) == 31: new_number_FR_H3="66,67,68,69,70,71,72,73,74,76,77,78,79,80,81,82,82A,82B,82C,83,84,85,86,87,88,89,90,91,92,93,94"
    elif len(FR_H3) == 32: new_number_FR_H3="66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,82A,82B,82C,83,84,85,86,87,88,89,90,91,92,93,94"
    else: print "ERROR: FR_H3 length",len(FR_H3),"is not between 30 and 32"

    # H3
    if   len(H3) ==  3: new_number_H3="95,96,97"
    elif len(H3) ==  4: new_number_H3="95,96,97,98"
    elif len(H3) ==  5: new_number_H3="95,96,97,98,99"
    elif len(H3) ==  6: new_number_H3="95,96,97,98,99,100"
    elif len(H3) ==  7: new_number_H3="95,96,97,98,99,101,102"
    elif len(H3) ==  8: new_number_H3="95,96,97,98,99,100,101,102"
    elif len(H3) ==  9: new_number_H3="95,96,97,98,99,100,100A,101,102"
    elif len(H3) == 10: new_number_H3="95,96,97,98,99,100,100A,100B,101,102"
    elif len(H3) == 11: new_number_H3="95,96,97,98,99,100,100A,100B,100C,101,102"
    elif len(H3) == 12: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,101,102"
    elif len(H3) == 13: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,101,102"
    elif len(H3) == 14: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,101,102"
    elif len(H3) == 15: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,101,102"
    elif len(H3) == 16: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,101,102"
    elif len(H3) == 17: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,101,102"
    elif len(H3) == 18: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,101,102"
    elif len(H3) == 19: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,101,102"
    elif len(H3) == 20: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,101,102"
    elif len(H3) == 21: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,101,102"
    elif len(H3) == 22: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,101,102"
    elif len(H3) == 23: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,101,102"
    elif len(H3) == 24: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,101,102"
    elif len(H3) == 25: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,101,102"
    elif len(H3) == 26: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,101,102"
    elif len(H3) == 27: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,101,102"
    elif len(H3) == 28: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,101,102"
    elif len(H3) == 29: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,100U,101,102"
    elif len(H3) == 30: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,100U,100V,101,102"
    elif len(H3) == 31: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,100U,100V,100W,101,102"
    elif len(H3) == 32: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,100U,100V,100W,100X,101,102"
    elif len(H3) == 33: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,100U,100V,100W,100X,100Y,101,102"
    elif len(H3) == 34: new_number_H3="95,96,97,98,99,100,100A,100B,100C,100D,100E,100F,100G,100H,100I,100J,100K,100L,100M,100N,100O,100P,100Q,100R,100S,100T,100U,100V,100W,100X,100Y,100Z,101,102"
    else: print "ERROR: H3 length",len(FR_H3),"is not between 3 and 34"

    # FR_H4
    new_number_FR_H4="103,104,105,106,107,108,109,110,111,112,113,114"

    try:
        (new_number_L1 and new_number_L2 and new_number_L3 and new_number_H1 and new_number_H2 and new_number_H3
         and new_number_FR_H1 and new_number_FR_H2 and new_number_FR_H3 and new_number_FR_H4
         and new_number_FR_L1 and new_number_FR_L2 and new_number_FR_L3 and new_number_FR_L4 )
    except:
        print "Numbering failed.  Exiting."
        sys.exit(1)

    # Converting all new_number_* vars in to a lists
    new_number_L1=new_number_L1.split(',');  new_number_L2=new_number_L2.split(',');  new_number_L3=new_number_L3.split(',');
    new_number_H1=new_number_H1.split(',');  new_number_H2=new_number_H2.split(',');  new_number_H3=new_number_H3.split(',');
    new_number_FR_L1=new_number_FR_L1.split(',');  new_number_FR_L2=new_number_FR_L2.split(',');  new_number_FR_L3=new_number_FR_L3.split(',');  new_number_FR_L4=new_number_FR_L4.split(',');
    new_number_FR_H1=new_number_FR_H1.split(',');  new_number_FR_H2=new_number_FR_H2.split(',');  new_number_FR_H3=new_number_FR_H3.split(',');  new_number_FR_H4=new_number_FR_H4.split(',');

    print_seq_FR_L1, print_seq_FR_L2, print_seq_FR_L3, print_seq_FR_L4, print_seq_FR_L4_extra = '', '', '', '', ''
    print_seq_FR_H1, print_seq_FR_H2, print_seq_FR_H3, print_seq_FR_H4, print_seq_FR_H4_extra = '', '', '', '', ''
    numbering_L, numbering_H = {}, {}

    # OUTPUT FOR LIGHT CHAIN. This should be save in the file 'numbering_L.txt'.
    for i, s in enumerate(FR_L1): #FR_L1
        numbering_L[ new_number_FR_L1[i] ] = s  #+='%s %s\n' % (s, new_number_FR_L1[i])
        v = int_(new_number_FR_L1[i])
        if v >= 10 and v <= 23: print_seq_FR_L1 += s

    for i, s in enumerate(L1): numbering_L[ new_number_L1[i] ] = s  # +='%s %s\n' % (s, new_number_L1[i])  # L1

    for i, s in enumerate(FR_L2):  #FR_L2
        numbering_L [ new_number_FR_L2[i] ] = s  #+='%s %s\n' % (s, new_number_FR_L2[i])
        v = int_(new_number_FR_L2[i])
        if (v >= 35 and v <= 38) or (v >= 45 and v <= 49): print_seq_FR_L2 += s

    for i, s in enumerate(L2): numbering_L[ new_number_L2[i] ] = s  #+='%s %s\n' % (s, new_number_L2[i])  # L2

    for i, s in enumerate(FR_L3):  #FR_L3
        numbering_L[ new_number_FR_L3[i] ] = s  # +='%s %s\n' % (s, new_number_FR_L3[i])
        v = int_(new_number_FR_L3[i])
        if (v >= 57 and v <= 66) or (v >= 71 and v <= 88): print_seq_FR_L3 += s

    for i, s in enumerate(L3): numbering_L[ new_number_L3[i] ] = s  #+='%s %s\n' % (s, new_number_L3[i])  # L3

    for i, s in enumerate(FR_L4):  #FR_L4
        numbering_L[ new_number_FR_L4[i] ] = s  #+='%s %s\n' % (s, new_number_FR_L4[i])
        v = int_(new_number_FR_L4[i])
        if v >= 98 and v <= 104: print_seq_FR_L4 += s
        if v >= 98 and v <= 101: print_seq_FR_L4_extra += s


    # OUTPUT FOR HEAVY CHAIN. This should be save in the file 'numbering_H.txt'.
    for i, s in enumerate(FR_H1):  #FR_H1
        numbering_H[ new_number_FR_H1[i] ] = s  #+='%s %s\n' % (s, new_number_FR_H1[i])
        v = int_(new_number_FR_H1[i])
        if v >= 10 and v <= 25: print_seq_FR_H1 += s

    for i, s in enumerate(H1): numbering_H[ new_number_H1[i] ] = s  #+='%s %s\n' % (s, new_number_H1[i])  # H1

    for i, s in enumerate(FR_H2):  #FR_H2
        numbering_H[ new_number_FR_H2[i] ] = s  #+='%s %s\n' % (s, new_number_FR_H2[i])
        v = int_(new_number_FR_H2[i])
        if (v >= 36 and v <= 39) or (v >= 46 and v <= 49): print_seq_FR_H2 += s

    for i, s in enumerate(H2): numbering_H[ new_number_H2[i] ] = s  #+='%s %s\n' % (s, new_number_H2[i])  # H2

    for i, s in enumerate(FR_H3):  #FR_H3
        numbering_H[ new_number_FR_H3[i] ] = s  #+='%s %s\n' % (s, new_number_FR_H3[i])
        v = int_(new_number_FR_H3[i])
        if v >= 66 and v <= 94: print_seq_FR_H3 += s

    for i, s in enumerate(H3): numbering_H[ new_number_H3[i] ] = s  #+='%s %s\n' % (s, new_number_H3[i])  # H3

    for i, s in enumerate(FR_H4):  #FR_H4
        numbering_H[ new_number_FR_H4[i] ] = s  #+='%s %s\n' % (s, new_number_FR_H4[i])
        v = int_(new_number_FR_H4[i])
        if v >= 103 and v <= 109: print_seq_FR_H4 += s
        if v >= 103 and v <= 106: print_seq_FR_H4_extra += s

    #if Options.verbose: print 'numbering_L:\n', numbering_L
    #if Options.verbose: print 'numbering_H:\n', numbering_H

    FRL = print_seq_FR_L1 + print_seq_FR_L2 + print_seq_FR_L3 + print_seq_FR_L4
    FRH = print_seq_FR_H1 + print_seq_FR_H2 + print_seq_FR_H3 + print_seq_FR_H4

    if len(FRL) != 58 and len(FRL) != 60:
        print "ERROR: Current DB does not cover the length of FRL of your query."
        print "ERROR: FRL length of your query:", len(FRL)
        print "ERROR: DB: 61 or 63"
        sys.exit(1)

    if len(FRH) != 63 and len(FRH) != 65:
        print "ERORR: Current DB does not cover the length of FRL of your query."
        print "ERORR: FRH length of your query:", len(FRH)
        print "ERORR: DB: 63 or 65"
        sys.exit(1)

    return dict(FRL = print_seq_FR_L1 + print_seq_FR_L2 + print_seq_FR_L3 + print_seq_FR_L4,
                FRH = print_seq_FR_H1 + print_seq_FR_H2 + print_seq_FR_H3 + print_seq_FR_H4,
                light = FR_L1 + L1 + FR_L2 + L2 + FR_L3 + L3 + print_seq_FR_L4_extra,
                heavy = FR_H1 + H1 + FR_H2 + H2 + FR_H3 + H3 + print_seq_FR_H4_extra,
                numbering_L=numbering_L, numbering_H=numbering_H)


def run_blast(cdr_query, prefix, blast, blast_database, verbose=False):
    print '\nRunning %s' % (blast)
    cdr_info, legend = {}, ''  # first reading cdr_info table
    for l in file( _script_path_ + '/info/cdr_info' ):
        if l.startswith('# '): legend = l[2:].split()
        elif len(l)>8: cdr_info[l.split()[0]] =  dict( zip(legend, l.split() ) )


    for line in file(  _script_path_ + '/info/fv_length' ): cdr_info[ line.split()[0] ]  [ line.split()[1]+'_length' ]  =  int( line.split()[2] )

    # cdr_info consistency check
    for name, i in cdr_info.items():
        if sum( [ i[k]!='none' and int(i[k+'_length'])!= len(i[k]) for k in ['L1', 'L2', 'L3', 'H1', 'H2', 'H3'] ]):
            print 'ERROR: cdr_info length info is inconsistent for line: %s' % name;
            for k in ['L1' , 'L2' , 'L3' , 'H1' , 'H2' , 'H3']:
                print k, i[k], i[k+'_length'], len(i[k])
            sys.exit(1)

    alignment = {'summary':[]}
    for k in _framework_names_:
        input_query = k + '.fasta'
        output = k + '.align'  # Options.prefix +

		# if CDRs, use length-depend DBs for BLAST
        len_cdr = len(cdr_query[k])
        if k.count('FR') or k.count('heavy') or k.count('light'): db = blast_database + '/database.%s' % k
        else: db = blast_database + '/database.%s.%s' % (k, len_cdr)

        # check that database file exists
        if (not os.path.isfile(db)):
            print '\nERROR: No %s templates of length %s (%s)\n' % (k, len_cdr, db)
            sys.exit(1)

        wordsize, e_value, matrix = (2, 0.00001, 'BLOSUM62') if k.count('FR') or k.count('heavy') or k.count('light') else (2, 2000, 'PAM30')

        commandline = 'cd %s && ' % '\ '.join(os.path.dirname(prefix).split()) + blast + \
                      (' -db %(db)s -query %(input_query)s -out %(output)s -evalue %(e_value)s -matrix %(matrix)s' + \
                       ' -word_size %(wordsize)s -outfmt 7 -max_target_seqs 600') % vars()

        res, output = commands.getstatusoutput(commandline)
        if verbose and not res: print commandline+'\n', output
        if res: print commandline+'\n', 'ERROR: Blast execution faild with code %s and message: %s' % (res, output);  sys.exit(1)

        #print 'Filtering results...'
        table, legend = [], ''
        for l in file(prefix +k + '.align'):
            if l.startswith('# Fields: '): legend = [ i.strip().replace('. ', '.').replace(' ', '-') for i in l[10:-1].split(',')]
            elif l.startswith('Query_1'): table.append( dict( zip(legend, l.split() ) ) )

        table_original = table[:]

        def sort_and_write_results(table, file_name):
            table.sort(key=lambda x: (-float(x['bit-score']), float(cdr_info[ x['subject-id'] ]['resolution'] ) ) )  # sort results first by bit-score then by resolution
            for r in table: r['resolution'] = cdr_info[ r['subject-id'] ]['resolution']
            with file(file_name, 'w') as f:
                f.write('# Fields: ' + ' '.join(legend) + '\n')
                for i in table: f.write( '\t'.join([i[k] for k in legend]) + '\n' )

                def try_to_float(v):
                    try: return int(v)
                    except (ValueError):
                        try: return float(v)
                        except (ValueError): return v

                for i, v in enumerate(table):
                    table[i] = dict( [(k, try_to_float(v[k])) for k in v] )


        for f in Filters:
            if getattr(Options, f.func_name):
                f(k, table, cdr_query, cdr_info)
                t = table_original[:];  f(k, t, cdr_query, cdr_info)
                sort_and_write_results([i for i in table_original if i not in t], prefix+'filtered-by.'+f.func_name[7:]+'.'+k+'.align')

        sort_and_write_results(table, prefix + 'filtered.' + k + '.align')  # Writing filtered results.
        alignment[k] = table;
        if table: alignment['summary'].append(dict(table[0]));  alignment['summary'][-1]['subject-id']=k

        custom_template = getattr(Options, k)

        if custom_template and not os.path.isfile(custom_template): custom_template = '/pdb%s_chothia.pdb' % custom_template
        if not custom_template:
            if table:
                custom_template = table[0]['subject-id']
            else:  # if there is no template... table is a list, which has a blast result
                for v in cdr_info.items():
                    check_length = '%s_length' % k
                    if len_cdr == int(v[1][check_length]):
                        pdb_random = v[0]
                        break
                print '\nWARNING: No template avaliable for %s after filtering! Using a random template of the same length as the query\n' % k
                custom_template = pdb_random
                #sys.exit(1)
            print "%s template: %s" % (k, custom_template)
        else: print 'Custom %s template: %s...' % (k, custom_template)
        shutil.copy(Options.antibody_database+'/'+custom_template, prefix+'/template.'+k+'.pdb')

    legend.remove('query-id'); legend.insert(1, 'resolution');  return alignment, legend

def create_virtual_template_pdbs(prefix):
    # Make a template PDB file, which has psuedo atoms, so that the query and a template can have the same length sequence.
    for chain in [ 'L', 'H' ]:
        with file(prefix+'/template.tmp.FR'+chain+'.pdb', 'w') as o:
            # Count the number of lines of a template PDB file
            cnt1=0
            for line2 in file(prefix+'/template.FR'+chain+'.pdb'):
                if line2[0:4] == 'ATOM':
                    chain_temp = line2[21:22]
                    if chain_temp == chain: cnt1+=1

            cnt_res=0
            for line in file(prefix+'/numbering_'+chain+'.txt'):
                cnt_res+=1
                check = 0
                res_num_q = line[2:].rstrip('\n')
                #res_aa    = line[:1]
                cnt2=0
                for line2 in file(prefix+'/template.FR'+chain+'.pdb'):
                    if line2[0:4] == 'ATOM':
                        chain_temp = line2[21:22]
                        res_aa  = line2[17:20]
                        if chain_temp == chain: cnt2+=1
                        res_num_temp = line2[23:27].strip()
                        iCode_tmp   = line2[27:27]
                        x_coord = line2[31:38]
                        y_coord = line2[39:46]
                        z_coord = line2[47:54]

                        # Identify the residue number of the last residue
                        if cnt2 == 1 and chain_temp == chain:
                            res_num_temp_first = res_num_temp
                        elif cnt2 == cnt1 and chain_temp == chain:
                            res_num_temp_last = res_num_temp

                        if chain == chain_temp and res_num_q == res_num_temp:
                            check = 1
                            o.write( line2 )

                x_coord_n  = '%8.3f' % ( float(x_coord) + 100 + cnt_res * 2)
                x_coord_ca = '%8.3f' % ( float(x_coord) + 110 + cnt_res * 2)
                x_coord_c  = '%8.3f' % ( float(x_coord) + 120 + cnt_res * 2)
                x_coord_o  = '%8.3f' % ( float(x_coord) + 130 + cnt_res * 2)

                y_coord_n  = '%8.3f' % ( float(y_coord) + 100 + cnt_res * 2)
                y_coord_ca = '%8.3f' % ( float(y_coord) + 110 + cnt_res * 2)
                y_coord_c  = '%8.3f' % ( float(y_coord) + 120 + cnt_res * 2)
                y_coord_o  = '%8.3f' % ( float(y_coord) + 130 + cnt_res * 2)

                z_coord_n  = '%8.3f' % ( float(z_coord) + 100 + cnt_res * 2)
                z_coord_ca = '%8.3f' % ( float(z_coord) + 110 + cnt_res * 2)
                z_coord_c  = '%8.3f' % ( float(z_coord) + 120 + cnt_res * 2)
                z_coord_o  = '%8.3f' % ( float(z_coord) + 130 + cnt_res * 2)

                #x_coord_n  = '%8.3f' % ( float(x_coord) + random.uniform(1, 200) )
                #x_coord_ca = '%8.3f' % ( float(x_coord) + random.uniform(1, 200) )
                #x_coord_c  = '%8.3f' % ( float(x_coord) + random.uniform(1, 200) )
                #x_coord_o  = '%8.3f' % ( float(x_coord) + random.uniform(1, 200) )

                #y_coord_n  = '%8.3f' % ( float(y_coord) + random.uniform(1, 200) )
                #y_coord_ca = '%8.3f' % ( float(y_coord) + random.uniform(1, 200) )
                #y_coord_c  = '%8.3f' % ( float(y_coord) + random.uniform(1, 200) )
                #y_coord_o  = '%8.3f' % ( float(y_coord) + random.uniform(1, 200) )

                #z_coord_n  = '%8.3f' % ( float(z_coord) + random.uniform(1, 200) )
                #z_coord_ca = '%8.3f' % ( float(z_coord) + random.uniform(1, 200) )
                #z_coord_c  = '%8.3f' % ( float(z_coord) + random.uniform(1, 200) )
                #z_coord_o  = '%8.3f' % ( float(z_coord) + random.uniform(1, 200) )

                if res_num_q.isdigit():
                    if int(res_num_q) < 100:
                        res_num_q_tmp = '  '+res_num_q
                    elif int(res_num_q) >= 100:
                        res_num_q_tmp = ' '+res_num_q
                    iCode_q_tmp   = ' '
                else:
                    res_num_q_tmp = res_num_q[:-1]
                    if int(res_num_q_tmp) < 100:
                        res_num_q_tmp = '  '+res_num_q_tmp
                    elif int(res_num_q_tmp) >= 100:
                        res_num_q_tmp = ' '+res_num_q_tmp
                    iCode_q_tmp   = res_num_q[-1:]

                #if Options.verbose: print chain, res_num_q, res_num_q_tmp, check, res_num_temp_first, res_num_temp_last
                if check == 0 and ( int(res_num_temp_first) <= int(res_num_q_tmp) <= int(res_num_temp_last) ):
                    r_n  = dict(tempFactor=' 25.00', chainID=chain, name=' N  ', altLoc=' ', occupancy='  1.00', element='  ', resSeq=res_num_q_tmp, charge='  ', iCode=iCode_q_tmp, resName='GLY', x=x_coord_n,  serial=' 0000', z=z_coord_n,  type='ATOM  ', y=y_coord_n)
                    r_ca = dict(tempFactor=' 25.00', chainID=chain, name=' CA ', altLoc=' ', occupancy='  1.00', element='  ', resSeq=res_num_q_tmp, charge='  ', iCode=iCode_q_tmp, resName='GLY', x=x_coord_ca, serial=' 0000', z=z_coord_ca, type='ATOM  ', y=y_coord_ca)
                    r_c  = dict(tempFactor=' 25.00', chainID=chain, name=' C  ', altLoc=' ', occupancy='  1.00', element='  ', resSeq=res_num_q_tmp, charge='  ', iCode=iCode_q_tmp, resName='GLY', x=x_coord_c,  serial=' 0000', z=z_coord_c,  type='ATOM  ', y=y_coord_c)
                    r_o  = dict(tempFactor=' 25.00', chainID=chain, name=' O  ', altLoc=' ', occupancy='  1.00', element='  ', resSeq=res_num_q_tmp, charge='  ', iCode=iCode_q_tmp, resName='GLY', x=x_coord_o,  serial=' 0000', z=z_coord_o,  type='ATOM  ', y=y_coord_o)

                    o.write( records_to_pdb_string(r_n)  + '\n' );
                    o.write( records_to_pdb_string(r_ca) + '\n' );
                    o.write( records_to_pdb_string(r_c)  + '\n' );
                    o.write( records_to_pdb_string(r_o)  + '\n' );

def thread_template_pdbs(CDRs, prefix):
    L_regions=dict(L1=(24,34), L2=(50,56), L3=(89,97));  H_regions=dict(H1=(26,35), H2=(50,65), H3=(95, 102) )
    for chain, k, numbering, regions in [ ( 'L', 'FRL', 'numbering_L', L_regions ), ( 'H', 'FRH', 'numbering_H', H_regions ) ]:
        with file(prefix+'/template.threaded.'+k+'.pdb', 'w') as o:
            for line in file(prefix+'/template.tmp.'+k+'.pdb'): # This is a template PDB in our database
            #for line in file(prefix+'/template.'+k+'.pdb'): # This is a template PDB in our database
                r = map_pdb_string_to_records(line)         # 'r' is a info of a template
                if r['type'] == "ATOM  ":
                    res_num = int(r['resSeq']);  res_num_icode = ( '%s%s' % (res_num, r['iCode']) ).strip()
                    if res_num_icode in CDRs[numbering]  and chain == r['chainID']:
                        for loop in regions:
                            if regions[loop][0] <= res_num <= regions[loop][1]:
                                #r['charge'] = '.' if r['resName'] == AA_Code[ CDRs[numbering][res_num_icode] ] else '!'
                                if r['name'] in [' N  ', ' CA ', ' C  ', ' O  ']:  # r['resName'] == AA_Code[ CDRs[numbering][res_num_icode] ]
                                    for coord in 'xyz':
                                        #print r[coord]
                                        r[coord] = '%8.3f' % ( float(r[coord]) + 100.0 )
                                    r['resName'] = AA_Code[ CDRs[numbering][res_num_icode] ]
                                    o.write( records_to_pdb_string(r) + '\n' );  break
                        else:
                            if r['resName'] == AA_Code[ CDRs[numbering][res_num_icode] ] or r['name'] in [' N  ', ' CA ', ' C  ', ' O  ']:
                                r['resName'] = AA_Code[ CDRs[numbering][res_num_icode] ]
                                o.write( records_to_pdb_string(r) + '\n' )

    for chain, numbering, regions in [ ( 'L', 'numbering_L', L_regions ), ( 'H', 'numbering_H', H_regions ) ]:
        for R in regions:
            with file(prefix+R+'.pdb', 'w') as o:
                for line in file(prefix+'/template.'+R+'.pdb'):
                    r = map_pdb_string_to_records(line)
                    if r['type'] == "ATOM  ":
                        res_num = int(r['resSeq']);  res_num_icode = ( '%s%s' % (res_num, r['iCode']) ).strip()
                        if res_num_icode in CDRs[numbering]  and chain == r['chainID']:
                            if regions[R][0]-4 <= res_num <= regions[R][1]+4:
                                if r['name'] in [' N  ', ' CA ', ' C  ', ' O  ']  or  r['resName'] == AA_Code[ CDRs[numbering][res_num_icode] ]:
                                    r['resName'] = AA_Code[ CDRs[numbering][res_num_icode] ]
                                    o.write( records_to_pdb_string(r) + '\n' );



profit_templates = { 'L': '''reference "%(prefix)s/template.light_heavy.pdb"
mobile "%(prefix)s/template.threaded.FRL.pdb"
atoms ca
zone L10-L23
zone L35-L49
zone L57-L66
zone L69-L88
zone L98-L100
fit
write "%(prefix)s/fitted.L.pdb"
quit''',
'H': '''reference "%(prefix)s/template.light_heavy.pdb"
mobile "%(prefix)s/template.threaded.FRH.pdb"
atoms ca
zone H10-H25
zone H36-H49
zone H66-H72
zone H74-H82
zone H83-H94
zone H103-H105
fit
write "%(prefix)s/fitted.H.pdb"
quit''' }

def superimpose_templates(CDRs, prefix):
    for chain in profit_templates:
        f_name = prefix + 'profit-%s' % chain
        with file(f_name+'.in', 'w') as f: f.write(profit_templates[chain] % vars() )

        f_name = '\ '.join(f_name.split())
        commandline = '%s < %s.in > %s.out' % (Options.profit, f_name, f_name)
        res, output = commands.getstatusoutput(commandline);
        if res: print commandline, output; sys.exit(1)

    pathPrefix = '\ '.join(vars()['prefix'].split())
    pathPrefix = pathPrefix[:-1] if pathPrefix.endswith('/') else pathPrefix
    res, output = commands.getstatusoutput('cat {0}/fitted.L.pdb {0}/fitted.H.pdb > {0}/FR.pdb'.format(pathPrefix))

    #res, output = commands.getstatusoutput('cat %(prefix)s/fitted.L.pdb %(prefix)s/fitted.H.pdb > %(prefix)s/FR.pdb' % vars())
    if res: print output;  sys.exit(1)


def run_rosetta(CDRs, prefix, rosetta_bin, rosetta_platform, rosetta_database):
    antibody_graft = rosetta_bin + '/antibody_graft.' + rosetta_platform
    if os.path.isfile( antibody_graft ):
        print '\nRunning antibody_graft'
        commandline = 'cd "%s/details" && "%s" -database %s -overwrite -s FR.pdb' % (os.path.dirname(prefix), antibody_graft, rosetta_database) + \
                      ' -antibody::graft_l1 -antibody::graft_l2 -antibody::graft_l3' + \
                      ' -antibody::graft_h1 -antibody::graft_h2 -antibody::graft_h3' + \
                      ' -antibody::h3_no_stem_graft'
        if Options.quick: commandline = commandline + ' -run:benchmark -antibody:stem_optimize false'
        res, output = commands.getstatusoutput(commandline)
        if Options.verbose or res: print commandline, output
        if res: print 'Rosetta run terminated with Error!'; sys.exit(1)
        model_file_prefix = 'grafted';  shutil.move(prefix+'details/FR_0001.pdb', prefix+model_file_prefix+'.pdb')
    else:
        print 'Rosetta executable %s was not found, skipping Rosetta run...' % antibody_graft
        return

    if Options.idealize:
        idealize_jd2 = rosetta_bin + '/idealize_jd2.' + rosetta_platform
        if os.path.isfile( idealize_jd2 ):
            print 'Running idealize_jd2'
            commandline = 'cd "%s" && "%s" -database %s -overwrite' % (os.path.dirname(prefix), idealize_jd2, rosetta_database) + \
                          ' -fast -s %s.pdb -ignore_unrecognized_res' % model_file_prefix
            res, output = commands.getstatusoutput(commandline)
            if Options.verbose or res: print commandline, output
            if res: print 'Rosetta run terminated with Error!  Commandline: %s' % commandline; sys.exit(1)
            shutil.move(prefix + model_file_prefix + '_0001.pdb', prefix + model_file_prefix + '.idealized.pdb');  model_file_prefix += '.idealized'
        else:
            print 'Rosetta executable %s was not found, skipping Rosetta run...' % idealize_jd2
            return

    if Options.relax:
        relax = rosetta_bin + '/relax.' + rosetta_platform
        if os.path.isfile( relax ):
            print 'Running relax with all-atom constraint'
            commandline = 'cd "%s" && "%s" -database %s -overwrite' % (os.path.dirname(prefix), relax, rosetta_database) + \
                          ' -s %s.pdb -ignore_unrecognized_res -relax:fast -relax:constrain_relax_to_start_coords' % model_file_prefix + \
                          ' -relax:coord_constrain_sidechains -relax:ramp_constraints false -ex1 -ex2 -use_input_sc'
            res, output = commands.getstatusoutput(commandline)
            if Options.verbose or res: print commandline, output
            if res: print 'Rosetta run terminated with Error!  Commandline: %s' % commandline; sys.exit(1)
            shutil.move(prefix + model_file_prefix + '_0001.pdb', prefix + model_file_prefix + '.relaxed.pdb');  model_file_prefix += '.relaxed'
        else:
            print 'Rosetta executable %s was not found, skipping Rosetta run...' % relax
            return

    shutil.copy(prefix + model_file_prefix + '.pdb', prefix+'model.pdb')
    make_cter_constraint(CDRs, prefix)



# Dihedral CA 220 CA 221 CA 222 CA 223 SQUARE_WELL2 0.523 0.698 200; KINK
# Dihedral CA 220 CA 221 CA 222 CA 223 SQUARE_WELL2 2.704 0.523 100; EXTEND
def make_cter_constraint(CDRs, prefix):
    print '\nPreparing cter_constraint file for H3 modeling'
    L36  = AA_Code[ CDRs['numbering_L']['36'] ]
    L46  = AA_Code[ CDRs['numbering_L']['46'] ]
    L49  = AA_Code[ CDRs['numbering_L']['49'] ]
    H93  = AA_Code[ CDRs['numbering_H']['93'] ]
    H94  = AA_Code[ CDRs['numbering_H']['94'] ]
    H99  = AA_Code[ CDRs['H3'][-4:-3] ]  # n-3
    H100 = AA_Code[ CDRs['H3'][-3:-2] ]  # n-2
    H101 = AA_Code[ CDRs['H3'][-2:-1] ]  # n-1
    len_h3 = len( CDRs['H3'] )
    print 'H3:', len_h3, CDRs['H3'],'\nKey residues: ', L36,L46,L49,H93,H94,H99,H100,H101

    # H3-rules
    if ( H93 == 'ARG' or H93 == 'LYS' ) and ( H94 == 'ARG' or H94 == 'LYS' ) and H101 == 'ASP': base = 'KINK'
    elif ( H94 == 'ARG' or H94 == 'LYS' ) and H101 == 'ASP':
        if ( L46 == 'ARG' or L46 == 'LYS' ) and L36 != 'TYR': base = 'EXTEND'
        else: base = 'KINK'
    elif ( H93 == 'ARG' or H93 == 'LYS' ) and H101 == 'ASP': base = 'KINK'
    elif H101 == 'ASP':
        if L49 == 'ARG' or L49 == 'LYS': base = 'KINK'
        elif (H100 == 'MET' or H100 == 'PHE') and ( H99 == 'ALA' or H99 == 'GLY' ): base = 'KINK'
        else: base = 'EXTEND'
    elif (H100 == 'MET' or H100 == 'PHE') and ( H99 == 'ALA' or H99 == 'GLY' ): base = 'KINK'
    elif H100 == 'ASP' or H100 == 'ASN' or H100 == 'LYS' or H100 == 'ARG': base = 'EXTEND'
    elif len_h3 == 7: base = 'EXTEND'
    else: base = 'KINK'

    cnt=0
    f=open(prefix+'cter_constraint', 'w')
    for line in file(prefix+'/model.pdb'):
        if line[0:4] == 'ATOM':
            chain = line[21:22]
            atom  = line[13:15]
            res_num = line[22:27].strip()

            if atom == 'CA':cnt+=1
            if chain == 'H' and atom == 'CA' and res_num == str(103):
                n1=cnt-3 # n-2 (H100X)
                n2=cnt-2 # n-1 (H101)
                n3=cnt-1 # n   (H102)
                n4=cnt   # n+1 (H103)

                if base == 'KINK': f.write( 'Dihedral CA '+str(n1)+' CA '+str(n2)+' CA '+str(n3)+' CA '+str(n4)+' SQUARE_WELL2 0.523 0.698 200' )
                elif base == 'EXTEND': f.write( 'Dihedral CA '+str(n1)+' CA '+str(n2)+' CA '+str(n3)+' CA '+str(n4)+' SQUARE_WELL2 2.704 0.523 100' )
    f.close()
    print 'Predicted base conformation is', base

# Various filter function
def filter_by_sequence_homologue(k, results, cdr_query, cdr_info):
    if Options.verbose: print 'filtering by sequence identity...'
    #print results

    for r in results[:]:
        pdb = r['subject-id']

        if k in ['L1','L2','L3','H1','H2','H3'] and sid_checker(cdr_query[k], cdr_info[pdb][k]) >= sid_cutoff_cdr:
            results.remove(r)
            if Options.verbose and r not in results: print 'Filter sequence_identity (%s%% cut-off), removing:%s %s_query:%s %s_info:%s %s%% identity' % (sid_cutoff_cdr, pdb, k, len(cdr_query[k]), k, len( cdr_info[pdb][k] ), round(sid_checker(cdr_query[k], cdr_info[pdb][k]), 2) )

        elif  k in ['FRL'] and sid_checker(cdr_query[k], frl_info[pdb][k]) >= sid_cutoff_fr:
            results.remove(r)
            if Options.verbose and r not in results: print 'Filter sequence_identity (%s%% cut-off), removing:%s %s_query:%s %s_info:%s %s%% identity' % (sid_cutoff_fr, pdb, k, len(cdr_query[k]), k, len( frl_info[pdb][k] ), round(sid_checker(cdr_query[k], frl_info[pdb][k]), 2) )

        elif  k in ['FRH'] and sid_checker(cdr_query[k], frh_info[pdb][k]) >= sid_cutoff_fr:
            results.remove(r)
            if Options.verbose and r not in results: print 'Filter sequence_identity (%s%% cut-off), removing:%s %s_query:%s %s_info:%s %s%% identity' % (sid_cutoff_fr, pdb, k, len(cdr_query[k]), k, len( frh_info[pdb][k] ), round(sid_checker(cdr_query[k], frh_info[pdb][k]), 2) )
        elif  k in ['light_heavy'] and float(r['%-identity'])  >= sid_cutoff_fr - 5:
            results.remove(r)
            if Options.verbose and r not in results: print 'Filter sequence_identity (%s%% cut-off), removing:%s %s %s%% identity' % (sid_cutoff_fr, pdb, k, r['%-identity'])


def sid_checker(seq_q, seq_t):
    seq_lenq  = len(seq_q)
    seq_lent = len(seq_t)

    #print seq_lenq, seq_lent

    if seq_lenq == seq_lent:
        num_res = 0
        for i in range(0, seq_lenq):
            #print seq_q[i], seq_t[i], cnt
            if seq_q[i] == seq_t[i]:
                num_res+=1
        ratio = float(num_res) / seq_lenq * 100
    else:
        ratio = 0

    return ratio


def filter_by_sequence_length(k, results, cdr_query, cdr_info):
    if Options.verbose: print 'filtering by sequence length...'

    for r in results[:]:
        pdb = r['subject-id']
        if   k in ['heavy','light_heavy']: return
        elif k in ['L1','L2','L3','H1','H2','H3'] and  len(cdr_query[k]) != len( cdr_info[pdb][k] ):
            results.remove(r)
            if Options.verbose and r not in results: print 'Filter sequence_length, removing:%s %s_query:%s %s_info:%s' % (pdb, k, len(cdr_query[k]), k, len( cdr_info[pdb][k] ))
        elif k == 'light' and 'light_lenght' in cdr_info[pdb]  and  len(cdr_query[k]) != int( cdr_info[pdb]['light_lenght'] ): results.remove(r)
        elif k == 'FRH':
            template_length = 61 if pdb == 'pdb2x7l_chothia.pdb' else 63
            if not len(cdr_query[k]) == template_length: results.remove(r)
            if Options.verbose and r not in results: print 'Filter sequence_length, removing:%s %s_query:%s %s_info:%s' % (pdb, k, len(cdr_query[k]), k, template_length)
            #elif k == 'FRH'  and  (not  len(cdr_query[k])-8 <= 67 <= len(cdr_query[k])+8): results.remove(r)
        elif k == 'FRL':
            template_length = 60 if pdb == 'pdb3h0t_chothia.pdb' else 58
            if not len(cdr_query[k]) == template_length: results.remove(r)
            if Options.verbose and r not in results: print 'Filter sequence_length, removing:%s %s_query:%s %s_info:%s' % (pdb, k, len(cdr_query[k]), k, template_length)
            #if not  len(cdr_query[k])-8 <= template_length <= len(cdr_query[k])+8: results.remove(r)


def filter_by_alignment_length(k, results, cdr_query, cdr_info):
    for r in results[:]:
        pdb = r['subject-id']
        if k == 'H3' and int(r['alignment-length']) < 0.10 * len( cdr_info[pdb]['H3'] ): results.remove(r)
        if k == 'H2' and int(r['alignment-length']) < 0.55 * len( cdr_info[pdb]['H2'] ): results.remove(r)
        if k in ['L1', 'L2', 'L3', 'H1']  and  int(r['alignment-length']) < 0.70 *  len( cdr_info[pdb][k] ): results.remove(r)
        if Options.verbose and r not in results: print 'Filter alignment_length removing:%s %s_query:%s alignment-length:%s ' % (pdb, k, len(cdr_info[pdb][k]), r['alignment-length'])


def filter_by_template_resolution(k, results, cdr_query, cdr_info):
    for r in results[:]:
        pdb = r['subject-id']
        if float(cdr_info[pdb]['resolution']) > 2.8: results.remove(r)
        if Options.verbose and r not in results: print 'Filter template_resolution, removing:%s resolution:%s' % (pdb, cdr_info[pdb]['resolution'])


def filter_by_template_bfactor(k, results, cdr_query, cdr_info):
    bfactor = {}

    if k in ['L1', 'L2', 'L3', 'H1', 'H2', 'H3']:
        for line in file( _script_path_ + '/info/list_bfactor50' ):
            for i, e in enumerate(['L1', 'L2', 'L3', 'H1', 'H2', 'H3']):
                if k == e: bfactor[ 'pdb'+line.split()[0]+'_chothia.pdb' ] = line.split()[i+1] == 'True'

        for r in results[:]:
            pdb = r['subject-id']
            if bfactor.get( pdb, False ): results.remove(r)
            if Options.verbose and r not in results: print 'Filter B-factor50, removing:%s' % pdb


def filter_by_outlier(k, results, cdr_query, cdr_info):
    outlier = {}
    for line in file( _script_path_ + '/info/outlier_list' ): outlier[tuple(line.split()[:2])] = line.split()[2] == 'true'
    for r in results[:]:
        pdb = r['subject-id']

        if outlier.get( (pdb, k), False ): results.remove(r)
        if Options.verbose and r not in results: print 'Filter outlier, removing:%s' % pdb


AA_Code = dict(A='ALA', V='VAL', L='LEU', I='ILE', P='PRO', W='TRP', F='PHE', M='MET', G='GLY', S='SER', T='THR', Y='TYR',
               C='CYS', N='ASN', Q='GLN', K='LYS', R='ARG', H='HIS', D='ASP', E='GLU')

# PDB Reader code
# PDB record maps, only ATOMS here for now
PDB_Records = {
    "ATOM  " : {
        "type"      : (1,  6),
        "serial"    : (7, 11),  # Integer
        "name"      : (13, 16), # Atom
        "altLoc"    : (17, 17), # Character
        "resName"   : (18, 20), # Residue name
		"chainID"   : (22, 22), # Character
		"resSeq"    : (23, 26), # Integer
		"iCode"     : (27, 27), # AChar
		"x"         : (31, 38), # Real(8.3)
		"y"         : (39, 46), # Real(8.3)
		"z"         : (47, 54), # Real(8.3)
		"occupancy" : (55, 60), # Real(6.2)
		"tempFactor": (61, 66), # Real(6.2)
		#///"segID",     Field(73, 76),
		"element"   : (77, 78), # LString(2)
		"charge"    : (79, 80)  # LString(2)
    },
    "UNKNOW" : {
        "type" : (1,  6),
        "info" : (7, 80)
    }
}

def map_pdb_string_to_records(s):
    rtype = s[:6]
    F = PDB_Records.get(rtype, PDB_Records["UNKNOW"])
    R = {}
    for f in F: R[f] = s[ F[f][0]-1: F[f][1] ]
    return R

def records_to_pdb_string(records):
    res = [' ']*80
    for f in records: res[ PDB_Records['ATOM  '][f][0]-1: PDB_Records['ATOM  '][f][1] ] = records[f]
    return ''.join(res)

def self_test():
    if os.path.isdir(Options.self_test_dir): print 'Removing old self-test-dir %s...' % Options.self_test_dir;  shutil.rmtree(Options.self_test_dir)
    os.makedirs(Options.self_test_dir)  # if not os.path.isdir( Options.self_test_dir ):

    tests = [d for d in os.listdir('test') if d != '.svn']
    print 'Checking [%s] targets: %s...' % (len(tests), tests)
    for t in tests:
        test_dir = Options.self_test_dir+t+'/';  os.makedirs(test_dir)
        if Options.verbose: print 'Testing target: %s...' % t
        light_chain = read_fasta_file('test/%s/query_l.fasta' % t)
        heavy_chain = read_fasta_file('test/%s/query_h.fasta' % t)
        answers = json.load( file('test/%s/%s.json' % (t, t) ) )
        answers['numbering_L'] = file('test/%s/numbering_L.txt' % t).read()
        answers['numbering_H'] = file('test/%s/numbering_H.txt' % t).read()

        if Options.verbose:
            print 'light_chain:', light_chain
            print 'heavy_chain:', heavy_chain

        CDRs = IdentifyCDRs(light_chain, heavy_chain)
        if Options.verbose: print 'CDR:', json.dumps(CDRs, sort_keys=True, indent=2)
        CDRs.update( Extract_FR_CDR_Sequences(**CDRs) )

        write_results(CDRs, prefix=test_dir)
        CDRs['numbering_L'] = file(test_dir+'/numbering_L.txt').read()
        CDRs['numbering_H'] = file(test_dir+'/numbering_H.txt').read()

        for a in sorted( answers.keys() ):  #['L1', 'L2', 'L3', 'H1', 'H2', 'H3']: #
            if answers[a] != CDRs[a]:
                print 'ERROR: target=%s field %s is not equal!!!\nexpected:%s\n     got:%s' % (t, a, answers[a], CDRs[a])
                sys.exit(1)
            elif Options.verbose: print '  %s: OK' % a



        '''
        run_blast(CDRs, prefix=test_dir, blast=Options.blast, blast_database=Options.blast_database)

        for a in ['FRL', 'FRH', 'L1', 'L2', 'L3', 'H1', 'H2', 'H3', 'light_heavy']:
            output_file = test_dir + a + '.align'
            answer_file = 'test/%s/%s.align' % (t, a)

            o_lines = [l for l in file(output_file) if not l.startswith('# ')]
            a_lines = [l for l in file(answer_file) if not l.startswith('# ')]

            if a_lines != o_lines:
                print 'ERROR: target=%s , Align files:%s and %s is not equal!!!' % (t, answer_file, output_file)
                sys.exit(1)
        '''

        if Options.verbose: print 'Testing target: %s... OK\n' % t


    print 'All tests passed!'
    return


if __name__ == "__main__": main(sys.argv)


'''def _convert_test_files():
    T = {}
    for file_name in ['CDR_ANSWERS', 'FRL_FRH_ANSWERS', 'LIGHT_HEAVY_ANSWERS']:
        lines = [l for l in file(file_name).read().split('\n') if l]
        legend = lines[0].split();  legend[0] = 'name'
        for l in lines[1:]:
            fields = l.split();  target = fields[0]
            for i, f in enumerate(fields[1:]):
                if target not in T: T[target] = {}
                T[target] [ legend[i+1] ] = f

    for t in T:
        f = file('test/%s/%s.json' % (t, t), 'w');  f.write( json.dumps(T[t], sort_keys=True, indent=2) );  f.close()
'''