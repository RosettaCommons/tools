#!/usr/bin/perl
###############################################################################
###############################################################################
# The purpose of this script is to generate input files required for the
# rosetta fragment picker and to run the picker to generate fragment libraries.
#
# This script takes a fasta file and runs the following:
#
#  1. psiblast to generate a PSSM (position-specific scoring matrix)
#  2. psipred, porter, and sam (predict-2nd) to predict secondary structure
#  3. the rosetta fragment picker to generate fragment libraries
#
#
#  Authors: Dylan Chivian, James Thompson, David E Kim
#
###############################################################################
# USER CONFIGURATION  -- ONLY CHANGE THIS SECTION.
#                        YOU MUST EDIT THESE BEFORE USING
###############################################################################

# BLAST install path (Requires non-blast+ NCBI version)
# recommended: blast-2.2.17
# ftp://ftp.ncbi.nih.gov/blast/executables/release/2.2.17/
local $BLAST_DIR = "/my_src_dir/blast-2.2.17_x64/blast-2.2.17";
local $BLAST_NUM_CPUS = 4;    # number of processors to use

# NR blast database filename
# Create blast db files using 'formatdb  -o T -i nr'
# ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
local $NR = "/my_database_dir/blast/nr";

# DEFAULT PDB STRUCTURE, PSSM, and SPARTA (chemical shift) DATABASE (VALL)
# '-vall_files' option overrides this
local $VALL = "/my_database_dir/vall.apr24.2008.extended";

# ROSETTA PICKER
local $PICKER =
  "/my_src_dir/mini/bin/picker.linuxgccrelease";
local $ROSETTA_DATABASE =
  "/my_database_dir/minirosetta_database";

# This is an optional script that you can provide to launch picker jobs on
# free nodes to run them in parallel. A separate picker job is run for each
# fragment size. The command passed to this script is the picker command.
# If left as an empty string, picker jobs will run serially.
local $SLAVE_LAUNCHER = "";

# THE FOLLOWING IS NOT REQUIRED IF YOU RUN SECONDARY STRUCTURE PREDICTIONS
# MANUALLY AND USE THE FOLLOWING OPTIONS:
# -psipredfile <psipred file>
# -porterfile  <porter file>
# -samfile     <sam file>
# The porter and sam files should be in psipred_ss2 format using
# 'ss_pred_converter.py'.

# PSIPRED install path (secondary structure prediction software)
# http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/
local $PSIPRED_DIR =
  "/my_src_dir/psipred";    # psipred installation directory

# pfilt filtered NR blast database filename used for PSIPRED (see PSIPRED readme)
local $PFILTNR = "/my_database_dir/blast/nr_pfilt";

# SAM install path
# http://compbio.soe.ucsc.edu/sam2src/
local $SAM_DIR = "/my_src_dir/sam3.5.x86_64-linux_current";

# SAM predict-2nd install path
# Secondary structure prediction software using SAM sequence alignment
# http://users.soe.ucsc.edu/~karplus/predict-2nd/
local $SAM_PREDICT_2ND_DIR =
  "/my_src_dir/predict-2nd";    # sam predict-2nd directory

# PORTER (secondary structure prediction software)
# http://distill.ucd.ie/porter/
local $PORTER = "/my_src_dir/porter_current/runPorter.pl";

### END OF USER CONFIGURATION #################################################
###############################################################################

use Cwd;
use bytes;

$DEBUG   = 0;

# tail of output fragment file names
local $TAIL = "_v1_3";

$| = 1;    # disable stdout buffering

local $BLASTPGP = "$BLAST_DIR/bin/blastpgp -a $BLAST_NUM_CPUS ";
local $MAKEMAT =
  "$BLAST_DIR/bin/makemat";    # makemat utility (part of NCBI tools)

local $PSIPRED  = "$PSIPRED_DIR/bin/psipred";   # psipred
local $PSIPASS2 = "$PSIPRED_DIR/bin/psipass2";  # psipass2 (part of psipred pkg)
local $PSIPRED_DATA = "$PSIPRED_DIR/data";  # dir containing psipred data files.

local $SAM           = "$SAM_DIR/bin/target2k";   #target99";     # sam target99
local $SAM_uniqueseq = "$SAM_DIR/bin/uniqueseq";  # sam uniqueseq

local $VALL_BLAST_DB =
  "$VALL.blast";    # blast database of VALL sequences for homolog detection

# This is not necessary because VALL PDB sequences should exist in $VALL.blast.
# PDB pdb_seqres.txt database filename for homolog search using '-nohoms' option.
# Create blast db files using 'formatdb  -o T -i pdb_seqres.txt'.
# ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt
local $PDB_SEQRES = "";

# default values
local @fragsizes = ( 3, 9 );
#  ( 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 );
local $n_frags      = 200;
local $n_candidates = 1000;

local @vall_files = ($VALL);

# argv
local %opts = &getCommandLineOptions();
local $file = $opts{f};
local $run_dir       = cwd();    # get the full path (needed for sam stuff)
local $no_homs       = 0;
local $no_frags      = 0;
local $psipred_iter  = 1;
local $psipred_hbias = 1;
local $psipred_sbias = 1;
local $runpsipred    = 1;
local $runporter     = 1;
local $runsam        = 1;
local $haspsipred    = 0;
local $hasporter     = 0;
local $hassam        = 0;
local $id            = "temp";
local $chain         = "_";
local $xx            = "aa";
local @checkpoint_matrix;
local $sequence;
local $num_ss_preds  = 0;
local $cleanup       = 1;
local @cleanup_files = ();

my @allss;

chop($run_dir) if ( substr( $run_dir, -1, 1 ) eq '/' );

# wait a random amount of time before proceeding, to avoid disk-choke problems
srand();
sleep( int( rand(10) ) + 1 );

if ( defined( $opts{verbose} ) ) {
    if ( $opts{verbose} == 1 ) {
        $DEBUG = 1;
        print "VERBOSE.\n";
    }
}

if ( defined( $opts{rundir} ) ) {
    $run_dir = $opts{rundir};
    &checkExist( 'd', $run_dir );
    chop($run_dir) if ( substr( $run_dir, -1, 1 ) eq '/' );

    ( !$DEBUG ) || print "RUN DIRECTORY: $run_dir\n";
}

if ( defined( $opts{frag_sizes} ) ) {
    @fragsizes = split( /,/, $opts{frag_sizes} );
    ( !$DEBUG ) || print "FRAGMENT SIZES: " . join( " ", @fragsizes ) . "\n";
}

if ( defined( $opts{n_frags} ) ) {
    $n_frags = $opts{n_frags};
    ( !$DEBUG ) || print "n_frags: $n_frags\n";
}

if ( defined( $opts{n_candidates} ) ) {
    $n_candidates = $opts{n_candidates};
    ( !$DEBUG ) || print "n_candidates: $n_candidates\n";
}

if ( defined( $opts{vall_files} ) ) {
    @vall_files = ();
    foreach my $vall ( split( /,/, $opts{vall_files} ) ) {
        if ( -s $vall ) {
            push( @vall_files, $vall );
        }
        else {
            warn "WARNING!! vall file $vall does not exist\n";
        }
    }
}

if ( defined( $opts{homs} ) ) {
    if ( $opts{homs} == 0 ) {
        $no_homs = 1;
        ( !$DEBUG ) || print "EXCLUDE HOMOLOGS.\n";
    }
}

if ( defined( $opts{frags} ) ) {
    if ( $opts{frags} == 0 ) {
        $no_frags = 1;
        ( !$DEBUG ) || print "don't make fragments.\n";
    }
}

if ( defined( $opts{psipred} ) ) {
    if ( $opts{psipred} == 0 ) {
        $runpsipred = 0;
        ( !$DEBUG ) || print "don't run psipred.\n";
    }
}

if ( defined( $opts{psipred_iter} ) ) {
    $psipred_iter = $opts{psipred_iter};
    ( !$DEBUG ) || print "psipred_iter = $psipred_iter\n";
}

if ( defined( $opts{psipred_hbias} ) ) {
    $psipred_hbias = $opts{psipred_hbias};
    ( !$DEBUG ) || print "psipred_hbias = $psipred_hbias\n";
}

if ( defined( $opts{psipred_sbias} ) ) {
    $psipred_sbias = $opts{psipred_sbias};
    ( !$DEBUG ) || print "psipred_sbias = $psipred_sbias\n";
}

if ( defined( $opts{porter} ) ) {
    if ( $opts{porter} == 0 ) {
        $runporter = 0;
        ( !$DEBUG ) || print "don't run porter.\n";
    }
}

if ( defined( $opts{sam} ) ) {
    if ( $opts{sam} == 0 ) {
        $runsam = 0;
        ( !$DEBUG ) || print "don't run sam.\n";
    }
}

if ( defined( $opts{xx} ) ) {
    $xx = substr( $opts{xx}, 0, 2 );
    ( !$DEBUG ) || print "picker xx code: $xx\n";
}

if ( defined( $opts{cleanup} ) ) {
    if ( $opts{cleanup} == 0 ) {
        $cleanup = 0;
    }
}

( !$DEBUG ) || print "FILENAME: $file\n";

if ( !defined( $opts{id} ) ) {
    ( !$DEBUG ) || print "no id specified. parse filename instead.\n";

    ($id) = $file =~ /(\w+\.\w+)$/;
    ( !$DEBUG ) || print "INTERMEDIATE: $id\n";
    ($id) = $id =~ /^(\w+)/;

    if ( length($id) != 5 ) {
        die(
"DANGER WILL ROBINSON! DANGER! Your fasta filename is more than/less than five letters!\n"
              . "You should either 1) rename your file to something of the form *****.fasta -or-\n"
              . "You should explicitly specify the id and chain with the -id option.\n"
        );
    }

    $chain = substr( $id, 4, 1 );
    $id    = substr( $id, 0, 4 );
    ( !$DEBUG ) || print "ID: $id CHAIN: $chain\n";

}
else {
    chomp $opts{id};

    ( !$DEBUG ) || print "id specified by user: $opts{id}\n";

    if ( length( $opts{id} ) != 5 ) {
        die("The id you specify must be 5 characters long.\n");
    }

    if ( $opts{id} =~ /\W+/ ) {
        die("Only alphanumeric characters and _ area allowed in the id.\n");
    }

    $id    = substr( $opts{id}, 0, 4 );
    $chain = substr( $opts{id}, 4, 1 );

    ( !$DEBUG ) || print "ID: $id CHAIN: $chain\n";
}

#########
# determine what ss predictions to run
#########

# default SS prediction file names if predictions are run by this script
local $psipred_file = "$run_dir/$id$chain.psipred_ss2";
local $sam_file     = "$run_dir/$id$chain.rdb_ss2";
local $porter_file  = "$run_dir/$id$chain.porter_ss2";

if ( defined( $opts{psipredfile} ) ) {
    $psipred_file = $opts{psipredfile};
    if ( !&fileExist($psipred_file) ) {
        print
"Specified psipred file ($psipred_file) not found!  Running psipred instead.\n";
        $psipred_file = "$run_dir/$id$chain.psipred_ss2";
        $runpsipred   = 1;
    }
}

if ( &fileExist($psipred_file) ) {
    $runpsipred = 0;
    print "Assuming $psipred_file is a psipred_ss2 file.\n";
    my $ssstr = &read_ss2($psipred_file);
    push( @allss, { type => 'psipred', ss => $sstr, file => $psipred_file } );
    $haspsipred = 1;
}

###############
if ( defined( $opts{porterfile} ) ) {
    $porter_file = $opts{porterfile};
    if ( !&fileExist($porter_file) ) {
        print
"Specified porter file ($porter_file) not found!  Running porter instead.\n";
        $porter_file = "$run_dir/$id$chain.porter_ss2";
        $runporter   = 1;
    }
}

if ( &fileExist($porter_file) ) {
    $runporter = 0;
    print "Assuming $porter_file is a porter file.\n";
    my $ssstr = &read_ss2($porter_file);
    push( @allss, { type => 'porter', ss => $sstr, file => $porter_file } );
    $hasporter = 1;
}
###############

if ( defined( $opts{samfile} ) ) {
    $sam_file = $opts{samfile};
    if ( !&fileExist($sam_file) ) {
        print
          "Specified sam file ($sam_file) not found!  Running sam instead.\n";
        $sam_file = "$run_dir/$id$chain.rdb_ss2";
        $runsam   = 1;
    }
}

if ( &fileExist($sam_file) ) {
    $runsam = 0;
    print "Assuming $sam_file is a sam file.\n";
    my $ssstr = &read_ss2($sam_file);
    push( @allss, { type => 'sam', ss => $sstr, file => $sam_file } );
    $hassam = 1;
}

###############################################################################
# main
###############################################################################

chdir($run_dir);

# get the sequence from the fasta file
open( SEQFILE, $file );

my $has_comment = 0;
my $has_eof     = 0;
my $eof;
while (<SEQFILE>) {
    $eof = $_;
    s/\s//g;
    if (/^>/) { $has_comment = 1; next; }
    chomp;
    $sequence .= $_;
}
$has_eof = 1 if ( $eof =~ /\n$/ );
( $has_comment && $has_eof )
  or die "fasta file must have a comment and end with a new line!\n";

( !$DEBUG ) || print "Sequence: $sequence\n";

# run blast
unless ( &fileExist("$id$chain.nr.chk") ) {
    if (
        !&try_try_again(
"$BLASTPGP -t 1 -i $file -F F -j2 -o $id$chain.nr.blast -d $NR -v10000 -b10000 -K1000 -h0.0009 -e0.0009 -C $id$chain.nr.chk -Q $id$chain.nr.pssm",
            2,
            ["$id$chain.nr.chk"],
            [
                "$id$chain.nr.chk",  "$id$chain.nr.blast",
                "$id$chain.nr.pssm", "error.log"
            ]
        )
      )
    {
        die("checkpoint psi-blast failed!\n");
    }
}

unless ( &fileExist("$id$chain.checkpoint") ) {

    # parse & fortran-ify the checkpoint matrix.
    @checkpoint_matrix = &parse_checkpoint_file("$id$chain.nr.chk");
    @checkpoint_matrix =
      &finish_checkpoint_matrix( $sequence, @checkpoint_matrix );
    &write_checkpoint_file( "$id$chain.checkpoint", $sequence,
        @checkpoint_matrix );
}

push( @cleanup_files,
    ( "$id$chain.nr.blast", "$id$chain.nr.pssm", "error.log" ) );

#############################################
# Secondary Structure Prediction methods
#############################################

#
# preliminaries -- run psi-blast for psipred
#
if ($runpsipred) {
    unless ( &fileExist("sstmp.chk")
        && &fileExist("sstmp.ascii")
        && &fileExist("ss_blast") )
    {
        ## PB - add -b -v args since ss_blast will be used by PROF
        if (
            !&try_try_again(
"$BLASTPGP -t 1 -b10000 -v10000 -j3 -h0.001 -d $PFILTNR -i $file -C sstmp.chk -Q sstmp.ascii -o ss_blast",
                2,
                [ "sstmp.chk", "sstmp.ascii" ],
                [ "sstmp.chk", "sstmp.ascii", "ss_blast" ]
            )
          )
        {
            die("psipred psi-blast failed!\n");
        }
        push( @cleanup_files, ( "ss_blast", "sstmp.chk", "sstmp.ascii" ) );
    }
}

#
# psipred
#
if ($runpsipred) {
    ( !$DEBUG ) || print "running psipred.\n";

    unless ( &fileExist("psitmp.sn") ) {
        &run( "echo $file > psitmp.sn", ("psitmp.sn") );
    }
    unless ( &fileExist("psitmp.pn") ) {
        &run( "echo sstmp.chk > psitmp.pn", ("psitmp.pn") );
    }
    unless ( &fileExist("sstmp.mtx") ) {
        if (
            !&try_try_again(
                "$MAKEMAT -P psitmp",
                2, ["sstmp.mtx"], ["sstmp.mtx"]
            )
          )
        {
            die("psipred: makemat failed.");
        }
    }
    unless ( &fileExist("psipred_ss") ) {
        if (
            !&try_try_again(
"$PSIPRED sstmp.mtx $PSIPRED_DATA/weights.dat $PSIPRED_DATA/weights.dat2 $PSIPRED_DATA/weights.dat3 $PSIPRED_DATA/weights.dat4 > psipred_ss",
                2,
                ["psipred_ss"],
                ["psipred_ss"]
            )
          )
        {
            die("psipred failed.");
        }
    }

    unless ( &fileExist("psipred_horiz") ) {
        if (
            !&try_try_again(
"$PSIPASS2 $PSIPRED_DATA/weights_p2.dat $psipred_iter $psipred_hbias $psipred_sbias psipred_ss2 psipred_ss > psipred_horiz",
                2,
                [ "psipred_ss2", "psipred_horiz" ],
                [ "psipred_ss2", "psipred_horiz" ]
            )
          )
        {
            die("psipred/psipass2 failed.");
        }
        rename( "psipred_horiz", "$id$chain.psipred" )
          || die("couldn't move psipred_horiz to $id$chain.psipred: $!\n");
        if ( !scalar(`grep 'PSIPRED VFORMAT' psipred_ss2`) ) {
            open( FILE, "psipred_ss2" );
            my @ss2 = <FILE>;
            close(FILE);
            open( FILE, ">$id$chain.psipred_ss2" );
            print FILE "# PSIPRED VFORMAT (PSIPRED V2.6 by David Jones)\n\n";
            print FILE @ss2;
            close(FILE);
        }
        else {
            rename( "psipred_ss2", "$id$chain.psipred_ss2" )
              || die(
                "couldn't move psipred_ss2 to $id$chain.psipred_ss2: $!\n");
        }
    }

    if ( &fileExist("$run_dir/$id$chain.psipred_ss2") ) {
        ( !$DEBUG ) || print "psipred file ok.\n";
        my $ssstr = &read_ss2("$run_dir/$id$chain.psipred_ss2");
        push(
            @allss,
            {
                type => 'psipred',
                ss   => $sstr,
                file => "$run_dir/$id$chain.psipred_ss2"
            }
        );
        $haspsipred = 1;
    }
    else {
        print "psipred run failed!\n";
    }
    push( @cleanup_files, ( <psitmp*>, "psipred_ss", "sstmp.mtx" ) );
}

if ($runporter) {
    ( !$DEBUG ) || print "running porter.\n";
    if (
        !&try_try_again(
            "$PORTER $file", 2,
            ["$id$chain.porter_ss2"], ["$id$chain.porter_ss2"]
        )
      )
    {
        die("porter failed!\n");
    }
    else {
        my $ssstr = &read_ss2("$run_dir/$id$chain.porter_ss2");
        push(
            @allss,
            {
                type => 'porter',
                ss   => $sstr,
                file => "$run_dir/$id$chain.porter_ss2"
            }
        );
        $hasporter = 1;
    }

# cleanup?
# 1utg_.chk  1utg_.dataset  1utg_.dataset.porter  1utg_.fasta  1utg_.flatblast  1utg_.msa0  1utg_.out0  1utg_.pssm  error.log  ModelsGlobalPaths.porter.txt
}

#
# sam -- target99 alignment and predict2nd with 6 state neural net - condensed output to 3 state
#
if ($runsam) {
    ( !$DEBUG ) || print "running sam.\n";

    my $target99_out      = $id . $chain . ".target99";
    my $target99_a2m_file = $target99_out . ".a2m";

    ## run target99 for a2m alignment
    if (
        !&try_try_again(
            "$SAM -seed $file -out " . $target99_out,
            2, [$target99_a2m_file], [$target99_out_a2m_file]
        )
      )
    {
        die "sam target99 failed!\n";
    }

    ## run uniqueseq
    my $uniqueseq_a2m_id   = $id . $chain . ".uniqueseq";
    my $uniqueseq_a2m_file = $uniqueseq_a2m_id . ".a2m";
    if (
        !&try_try_again(
            "$SAM_uniqueseq $uniqueseq_a2m_id -percent_id 0.9 -alignfile "
              . $target99_a2m_file,
            2,
            [$uniqueseq_a2m_file],
            [$uniqueseq_a2m_file]
        )
      )
    {
        die "sam uniqueseq failed!\n";
    }
    ## run predict-2nd
    chop($SAM_PREDICT_2ND_DIR)
      if ( substr( $SAM_PREDICT_2ND_DIR, -1, 1 ) eq '/' );

    # create samscript
    my $sam_6state        = "$run_dir/$id$chain.sam_6state";
    my $sam_ebghtl        = "$run_dir/$id$chain.sam_ebghtl";
    my $sam_log           = "$run_dir/$id$chain.sam_log";
    my $samscript_txt_buf = qq{ReadAlphabet $SAM_PREDICT_2ND_DIR/std.alphabet
	       ReadAlphabet $SAM_PREDICT_2ND_DIR/DSSP.alphabet
		   ReadNeuralNet $SAM_PREDICT_2ND_DIR/overrep-3617-IDaa13-7-10-11-10-11-7-7-ebghtl-seeded3-stride-trained.net
		       ReadA2M $run_dir/$uniqueseq_a2m_file
			   PrintRDB $sam_6state
			       PrintPredictionFasta $sam_ebghtl
			       };
    my $samscript_txt = "$run_dir/$id$chain.samscript.txt";
    open( SAMSCRIPT, '>' . $samscript_txt );
    print SAMSCRIPT $samscript_txt_buf;
    close(SAMSCRIPT);

    # get into $SAM_PREDICT_2ND_DIR so sam can read file recode3.20comp
    chdir $SAM_PREDICT_2ND_DIR;
    if (
        !&try_try_again(
"$SAM_PREDICT_2ND_DIR/predict-2nd -noalph < $samscript_txt >& $sam_log",
            2,
            [$sam_6state],
            [$sam_6state]
        )
      )
    {
        die "sam predict-2nd failed!\n";
    }

    chdir $run_dir;

    # condense to 3 state prediction
    # formerly a separate perl script: condense_6state.pl
    # integrated here 8/1/03   C. Rohl

    my $rows;
    open IN, "<$sam_6state" || die "Cannot open file $sam_6state: $!\n";
    @{$rows} = <IN>;
    close(IN);
    $rows = &condense_rdb6_rdb($rows);
    open OUT, ">$id$chain.rdb" || die "Cannot open file $id$chain.rdb: $!\n";
    print OUT @{$rows};
    close(OUT);

    # convert rdb to rdb_ss2
    &convert_sam_ss("$id$chain.rdb");
    ( -s "$id$chain.rdb_ss2" ) or die("SAM failed!\n");
    my $ssstr = &read_ss2("$run_dir/$id$chain.rdb_ss2");
    push( @allss,
        { type => 'sam', ss => $sstr, file => "$run_dir/$id$chain.rdb_ss2" } );
    $hassam = 1;

    ( !$DEBUG ) || print "sam file ok.\n";

    push(
        @cleanup_files,
        (
            $samscript_txt, $sam_log,           $sam_6state,
            $sam_ebghtl,    $target99_a2m_file, $uniqueseq_a2m_file,
            "$target99_out.cst"
        )
    );

}
if ($DEBUG) {
    print "\n";
    print "              $sequence\n";
    foreach my $ss (@allss) {
        printf "SS %10s $ss->{ss} $ss->{file}\n", $ss->{type};
    }
    print "\n";
}

#############################################
# Vall and homologue searches
#############################################

if ($no_homs) {
    foreach my $vall (@vall_files) {
        my $vallname = $vall;
        $vallname =~ s/^.*\/([^\/]+)\s*$/$1/gs;
        unless ( &fileExist("$vallname.blast.outn") ) {
            if (
                !&try_try_again(
"$BLASTPGP -t 1 -i $file -j 1 -R $id$chain.nr.chk -o $vallname.blast.outn -e 0.05 -d $vall.blast",
                    2,
                    ["$vallname.blast.outn"],
                    ["$vallname.blast.outn"]
                )
              )
            {
                die("homolog blast failed!\n");
            }
        }
        &exclude_outn( $id, $chain, "$vallname.blast.outn" );
    }

#    if (!-s $PDB_SEQRES && $DEBUG) {
#        print "WARNING!! pdb_seqres.txt does not exist so skipping homolog check from PDB database.\n";
#    }
    unless ( !-s $PDB_SEQRES || &fileExist("$id$chain.pdb_blast") ) {
        if (
            !&try_try_again(
"$BLASTPGP -t 1 -i $file -j 1 -R $id$chain.nr.chk -o $id$chain.pdb_blast -e 0.05 -d $PDB_SEQRES",
                2,
                ["$id$chain.pdb_blast"],
                ["$id$chain.pdb_blast"]
            )
          )
        {
            die("homolog pdb blast failed!\n");
        }
    }

    &exclude_pdbblast( $id, $chain, "$id$chain.pdb_blast" );
    &exclude_blast( $id, $chain, "$id$chain.blast" );
}

#############################################
# picker
#############################################

# don't run picker?
goto CLEANUP if ($no_frags);

###############################################################################
###############################################################################
#
# PICKER CONFIGURATION
###############################################################################

local $PICKER_QUOTA_DEF = "";
my $ss_pred_str = "";
if ( $haspsipred && $hassam && $hasporter ) {
    $PICKER_QUOTA_DEF = <<QUOTA_DEF;
#pool_id        pool_name       fraction
1               psipred         0.50
2               porter          0.25
3               sam             0.25
QUOTA_DEF
    $ss_pred_str = "$psipred_file psipred $porter_file porter $sam_file sam";
}
elsif ( $haspsipred && $hassam ) {
    $PICKER_QUOTA_DEF = <<QUOTA_DEF;
#pool_id        pool_name       fraction
1               psipred         0.50
2               sam             0.50
QUOTA_DEF
    $ss_pred_str = "$psipred_file psipred $sam_file sam";
}
elsif ( $haspsipred && $hasporter ) {
    $PICKER_QUOTA_DEF = <<QUOTA_DEF;
#pool_id        pool_name       fraction
1               psipred         0.50
2               porter          0.50
QUOTA_DEF
    $ss_pred_str = "$psipred_file psipred $porter_file porter";
}
elsif ( $hasporter && $hassam ) {
    $PICKER_QUOTA_DEF = <<QUOTA_DEF;
#pool_id        pool_name       fraction
1               porter          0.50
2               sam             0.50
QUOTA_DEF
    $ss_pred_str = "$porter_file porter $sam_file sam";
}
elsif ($haspsipred) {
    $PICKER_QUOTA_DEF = <<QUOTA_DEF;
#pool_id        pool_name       fraction
1               psipred         1.00
QUOTA_DEF
    $ss_pred_str = "$psipred_file psipred";
}
elsif ($hasporter) {
    $PICKER_QUOTA_DEF = <<QUOTA_DEF;
#pool_id        pool_name       fraction
1               porter          1.00
QUOTA_DEF
    $ss_pred_str = "$porter_file porter";
}
elsif ($hassam) {
    $PICKER_QUOTA_DEF = <<QUOTA_DEF;
#pool_id        pool_name       fraction
1               sam             1.00
QUOTA_DEF
    $ss_pred_str = "$sam_file sam";
}
else {
    die "ERROR! you must use at least one secondary prediction\n";
}

my $ss_score_str = "";
$ss_score_str .= "SecondarySimilarity     200     0.5     -       psipred\n"
  if ($haspsipred);
$ss_score_str .= "SecondarySimilarity     200     0.5     -       porter\n"
  if ($hasporter);
$ss_score_str .= "SecondarySimilarity     200     0.5     -       sam\n"
  if ($hassam);
local $PICKER_SCORE_CFG = <<PICKER_SCORE_CFG;
# score name          priority  wght   min_allowed  extras
$ss_score_str\ProfileScoreL1          300     1.0     -
RamaScore               100     2.0     -
PICKER_SCORE_CFG

open( F, ">$id$chain\_Quota.def" )
  or die "ERROR! cannot open $id$chain\_Quota.def: $!\n";
print F $PICKER_QUOTA_DEF;
close(F);
open( F, ">$id$chain\_Scores.cfg" )
  or die "ERROR! cannot open $id$chain\_Scores.cfg: $!\n";
print F $PICKER_SCORE_CFG;
close(F);

my $valls_str = join( " ", @vall_files );
print "VALL files: $valls_str\n" if $DEBUG;

foreach my $size (@fragsizes) {

    open( PATH_DEFS, ">$xx$id$chain\_picker_cmd_size$size.txt" );
    my $cmdtxt = <<CMDTXT;
-in:file:fasta		$file
-in::path::database     $ROSETTA_DATABASE
-in::file::vall         $valls_str
-frags::n_candidates	$n_candidates
-frags::n_frags		$n_frags
-frags::frag_sizes	$size
-frags::picking::quota_config_file  ./$id$chain\_Quota.def
#-frags::describe_fragments test.fsc.score
-out::file::frag_prefix   $xx$id$chain
-frags::scoring::config   ./$id$chain\_Scores.cfg
-out:level 200
-in::file::checkpoint     $id$chain.checkpoint
-frags::ss_pred           $ss_pred_str
CMDTXT

    if ($no_homs) {
        $cmdtxt .= "-frags:denied_pdb         $id$chain.homolog\n";
    }

    print PATH_DEFS $cmdtxt;
    close(PATH_DEFS);
}

my $attempts        = 5;
my $pickerexecregex = $PICKER;
$pickerexecregex =~ s/(\W)/\\$1/g;

## This is an ugly way to run the picker in parallel, sorry for the hackiness - dk
my $slave_rundir = `pwd`;
chomp $slave_rundir;
$slave_rundir =~ s/\/$//;
foreach my $attempt ( 1 .. $attempts ) {

    if ( !$SLAVE_LAUNCHER ) {    # run in series
        foreach my $size (@fragsizes) {
            my $sizestr = sprintf( "%2.2d", $size );
            my ($resultfile1) = <$xx$id$chain$sizestr*$TAIL>;
            my ($resultfile2) = <$xx$id$chain*\.$size\mers>;
            unless ( -s $resultfile1 || -s $resultfile2 ) {
                my $shell =
"$PICKER \@$xx$id$chain\_picker_cmd_size$size.txt >& $xx$id$chain\_picker_cmd_size$size.log";
                print "shell: $shell\n" if ($DEBUG);
                ( system($shell) == 0 ) || die "ERROR! $shell command failed\n";
            }
        }
    }
    else {    # run in parallel
        print "Attempt $attempt\n" if ($DEBUG);

        my %submitted;
        my @ps = `ps ux`;
        foreach my $size (@fragsizes) {
            my $sizestr = sprintf( "%2.2d", $size );
            my ($resultfile1) = <$xx$id$chain$sizestr*$TAIL>;
            my ($resultfile2) = <$xx$id$chain*\.$size\mers>;
            my $pickerregex   = "$xx$id$chain\_picker_cmd_size$size.txt";
            $pickerregex =~ s/(\W)/\\$1/g;
            my $pwdregex = $slave_rundir;
            $pwdregex =~ s/(\W)/\\$1/g;
            unless (
                   -s $resultfile1
                || -s $resultfile2
                || scalar(
                    grep { /$pickerregex/ && /$pickerexecregex/ && /$pwdregex/ }
                      @ps
                )
              )
            {
                my $shell =
"$SLAVE_LAUNCHER $PICKER \@$xx$id$chain\_picker_cmd_size$size.txt >& $xx$id$chain\_picker_cmd_size$size.log &";
                print "shell: $shell\n" if ($DEBUG);
                system($shell);
                sleep(2);    # give the disks a little break
                $submitted{$size} = 1;
            }
            else {
                $submitted{$size} = 1
                  if (
                    scalar(
                        grep {
                            /$pickerregex/ && /$pickerexecregex/ && /$pwdregex/
                          } @ps
                    )
                  );
            }
        }

        # wait for jobs to complete
        my $time = time();
        my %done;
        my %noprocess_nofile;
        if ( scalar( keys %submitted ) ) {
            while (1) {
                @ps = `ps ux`;
                foreach my $size ( keys %submitted ) {
                    next if ( $done{$size} );
                    my $pickerregex = "$xx$id$chain\_picker_cmd_size$size.txt";
                    $pickerregex =~ s/(\W)/\\$1/g;
                    my $pwdregex = $slave_rundir;
                    $pwdregex =~ s/(\W)/\\$1/g;
                    my ($resultfile) = <$xx$id$chain*\.$size\mers>;
                    if (
                        $resultfile && $resultfile !~ /^\s*$/gs && !scalar(
                            grep {
                                     /$pickerregex/
                                  && /$pickerexecregex/
                                  && /$pwdregex/
                              } @ps
                        )
                      )
                    {
                        sleep(30);
                        if ( -s $resultfile ) {
                            $done{$size} = 1;
                            print "Size $size done!\n" if ($DEBUG);
                        }
                        else {
                            $noprocess_nofile{$size}++;
                            if ( $noprocess_nofile{$size} > 20 ) {
                                warn
"WARNING!! giving up on Size $size, process and result file do not exist, will probably try to run again.\n";
                                $done{$size} = 1;
                            }
                        }
                    }
                    elsif (
                        !scalar(
                            grep {
                                     /$pickerregex/
                                  && /$pickerexecregex/
                                  && /$pwdregex/
                              } @ps
                        )
                      )
                    {
                        $noprocess_nofile{$size}++;
                        if ( $noprocess_nofile{$size} > 20 ) {

                            # give up on this one
                            warn
"WARNING!! giving up on Size $size, process and result file do not exist, will probably try to run again.\n";
                            $done{$size} = 1;
                        }
                    }
                }
                my $diff = time() - $time;
                last
                  if ( scalar( keys %done ) == scalar( keys %submitted )
                    || $diff > 5 * 60 * 60 );    # max wait of 5 hours
                sleep(120);
            }
        }

    }

    # convert files to old style names
    my @frags = `ls -1 $xx$id$chain*mers`;
    foreach my $frag (@frags) {
        chomp $frag;
        if ( $frag =~ /^$xx$id$chain.*\.(\d+)mers$/ ) {
            eval {
                ## check frag file format!
                &check_fragformat( $file, $frag, $DEBUG );
            };
            if ($@) {

                # format is bad so remove it and try again
                unlink($frag);
            }
            else {
                my $size = "00$1";
                $size =~ s/^.*(\d{2})$/$1/gc;
                my $newfrag = "$xx$id$chain$size\_05.200$TAIL";
                system("mv $frag $newfrag");
                ( -s $newfrag ) or die "picker failed!\n";
            }
        }
    }

    my %checkdone;
    foreach my $size (@fragsizes) {
        my $sizestr = sprintf( "%2.2d", $size );
        my ($resultfile1) = <$xx$id$chain$sizestr*$TAIL>;
				if ( $resultfile1 && -s $resultfile1 ) {
            $checkdone{$resultfile1} = 1;
            print "$resultfile1\n" if $DEBUG;
        }
    }

    if ( scalar( keys %checkdone ) == scalar(@fragsizes) ) {
        last;
    }
    else {
        warn
"WARNING! all fragment files do not exist after attempt $attempt, attempting to run picker again\n";
        if ( $attempt == $attempts ) {
            die "ERROR! picker failed after $attempt attempts\n";
        }
    }
}

push(
    @cleanup_files,
    (
        <zscore_*_$TAIL_$id$chain>, <status.*$TAIL_$xx$id>,
        <names.*$TAIL$xx$id$chain>
    )
);

#############################################
# cleanup
#############################################

CLEANUP:

if ($cleanup) {
    my $file;

    foreach $file (@cleanup_files) {
        unlink($file);
    }
}

print "Done!\n" if ($DEBUG);

# done
exit 0;

###############################################################################
# util
###############################################################################

sub check_fragformat {
    my ( $fasta, $frag, $verbose ) = @_;
    chomp $fasta;
    chomp $frag;
    my $seq;
    print "Checking frag file format: $frag\n" if ($verbose);
    open( F, $fasta ) or die "ERROR! cannot open $fasta: $!\n";
    while (<F>) {
        next if (/^>/);
        $seq .= $_;
    }
    close(F);
    $seq =~ s/\s+//gs;
    my ( @n, $pos, $nbrs, $nbrscnt, $prevnbrs, $fragsize );
    my $poscnt = 0;
    open( F, $frag ) or die "ERROR! cannot open $frag: $!\n";
    while (<F>) {
        if (/^\s*position:\s+(\d+)\s+neighbors:\s+(\d+)\s*$/) {
            $prevnbrs = $nbrs;
            $pos      = $1;
            $nbrs     = $2;
            $poscnt++;
            if ( $prevnbrs != $nbrscnt ) {
                die "ERROR! neighbor truncated\n";
            }
            $nbrscnt = 0;
        }
        elsif (/^\s+\w\w\w\w\s+[\w\-]\s+\d+/) {
            push( @n, $_ );
        }
        elsif (/^\s*$/) {
            if ( scalar @n ) {
                if ( !$nbrscnt ) {
                    $fragsize = scalar @n;
                }
                else {
                    if ( scalar @n != $fragsize ) {
                        die "ERROR! neighbor count off\n";
                    }
                }
                $nbrscnt++;
            }
            @n = ();
        }
    }
    close(F);

    if ( $nbrs != $nbrscnt ) {
        die "ERROR! neighbor truncated\n";
    }

    if ( length($seq) != $pos + $fragsize - 1 ) {
        my $len    = length($seq);
        my $poslen = $pos + $fragsize - 1;
        die "ERROR! fasta seq length $len does not match positions $poslen\n";
    }
    print "Format okay!\n" if ($verbose);
}

# getCommandLineOptions()
#
#  rets: \%opts  pointer to hash of kv pairs of command line options
#
sub getCommandLineOptions {
    use Getopt::Long;

    my $fragsizes = join( ",", @fragsizes );
    my $usage = qq{usage: $0
		[
		   \t-rundir <full path to run directory (default = ./)>
		   \t-verbose  specify for chatty output
		   \t-id  <5 character pdb code/chain id, e.g. 1tum_>
		   \t-xx  <silly little 2-letter code, if you care>
		   \t-nohoms  specify to omit homologs from search
		   \t-nopsipred  don\'t run psipred.
		   \t-psipred_iter  \# of psipred iterations (default = 1)
		   \t-psipred_hbias  psipred helix bias (default = 1)
		   \t-psipred_sbias  psipred strand bias (default = 1)
		   \t-psipredfile  <path to file containing psipred ss prediction>
		   \t-nosam   don\'t run sam.
		   \t-samfile  <path to file containing sam ss prediction>
		   \t-noporter  don\'t run porter.
		   \t-porterfile  <path to file containing porter ss prediction>
		   \t-nofrags  specify to not make fragments and just run SS predictions
		   \t-vall_files <vall1,vall2,...>  (default = $VALL)
		   \t-frag_sizes <size1,size2,...n> (default = $fragsizes)
                   \t-n_frags <number of fragments> (default = $n_frags)
                   \t-n_candidates <number of candidates> (default = $n_candidates)
		   \t-nocleanup  specify to keep all the temporary files
		]\n
		<fasta file>\n};

    # Get args
    my %opts = ();
    &GetOptions(
        \%opts,            "psipred!",
        "psipred_iter=f",  "psipred_hbias=f",
        "psipred_sbias:f", "sam!",
        "homs!",           "porter!",
        "psipredfile=s",   "samfile=s",
        "porterfile=s",    "xx=s",
        "verbose!",        "rundir=s",
        "id=s",            "frags!",
        "cleanup!",        "frag_sizes=s",
        "n_frags=i",       "n_candidates=i",
        "vall_files=s"
    );

    if ( scalar(@ARGV) != 1 ) {
        die "$usage\n";
    }

    $opts{f} = $ARGV[0];

    if ( $opts{f} =~ /[^\w\.\/]/ ) {
        die(
"Only alphanumeric characters, . and _ are allowed in the filename.\n"
        );
    }

    &checkExist( "f", $opts{f} );

    if ( defined $opts{frag_sizes} && $opts{frag_sizes} !~ /^[\d\,]+$/ ) {
        die("Fragment sizes are invalid\n");
    }

    return %opts;
}

# checkExist()
#
sub checkExist {
    my ( $type, $path ) = @_;
    if ( $type eq 'd' ) {
        if ( !-d $path ) {
            print STDERR "$0: dirnotfound: $path\n";
            exit -3;
        }
    }
    elsif ( $type eq 'f' ) {
        if ( !-f $path ) {
            print STDERR "$0: filenotfound: $path\n";
            exit -3;
        }
        elsif ( !-s $path ) {
            print STDERR "$0: emptyfile: $path\n";
            exit -3;
        }
    }
}

sub fileExist {
    my $path = shift;

    return 0 if ( !-f $path );
    return 0 if ( -z $path );

    return 1;
}

sub run {
    my ( $cmd, @files ) = @_;
    my $pid;

    if ($DEBUG) {
        print "cmd is: $cmd\n";
    }

  FORK: {
        if ( $pid = fork ) {

            # parent
            local $SIG{TERM} = sub {
                kill 9, $pid;
                `rm -f @files`;
                exit;
            };

            local $SIG{INT} = sub {
                kill 9, $pid;
                `rm -f @files`;
                exit;
            };

            wait;

            return $?;

        }
        elsif ( defined $pid ) {

            #child
            exec($cmd);
        }
        elsif ( $! =~ /No more process/ ) {

            # recoverable error
            sleep 5;
            redo FORK;

        }
        else {

            # unrecoverable error
            die("couldn't fork: $!\n");
        }
    }
}

sub try_try_again {
    my ( $cmd, $max_tries, $success_files, $cleanup_files ) = @_;
    my $try_count     = 0;
    my $missing_files = 1;
    my $f;

    # if at first you don't succeed in running $cmd, try, try again.
    while ( ( $try_count < $max_tries ) && ( $missing_files > 0 ) ) {
        sleep 10;
        &run( $cmd, @$cleanup_files );
        ++$try_count;

        $missing_files = 0;

        foreach $f (@$success_files) {
            if ( !&fileExist($f) ) {
                ++$missing_files;
            }
        }
    }

    # but if you've tried $max_tries times, give up.
    # there's no use being a damn fool about things.
    if ( $missing_files > 0 ) {
        return 0;
    }

    return 1;
}

sub exclude_outn {
    my ( $pdb, $chain, $outn ) = @_;
    my @hits;
    my $hit;
    my $hit_pdb;
    my $hit_chain;
    my %uniq_hits;

    open( OUTN,  "$outn" );
    open( EXCL,  ">>$pdb$chain.homolog_vall" );    # append
    open( EXCL2, ">>$pdb$chain.homolog" );

    @hits = grep( /^>/, <OUTN> );

    foreach $hit (@hits) {
        $uniq_hits{$hit} = 1;
    }

    foreach $hit ( sort keys %uniq_hits ) {
        ($hit_pdb) = $hit =~ /^>\s*(\w+)/;
        $hit_chain = substr( $hit_pdb, 4, 1 );
        $hit_pdb   = substr( $hit_pdb, 0, 4 );

        $hit_pdb =~ tr/[A-Z]/[a-z]/;
        $hit_chain = '_' if ( $hit_chain eq ' ' );
        $hit_chain = '_' if ( $hit_chain eq '0' );

        print EXCL "$pdb$chain $hit_pdb$hit_chain\n";
        print EXCL2 "$pdb$chain $hit_pdb$hit_chain\n";
    }
    close(OUTN);
    close(EXCL);
    close(EXCL2);
}

sub exclude_pdbblast {
    my ( $pdb, $chain, $blast ) = @_;
    my @hits;
    my $hit;
    my $hit_pdb;
    my $hit_chain;
    my %uniq_hits;
    return if ( !-s $blast );
    open( BLAST, $blast );
    open( EXCL,  ">$pdb$chain.homolog_pdb" );
    open( EXCL2, ">>$pdb$chain.homolog" );

    @hits = grep( /^>/, <BLAST> );

    foreach $hit (@hits) {
        $uniq_hits{$hit} = 1;
    }

    foreach $hit ( sort keys %uniq_hits ) {
        ($hit_pdb) = $hit =~ /^>\s*(\w+)/;
        $hit_chain = substr( $hit_pdb, 4, 1 );
        $hit_pdb   = substr( $hit_pdb, 0, 4 );

        $hit_pdb =~ tr/[A-Z]/[a-z]/;
        $hit_chain = '_' if ( $hit_chain eq ' ' );
        $hit_chain = '_' if ( $hit_chain eq '0' );

        print EXCL "$pdb$chain $hit_pdb$hit_chain\n";
        print EXCL2 "$pdb$chain $hit_pdb$hit_chain\n";
    }

    close(EXCL);
    close(EXCL2);
    close(BLAST);
}

sub exclude_blast {
    my ( $pdb, $chain, $blast ) = @_;
    my @hits;
    my $hit;
    my $hit_pdb;
    my $hit_chain;
    my %uniq_hits;

    open( BLAST, "$blast" );
    open( EXCL,  ">$pdb$chain.homolog_nr" );
    open( EXCL2, ">>$pdb$chain.homolog" );

    @hits = grep( /pdb\|/, <BLAST> );

    foreach $hit (@hits) {
        $uniq_hits{$hit} = 1;
    }

    foreach $hit ( sort keys %uniq_hits ) {
        while ( $hit =~ s/pdb\|(\w{4})\|(.?)// ) {
            $hit_pdb   = $1;
            $hit_chain = $2;

            $hit_pdb =~ tr/[A-Z]/[a-z]/;
            $hit_chain = '_' if ( $hit_chain =~ /^\s*$/ );

            print EXCL "$pdb$chain $hit_pdb$hit_chain\n";
            print EXCL2 "$pdb$chain $hit_pdb$hit_chain\n";
        }
    }
    close(BLAST);
    close(EXCL);
    close(EXCL2);
}

## parse_checkpoint_file -- parses a PSI-BLAST binary checkpoint file.
#
# args:  filename of checkpoint file.
# rets:  N x 20 array containing checkpoint weight values, where N
#        is the size of the protein that BLAST thought it saw...

sub parse_checkpoint_file {
    my $filename = shift;
    my $buf;
    my $seqlen;
    my $seqstr;
    my $i;
    my $j;
    my @aa_order = split( //, 'ACDEFGHIKLMNPQRSTVWY' );
    my @altschul_mapping =
      ( 0, 4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15, 16, 19, 17, 18 );
    my @w;
    my @output;

    open( INPUT, $filename ) || die("Couldn't open $filename for reading.\n");

    read( INPUT, $buf, 4 ) || die("Couldn't read $filename!\n");
    $seqlen = unpack( "i", $buf );

    ( !$DEBUG ) || print "Sequence length: $seqlen\n";

    read( INPUT, $buf, $seqlen ) || die("Premature end: $filename.\n");
    $seqstr = unpack( "a$seqlen", $buf );

    for ( $i = 0 ; $i < $seqlen ; ++$i ) {
        read( INPUT, $buf, 160 ) || die("Premature end: $filename, line: $i\n");
        @w = unpack( "d20", $buf );

        for ( $j = 0 ; $j < 20 ; ++$j ) {
            $output[$i][$j] = $w[ $altschul_mapping[$j] ];
        }
    }

    return @output;
}

## finish_checkpoint_matrix -- "finishes" the parsed PSI-BLAST checkpoint matrix,
##                             by adding pseudo-counts to any empty columns.
#
# args:  1) sequence string  2) array returned by parse_checkpoint_file
# rets:  "finished" array.  suitable for printing, framing, etc.

sub finish_checkpoint_matrix {
    my ( $s, @matrix ) = @_;
    my @sequence = split( //, $s );
    my $i;
    my $j;
    my $line;
    my $sum;
    my @words;
    my @b62;
    my @blos_aa =
      ( 0, 14, 11, 2, 1, 13, 3, 5, 6, 7, 9, 8, 10, 4, 12, 15, 16, 18, 19, 17 );
    my %aaNum = (
        A => 0,
        C => 1,
        D => 2,
        E => 3,
        F => 4,
        G => 5,
        H => 6,
        I => 7,
        K => 8,
        L => 9,
        M => 10,
        N => 11,
        P => 12,
        Q => 13,
        R => 14,
        S => 15,
        T => 16,
        V => 17,
        W => 18,
        Y => 19,
        X => 0     ### cheep fix for now ##
    );

    ( length($s) == scalar(@matrix) )
      || die("Length mismatch between sequence and checkpoint file!\n");

    my $blosum62 = <<BLOSUM;
#  BLOSUM Clustered Target Frequencies=qij
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
   A      R      N      D      C      Q      E      G      H      I      L      K      M      F      P      S      T      W      Y      V
0.0215
0.0023 0.0178
0.0019 0.0020 0.0141
0.0022 0.0016 0.0037 0.0213
0.0016 0.0004 0.0004 0.0004 0.0119
0.0019 0.0025 0.0015 0.0016 0.0003 0.0073
0.0030 0.0027 0.0022 0.0049 0.0004 0.0035 0.0161
0.0058 0.0017 0.0029 0.0025 0.0008 0.0014 0.0019 0.0378
0.0011 0.0012 0.0014 0.0010 0.0002 0.0010 0.0014 0.0010 0.0093
0.0032 0.0012 0.0010 0.0012 0.0011 0.0009 0.0012 0.0014 0.0006 0.0184
0.0044 0.0024 0.0014 0.0015 0.0016 0.0016 0.0020 0.0021 0.0010 0.0114 0.0371
0.0033 0.0062 0.0024 0.0024 0.0005 0.0031 0.0041 0.0025 0.0012 0.0016 0.0025 0.0161
0.0013 0.0008 0.0005 0.0005 0.0004 0.0007 0.0007 0.0007 0.0004 0.0025 0.0049 0.0009 0.0040
0.0016 0.0009 0.0008 0.0008 0.0005 0.0005 0.0009 0.0012 0.0008 0.0030 0.0054 0.0009 0.0012 0.0183
0.0022 0.0010 0.0009 0.0012 0.0004 0.0008 0.0014 0.0014 0.0005 0.0010 0.0014 0.0016 0.0004 0.0005 0.0191
0.0063 0.0023 0.0031 0.0028 0.0010 0.0019 0.0030 0.0038 0.0011 0.0017 0.0024 0.0031 0.0009 0.0012 0.0017 0.0126
0.0037 0.0018 0.0022 0.0019 0.0009 0.0014 0.0020 0.0022 0.0007 0.0027 0.0033 0.0023 0.0010 0.0012 0.0014 0.0047 0.0125
0.0004 0.0003 0.0002 0.0002 0.0001 0.0002 0.0003 0.0004 0.0002 0.0004 0.0007 0.0003 0.0002 0.0008 0.0001 0.0003 0.0003 0.0065
0.0013 0.0009 0.0007 0.0006 0.0003 0.0007 0.0009 0.0008 0.0015 0.0014 0.0022 0.0010 0.0006 0.0042 0.0005 0.0010 0.0009 0.0009 0.0102
0.0051 0.0016 0.0012 0.0013 0.0014 0.0012 0.0017 0.0018 0.0006 0.0120 0.0095 0.0019 0.0023 0.0026 0.0012 0.0024 0.0036 0.0004 0.0015 0.0196
BLOSUM

    $i = 0;
    my @blosum62 = split( /\n/, $blosum62 );

    # read/build the blosum matrix
    foreach (@blosum62) {
        next if ( $_ !~ /^\d/ );
        chomp;
        @words = split(/\s/);
        for ( $j = 0 ; $j <= $#words ; ++$j ) {
            $b62[ $blos_aa[$i] ][ $blos_aa[$j] ] = $words[$j];
            $b62[ $blos_aa[$j] ][ $blos_aa[$i] ] = $words[$j];
        }

        ++$i;
    }

    # normalize the blosum matrix so that each row sums to 1.0
    for ( $i = 0 ; $i < 20 ; ++$i ) {
        $sum = 0.0;

        for ( $j = 0 ; $j < 20 ; ++$j ) {
            $sum += $b62[$i][$j];
        }

        for ( $j = 0 ; $j < 20 ; ++$j ) {
            $b62[$i][$j] = ( $b62[$i][$j] / $sum );
        }
    }

    # substitute appropriate blosum row for 0 rows
    for ( $i = 0 ; $i <= $#matrix ; ++$i ) {
        $sum = 0;

        for ( $j = 0 ; $j < 20 ; ++$j ) {
            $sum += $matrix[$i][$j];
        }

        if ( $sum == 0 ) {
            for ( $j = 0 ; $j < 20 ; ++$j ) {
                $matrix[$i][$j] = $b62[ $aaNum{ $sequence[$i] } ][$j];
            }
        }
    }

    return @matrix;
}

sub write_checkpoint_file {
    my ( $filename, $sequence, @matrix ) = @_;
    my $row;
    my $col;
    my @seq = split( //, $sequence );

    open( OUTPUT, ">$filename" );

    die("Length mismatch between sequence and checkpoint matrix!\n")
      unless ( length($sequence) == scalar(@matrix) );

    print OUTPUT scalar(@matrix), "\n";

    for ( $row = 0 ; $row <= $#matrix ; ++$row ) {
        print OUTPUT "$seq[$row] ";
        for ( $col = 0 ; $col < 20 ; ++$col ) {
            printf OUTPUT "%6.4f ", $matrix[$row][$col];
        }
        print OUTPUT "\n";
    }

    print OUTPUT "END";

    close(OUTPUT);
}

## verify_existing_checkpoint
#
# Does what it says...it verifies a checkpoint file to make sure it's kosher.
#
# arg: filename of checkpoint file.
# rets: 1 if good, 0 if bad.

sub verify_existing_checkpoint {
    my $filename = shift;
    my @file;

    open( CPT, $filename ) || return 0;

    ( !$DEBUG ) || print "opened checkpoint file.\n";

    @file = <CPT>;

    ( $file[0] =~ /^(\d+)$/ ) || return 0;

    ( !$DEBUG ) || print "starting line ok.\n";

    ( ( scalar(@file) - 2 ) == $1 ) || return 0;

    ( !$DEBUG ) || print "length ok.\n";

    ( $file[$#file] =~ /^END/ ) || return 0;

    ( !$DEBUG ) || print "end ok.\n";

    close(CPT);
    return 1;
}

## convert 6 state rdb file to 3 state
sub condense_rdb6_rdb {
    my ($rows) = @_;
    my ( $outrows, @parts, $line );
    for ( @{$rows} ) {
        if (/^#/) {
            push( @{$outrows}, $_ );
        }
        elsif (/^Pos/) {
            push( @{$outrows}, "Pos\tAA\tE\tH\tL\n" );
        }
        elsif (/^10N/) {
            push( @{$outrows}, "10N\t1S\t5N\t5N\t5N\n" );
        }
        elsif (/^[0-9]/) {
            @parts = split;
            die "Wrong number of columns in results row: $_\n" if $#parts != 7;
            $line = sprintf "%8s\t%1s\t%8s\t%8s\t%8s\n", $parts[0] - 1,
              $parts[1], $parts[2] + $parts[3], $parts[4] + $parts[5],
              $parts[6] + $parts[7];
            push( @{$outrows}, $line );
        }
        else {
            die "Unrecognized row: $_\n";
        }
    }
    return $outrows;
}

sub convert_sam_ss {
    my $ss_pred = shift;
    my ( $e, $h, $c, @ss );
    my $plusOne = 0;
    open( SAM, $ss_pred ) or die "ERROR! cannot open $ss_pred: $!\n";
    open( SAMCONV, ">$ss_pred\_ss2" )
      or die "ERROR! cannot open $ss_pred\_ss2: $!\n";
    print SAMCONV "# PSIPRED VFORMAT (PSIPRED V2.6 by David Jones)\n";
    while (<SAM>) {
        $_ =~ s/^\s+|\s+$//;
        next if ( $_ =~ /^#/ || $_ =~ /^\s*$/ );
        my @cols = split( /\s+/, $_ );
        next if ( scalar @cols != 5 || $cols[0] !~ /^\d+$/ );
        $e       = $cols[2];
        $h       = $cols[3];
        $c       = $cols[4];
        $plusOne = 1 if ( $cols[0] eq '0' );

        if ( $e >= $h && $e >= $c ) {
            printf SAMCONV "%4d $cols[1] E   %5.3f  %5.3f  %5.3f\n",
              $cols[0] + $plusOne, $h, $e, $c;
        }
        elsif ( $h > $e && $h >= $c ) {
            printf SAMCONV "%4d $cols[1] H   %5.3f  %5.3f  %5.3f\n",
              $cols[0] + $plusOne, $h, $e, $c;
        }
        elsif ( $c > $e && $c > $h ) {
            printf SAMCONV "%4d $cols[1] C   %5.3f  %5.3f  %5.3f\n",
              $cols[0] + $plusOne, $h, $e, $c;
        }
    }
    close(SAMCONV);
    close(SAM);
}

sub read_ss2 {
    my $ss2 = shift;
    $sstr = "";
    open( SS2, $ss2 ) or die "ERROR! cannot open $ss2: $!\n";
    while (<SS2>) {
        if (/^\s*\d+\s+\w\s+(\w)\s+\d/) {
            $sstr .= $1;
        }
    }
    close(SS2);
    if (!$sstr) {
        die "ERROR! cannot extract secondary structure prediction from $ss2: psipred_ss2 format required.\n";
    }
    return $sstr;
}




