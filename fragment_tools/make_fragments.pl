#!/usr/bin/perl -w
use strict;
# Resolve symlinks to find source deploymnet
use FindBin qw( $RealBin );
my $Bin = $RealBin;
use constant VERSION => 3.10;
$| = 1; # disable stdout buffering

# Install and build Rosetta and that's it. See fragments.README.
#
# This script initially installs all other dependencies. The installation takes
# a long time due to downloading and formatting the NCBI non-redundant sequence
# database. If you have any dependencies already installed, you can optionally
# modify the configuration below. The installation requires ~73+ Gigs of free
# disk space.

###############################################################################
# OPTIONAL USER CONFIGURATION  -- ONLY CHANGE THIS SECTION IF NECESSARY
###############################################################################

# REQUIRED
#  Rosetta				https://www.rosettacommons.org
#
#  The script "install_dependencies.pl" installs external dependencies:
#    ./blast               ftp://ftp.ncbi.nih.gov/blast/executables/release/2.2.17/
#    ./databases/nr        ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
#    ./psipred             http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/
#    ./databases/nr_pfilt  from ./databases/nr using ./psipred/bin/pfilt
#    ./sparks-x            http://sparks-lab.org/yueyang/server/SPARKS-X/
#    ./csblast             https://github.com/cangermueller/csblast

# ROSETTA
my $FRAGMENT_PICKER = "$Bin/../../main/source/bin/fragment_picker.boost_thread.linuxgccrelease";
my $FRAGMENT_PICKER_NUM_CPUS = 8;    # number of processors to use
my $ROSETTA_DATABASE = "$Bin/../../main/database"; # rosetta database
my $VALL = "$Bin/vall.jul19.2011"; # template database

# BLAST path (Requires non-blast+ NCBI version)
my $BLAST_DIR = "$Bin/blast";
my $BLAST_NUM_CPUS = 8;    # number of processors to use (blastpgp -a option)

# NR database path
my $NR = "$Bin/databases/nr";

# spine-x/sparks (for phi, psi, and solvent accessibility predictions)
my $SPARKS = "$Bin/sparks-x/bin/buildinp_query.sh";

# PSIPRED (for secondary structure prediction)
my $PSIPRED_DIR = "$Bin/psipred";
my $PSIPRED_USE_weights_dat4 = 0;    # set to 0 if using psipred version 3.2+

# CSBLAST/CSBUILD (for de-novo sequence profile generation)
my $CSBLAST_DIR = "$Bin/csblast";

# pfilt filtered NR database used for PSIPRED (see PSIPRED readme)
# $NR will be used if empty
my $PFILTNR = "$Bin/databases/nr_pfilt";

my $INTERNET_HOST = "localhost";

### EXTRA OPTIONAL FEATURES ###################################################
###############################################################################

# This is an optional script that YOU must provide to launch parallel jobs on
# your cluster. The script should take any command as the argument(s).
# If the script does not exist, jobs will run serially.
# Requires: http://search.cpan.org/CPAN/authors/id/D/DL/DLUX/Parallel-ForkManager-0.7.5.tar.gz
my $SLAVE_LAUNCHER = "";
my $SLAVE_LAUNCHER_MAX_JOBS = 40;    # depends on your available machines/cpus
## for SLAVE_LAUNCHER parallel jobs
my $SLAVE_MAX_WAIT     = 96 * 60 * 60;
my $SLAVE_MAX_ATTEMPTS = 2;

# pdb2vall.py script for adding specific PDBs to the vall (-add_pdbs_to_vall)
#  --no_structure_profile option is added to reduce the run time
# This feature is not supported yet in the Rosetta release
my $PDB2VALL = "$Bin/pdb2vall/pdb2vall.py --no_structure_profile";
my $PDB2VALL_IGNORE_ERRORS = 1;    # ignore pdb2vall jobs that fail

# The following can be ignored unless you want to use the secondary structure prediction
# quota system with Psipred, SAM, and Porter. The porter should be in psipred_ss2 format
# using 'ss_pred_converter.py'.

# SAM install path
# http://compbio.soe.ucsc.edu/sam2src/
my $SAM_DIR = "";

# SAM predict-2nd install path
# Secondary structure prediction software using SAM sequence alignment
# http://users.soe.ucsc.edu/~karplus/predict-2nd/
my $SAM_PREDICT_2ND_DIR = "";

# PORTER (secondary structure prediction software)
# http://distill.ucd.ie/porter/
my $PORTER = "";

my $INSTALL_DEPENDENCIES = "standard"; # "overwrite";
my $INSTALL_DEPENDENCIES_DATABASE = "nr"; # "uniref50" # "uniref90"

### YOU CAN IGNORE THE REST ###################################################
###############################################################################


use File::Path;
use File::Copy qw/ copy /;
use File::Basename;
use Time::Local;


use Cwd qw/ cwd abs_path /;
use bytes;

my %options;

# initialize options
my %opts = &getCommandLineOptions();
$options{fastafile} = abs_path( $opts{f} );
$options{rundir}      = cwd();     # get the full path (needed for sam stuff)
$options{homs}        = 1;
$options{frags}  = 1;
$options{csbuild_profile} = 0;
$options{psipredfile} = "";
$options{samfile}     = "";
$options{porterfile}  = "";
$options{psipred}     = 1;         # use psipred by default
$options{porter}      = 0;         # skip porter by default
$options{sam}         = 0;         # skip sam by default
$options{id}          = "temp";
$options{chain}       = "_";
$options{runid}       = "temp_";
$options{cleanup}     = 1;
$options{torsion_bin} = 0;
$options{exclude_homologs_by_pdb_date} = 0;
$options{old_name_format}              = 0;
$options{add_pdbs_to_vall}             = "";

my @cleanup_files  = ();
my @fragsizes      = ( 3, 9 );  #4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20 );
my @add_vall_files = ();
my @use_vall_files = ();
my @pdbs_to_vall   = ();
$options{n_frags}      = 200;
$options{n_candidates} = 1000;

foreach my $key ( keys %opts ) {
    $options{$key} = $opts{$key};
}

$options{DEBUG} = 1 if ( $options{verbose} );

### check for dependencies and install if necessary
$BLAST_DIR = "$Bin/blast" if (!-d "$BLAST_DIR/bin" || !-d "$BLAST_DIR/data");
$SPARKS = "$Bin/sparks-x/bin/buildinp_query.sh" if (!-s $SPARKS);
$PSIPRED_DIR = "$Bin/psipred" if (!-d "$PSIPRED_DIR/bin" || !-d "$PSIPRED_DIR/data");

my @POSSIBLE_NR_LOCATIONS = (
	$NR,
	"$Bin/databases/nr",
	"$Bin/../../../databases/nr/nr",
	"/scratch/robetta/local_db/nr/nr",
	"/work/robetta/databases/local_db/nr/nr"
);

my $skip_nr = "";
foreach my $n (@POSSIBLE_NR_LOCATIONS) {
	if (-s "$n.pal") {
		$NR = $n;
		$skip_nr = "skip_nr";
		last;
	}
}
foreach my $n (@POSSIBLE_NR_LOCATIONS) {
	if (-s $n."_pfilt.pal") {
		$PFILTNR = $n."_pfilt";
		last;
	}
}

my $must_install = 0;
if (!-d "$BLAST_DIR/bin" || !-d "$BLAST_DIR/data") {
	print "Dependency 'blast' does not exist!\n";
	$must_install = 1;
} elsif (!-s $SPARKS) {
	print "Dependency 'sparks-x' does not exist!\n";
	$must_install = 1;
} elsif (!-d "$PSIPRED_DIR/bin" || !-d "$PSIPRED_DIR/data") {
	print "Dependency 'psipred' does not exist!\n";
	$must_install = 1;
}

if (!$options{csbuild_profile}) {
  if (!-s "$NR.pal") {
    print "Dependency 'nr' does not exist!\n";
    $NR = "$Bin/databases/nr";
    $must_install = 1;
  } elsif (!-s "$PFILTNR.pal") {
    print "Dependency 'nr_pfilt' does not exist!\n";
    $PFILTNR = "$Bin/databases/nr_pfilt";
    $skip_nr = "";
    $must_install = 1;
  }
} else {
  $skip_nr = "skip_nr";
}

if ($must_install && $INSTALL_DEPENDENCIES) {
	print "\n";
	print "Running install_dependencies.pl\n";
	print "Note: the NCBI non-redundant (nr) sequence database is very large.\n";
	print "Please be patient.....\n\n";
	system("$Bin/install_dependencies.pl $INSTALL_DEPENDENCIES $INSTALL_DEPENDENCIES_DATABASE $skip_nr");
	print "\n";
}

## THESE FILE LOCATIONS DEPEND ON THE SOFTWARE PACKAGE
my $PSIBLAST = "$BLAST_DIR/bin/blastpgp -a $BLAST_NUM_CPUS ";
my $MAKEMAT = "$BLAST_DIR/bin/makemat";   # makemat utility (part of NCBI tools)
my $PSIPRED = "$PSIPRED_DIR/bin/psipred"; # psipred
my $PSIPASS2 = "$PSIPRED_DIR/bin/psipass2";    # psipass2 (part of psipred pkg)
my $PSIPRED_DATA = "$PSIPRED_DIR/data";    # dir containing psipred data files.
my $SAM           = "$SAM_DIR/bin/target99";     # sam target99
my $SAM_uniqueseq = "$SAM_DIR/bin/uniqueseq";    # sam uniqueseq
my $CSBUILD = "$CSBLAST_DIR/bin/csbuild";

# check for picker
if (!-s $FRAGMENT_PICKER) {
  warn "WARNING! $FRAGMENT_PICKER does not exist. Trying default.\n";
  $FRAGMENT_PICKER =~ s/fragment_picker[^\/]*(\.[^\.]+)$/fragment_picker$1/;
	if (!-s $FRAGMENT_PICKER) {
		$FRAGMENT_PICKER =~ s/fragment_picker(\.[^\.]+)$/fragment_picker.default$1/;
	}
  die "ERROR! $FRAGMENT_PICKER does not exist.\n" if (!-s $FRAGMENT_PICKER);
}

# check nr database
if (!$options{csbuild_profile}) {
  if (!-s $PFILTNR) {
    $PFILTNR = $NR;
    print_debug("nr_pfilt database missing so using nr: $NR");
  }

  (-s $NR) or die "ERROR! $NR does not exist.\n";
}


# for homolog detection
my $VALL_BLAST_DB = "$VALL.blast";
$VALL_BLAST_DB =~ s/\.gz\.blast$/\.blast/;
my $PDB_SEQRES      = "$Bin/pdb_seqres.txt";
if ( !$options{homs} ) {
	if (!-s "$VALL_BLAST_DB.phr") {
		system("gunzip $VALL_BLAST_DB.gz") if (-s "$VALL_BLAST_DB.gz");
		system("$BLAST_DIR/bin/formatdb -i $VALL_BLAST_DB") if (-s $VALL_BLAST_DB);
	}
	if (!-s $PDB_SEQRES) {
		if ($INTERNET_HOST) {
			my $pwd = cwd();
			system("ssh $INTERNET_HOST 'cd $pwd; wget ftp://ftp.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt'");
		} else {
			system("wget ftp://ftp.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt");
		}
		if (!-s "pdb_seqres.txt") { die "ERROR! wget ftp://ftp.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt failed\n"; }
		system("mv pdb_seqres.txt $PDB_SEQRES") if (!-s $PDB_SEQRES);
	}
	if (!-s "$PDB_SEQRES.phr") {
		warn "WARNING! $PDB_SEQRES does not seem to be formated. Trying to format.\n";
		system("$BLAST_DIR/bin/formatdb -i $PDB_SEQRES");
	}
}

if ( defined( $opts{frag_sizes} ) ) {
    @fragsizes = split( /,/, $opts{frag_sizes} );
    print_debug( "fragment sizes: " . join( " ", @fragsizes ) );
}

if ( exists $opts{add_vall_files} ) {
    foreach my $vall ( split( /,/, $opts{add_vall_files} ) ) {
        push( @add_vall_files, $vall ) if ( -s $vall );
    }
}

if ( exists $opts{use_vall_files} ) {
    foreach my $vall ( split( /,/, $opts{use_vall_files} ) ) {
        push( @use_vall_files, $vall ) if ( -s $vall );
    }
}

if ( defined $opts{add_pdbs_to_vall} ) {
    if ( !$PDB2VALL ) {
        print_debug(
            "ignoring -add_pdbs_to_vall option: PDB2VALL not configured");
    }
    else {
        @pdbs_to_vall = split( /,/, $opts{add_pdbs_to_vall} );
	# make sure the PDB code is lower case
	for (my $i=0;$i<=$#pdbs_to_vall;$i++) {
		substr($pdbs_to_vall[$i],0,4) = lc substr($pdbs_to_vall[$i],0,4);
	}
        print_debug( "pdbs to vall: " . join( " ", @pdbs_to_vall ) );
    }
}
mkpath( $options{rundir} );
$options{rundir} = abs_path( $options{rundir} );
chop( $options{rundir} ) if ( substr( $options{rundir}, -1, 1 ) eq '/' );
&checkExist( 'd', $options{rundir} );

if ( !defined( $opts{id} ) ) {
    print_debug("no id specified. parsing filename instead.");

    ( $options{id} ) = $options{fastafile} =~ /(\w+\.\w+)$/;
    ( $options{id} ) = $options{id} =~ /^(\w+)/;

    if ( length( $options{id} ) != 5 ) {
	print_debug("cannot parse id from filename so using 't001_'");
        $options{id} = 't001_';
    }

    $options{chain} = substr( $options{id}, 4, 1 );
    $options{id}    = substr( $options{id}, 0, 4 );
    print_debug("ID: $options{id} CHAIN: $options{chain}");
}
else {
    chomp $opts{id};

    print_debug("id specified by user: $opts{id}");

    if ( length( $opts{id} ) != 5 ) {
        die("The id you specify must be 5 characters long.\n");
    }

    if ( $opts{id} =~ /\W+/ ) {
        die("Only alphanumeric characters and _ area allowed in the id.\n");
    }

    $options{id}    = substr( $opts{id}, 0, 4 );
    $options{chain} = substr( $opts{id}, 4, 1 );

    if ( -s "$options{id}$options{chain}.fasta" ) {
        my @diff = `diff $options{id}$options{chain}.fasta $options{fastafile}`;
        if ( scalar @diff ) {
            die
"ERROR! $options{id}$options{chain}.fasta already exists but does not match $options{fastafile}: $options{id}$options{chain}.fasta is autogenerated and should not exist\n";
        }
    }
    print_debug("using $options{fastafile} as query fasta");
    print_debug("ID: $options{id} CHAIN: $options{chain}");
}

$options{runid} = "$options{id}$options{chain}";

if ( abs_path( $options{fastafile} ) ne abs_path("$options{runid}.fasta") ) {
    copy( $options{fastafile}, "$options{runid}.fasta" );
    $options{fastafile} = "$options{runid}.fasta";
}

# determine what ss predictions to run
foreach my $ss_pred (qw/ porter sam psipred /) {
    if ( $options{$ss_pred} ) {
        my $fn_key = join '', ( $ss_pred, 'file' );
        my $fn = $options{$fn_key};
        $options{$ss_pred} = file_overrides_option( $fn, $ss_pred );
    }
}

my $abs_path_fasta = abs_path( $options{fastafile} );
$options{fastafile} = basename( $options{fastafile} );
if ( $abs_path_fasta ne abs_path("$options{rundir}/$options{fastafile}") ) {
    copy( $options{fastafile}, "$options{rundir}/" )
      or die "Error copying $options{fastafile} into $options{rundir}!\n";
}

print "picking fragments with options:\n", options_to_str( \%options ), "\n";
print_debug("FILENAME: $options{fastafile}");

# main

chdir( $options{rundir} );

if ( -f "$options{runid}.make_fragments.success" ) {
    print_debug("done! $options{runid}.make_fragments.success exists");
    exit(0);
}

# get the sequence from the fasta file
my $sequence = read_fasta( $options{fastafile} );
print_debug("Sequence: $sequence");

# Run csbuild to generate sequence profile (.check) and pssm (.pssm)
# This bypasses psiblast calls for profile generation in sparks and psipred
if ($options{csbuild_profile}) {
  print_debug("Generating structure profile & pssm via csbuild.");
  system("$CSBUILD -i $options{fastafile} -I fas -D $CSBLAST_DIR/data/K4000.crf -o $options{runid}.check -O chk");
  system("cp $options{runid}.check sstmp.chk");
  system("echo sstmp.chk > sstmp.pn");
  system("echo $options{fastafile} > sstmp.sn");
  system("$MAKEMAT -P sstmp");
  (-s "sstmp.mtx") or die "ERROR! Failed to create .mtx file: makemat failed!\n";

  my $placeholder = "$Bin/pdb2vall/structure_profile_scripts/placeholder_seqs/placeholder_seqs";
  if (!-s "$placeholder.phr") {
    system("$BLAST_DIR/bin/formatdb -o T -i $placeholder");
  }

  my $blast = "$options{runid}.blast";
  open(F, ">$blast") or die "ERROR! cannot open $blast: $!\n";
  my $aacnt = 0;
  my $totalaacnt = 0;
  foreach my $aa (split(//,$sequence)) {
          $aacnt++;
          $totalaacnt++;
          if ($aacnt == 1) {
                  printf F "%-14.14s", $options{runid};
          }
          if ($aacnt >= 101) {
                  print F "$aa\n\n";
                  $aacnt = 0;
                  next;
          }
          print F $aa;
          print F "\n\n" if ($totalaacnt >= length($sequence));
  }
  close(F);

  system("$PSIBLAST -i $options{fastafile} -B $blast -Q $options{runid}.pssm -t 1 -j 1 -h 0.001 -e 0.001 -b 0 -k 0 -d $placeholder");
  (-s "$options{runid}.pssm") or die "ERROR! failed to create single sequence pssm file: blastpgp failed!\n";
  system("cp $options{runid}.pssm $options{fastafile}.pssm");
}

# run sparks
if ($SPARKS) {
    unless ( &nonempty_file_exists( $options{fastafile} . ".phipsi" ) ) {
        print_debug(
            "Running sparks for phi, psi, and solvent accessibility predictions"
        );
        my $sparks_result = system("$SPARKS $options{fastafile}");
        if (!($sparks_result == 0 && -s $options{fastafile} . ".phipsi" )) {
          if (length($SPARKS) > 100) {
            print "sparks path length > 100 characters, which may cause internal errors. Try moving to a shorter path prefix.\n";
          }
          die("sparks failed!\n");
        }
    }
}

# run blast
unless ( &nonempty_file_exists("$options{runid}.check") ) {
    print_debug("Running psiblast for sequence profile");
    print_debug("Using nr: $NR");
    if (
        !&try_try_again(
"$PSIBLAST -t 1 -i $options{fastafile} -F F -j2 -o $options{runid}.blast -d $NR -v10000 -b10000 -K1000 -h0.0009 -e0.0009 -C $options{runid}.check -Q $options{runid}.pssm",
            2,
            ["$options{runid}.check"],
            [
                "$options{runid}.check", "$options{runid}.blast",
                "$options{runid}.pssm",  "error.log"
            ]
        )
      )
    {
        die("checkpoint psi-blast failed!\n");
    }
}

unless ( &nonempty_file_exists("$options{runid}.checkpoint") ) {

    # parse & fortran-ify the checkpoint matrix.
    my @checkpoint_matrix;
    @checkpoint_matrix = &parse_checkpoint_file("$options{runid}.check");
    @checkpoint_matrix =
      &finish_checkpoint_matrix( $sequence, @checkpoint_matrix );
    &write_checkpoint_file( "$options{runid}.checkpoint", $sequence,
        @checkpoint_matrix );
}

push( @cleanup_files, ( "$options{runid}.pssm", "error.log" ) );

# Secondary Structure Prediction methods
if ( $options{psipred}
    || ( $options{psipredfile} && -s $options{psipredfile} ) )
{
    if ( $options{psipredfile} && -s $options{psipredfile} ) {
        $options{psipred} = 1;
        system("cp $options{psipredfile} $options{runid}.psipred_ss2");
    }
    else {

        # run psi-blast for psipred
        unless ( &nonempty_file_exists("sstmp.chk") )
        {
            print_debug("Running psi-blast for psipred, using $PFILTNR.");
            my $psiblast_success = &try_try_again(
              "$PSIBLAST -t 1 -b10000 -v10000 -j3 -h0.001 -d $PFILTNR -i $options{fastafile} -C sstmp.chk -Q sstmp.ascii -o ss_blast",
              2, [ "sstmp.chk", "sstmp.ascii" ], [ "sstmp.chk", "sstmp.ascii", "ss_blast" ]);
            $psiblast_success or die("psipred psi-blast failed!\n");

            push( @cleanup_files, ( "ss_blast", "sstmp.chk", "sstmp.ascii" ) );
        } else {
            print_debug("Using existing sstmp.chk");
        }

        unless ( &nonempty_file_exists("sstmp.mtx") ) {
            &run( "echo $options{fastafile} > psitmp.sn", ("psitmp.sn") );
            &run( "echo sstmp.chk > psitmp.pn", ("psitmp.pn") );
            if (
                !&try_try_again(
                    "$MAKEMAT -P psitmp",
                    2, ["sstmp.mtx"], ["sstmp.mtx"]
                )
              )
            {
                die("psipred: makemat failed.");
            }
        } else {
            print_debug("Using existing sstmp.mtx");
          }

        unless ( &nonempty_file_exists("psipred_ss") ) {
            my $psipredcmd =
              ($PSIPRED_USE_weights_dat4)
              ? "$PSIPRED sstmp.mtx $PSIPRED_DATA/weights.dat $PSIPRED_DATA/weights.dat2 $PSIPRED_DATA/weights.dat3 $PSIPRED_DATA/weights.dat4 > psipred_ss"
              : "$PSIPRED sstmp.mtx $PSIPRED_DATA/weights.dat $PSIPRED_DATA/weights.dat2 $PSIPRED_DATA/weights.dat3 > psipred_ss";
            if (
                !&try_try_again(
                    $psipredcmd, 2, ["psipred_ss"], ["psipred_ss"]
                )
              )
            {
                die("psipred failed.");
            }
        }

        unless ( &nonempty_file_exists("$options{runid}.psipred_ss2") ) {
            if (
                !&try_try_again(
"$PSIPASS2 $PSIPRED_DATA/weights_p2.dat 1 1.0 1.0 psipred_ss2 psipred_ss > psipred_horiz",
                    2,
                    [ "psipred_ss2", "psipred_horiz" ],
                    [ "psipred_ss2", "psipred_horiz" ]
                )
              )
            {
                die("psipred/psipass2 failed.");
            }

            rename( "psipred_horiz", "$options{runid}.psipred" )
              or die(
                "couldn't move psipred_horiz to $options{runid}.psipred: $!\n");

            if ( !scalar(`grep 'PSIPRED VFORMAT' psipred_ss2`) ) {
                open( FILE, "psipred_ss2" );
                my @ss2 = <FILE>;
                close(FILE);
                open( FILE, ">$options{runid}.psipred_ss2" );
                print FILE
                  "# PSIPRED VFORMAT (PSIPRED V2.6 by David Jones)\n\n";
                print FILE @ss2;
                close(FILE);
            }
            else {
                rename( "psipred_ss2", "$options{runid}.psipred_ss2" )
                  or die(
"couldn't move psipred_ss2 to $options{runid}.psipred_ss2: $!\n"
                  );
            }
        }

    }

    if ( &nonempty_file_exists("$options{runid}.psipred_ss2") ) {
        print_debug("psipred file ok.");
    }
    else {
        print "psipred run failed!\n";
    }

    push( @cleanup_files, ( glob("psitmp*"), "psipred_ss", "sstmp.mtx" ) );
}

if ( $options{porter} || ( $options{porterfile} && -s $options{porterfile} ) ) {
    if ( $options{porterfile} && -s $options{porterfile} ) {
        $options{porter} = 1;
				if (&is_ss2_format($options{porterfile})) {
          system("cp $options{porterfile} $options{runid}.porter_ss2");
				} else {
					(system("$Bin/ss_pred_converter.py -p $options{porterfile} > $options{runid}.porter_ss2") == 0 && -s "$options{runid}.porter_ss2") or
						die "ERROR! $options{porterfile} conversion to ss2 format failed. Check porter file format.\n";
				}
    }
    if ( &nonempty_file_exists("$options{runid}.porter_ss2") ) {
        print_debug("porter file ok.");
    }
    else {
	print_debug("Running porter.");
			if (
        !&try_try_again(
            "$PORTER $options{fastafile}", 2,
            ["$options{runid}.porter_ss2"], ["$options{runid}.porter_ss2"]
        )
      ) { die("porter failed!\n"); }
		}
}

# sam -- target99 alignment and predict2nd with 6 state neural net - condensed output to 3 state
if ( $options{sam} || ( $options{samfile} && -s $options{samfile} ) ) {
    if ( $options{samfile} && -s $options{samfile} ) {
        $options{sam} = 1;
        if (&is_ss2_format($options{samfile})) {
          system("cp $options{samfile} $options{runid}.rdb_ss2");
        } elsif (&is_sam_6state($options{samfile})) {
          # condense to 3 state prediction
          my $rows;
          open IN, "<$options{samfile}" or die "Cannot open file $options{samfile}: $!\n";
          @{$rows} = <IN>;
          close(IN);
          $rows = &condense_rdb6_rdb($rows);
          open OUT, ">$options{runid}.rdb"
            or die "Cannot open file $options{runid}.rdb: $!\n";
          print OUT @{$rows};
          close(OUT);

          # convert rdb to rdb_ss2
          &convert_sam_ss("$options{runid}.rdb");
          die("SAM failed!\n") if ( !-s "$options{runid}.rdb_ss2" );
				} else {
          (system("$Bin/ss_pred_converter.py -s $options{samfile} > $options{runid}.rdb_ss2") == 0 && -s "$options{runid}.rdb_ss2") or
            die "ERROR! $options{samfile} conversion to ss2 format failed. Check sam file format.\n";
				}
    }
    if ( &nonempty_file_exists("$options{runid}.rdb_ss2") ) {
        print_debug("sam file ok.");
    }
    else {
	print_debug("Running sam.");
        my $target99_out      = "$options{runid}.target99";
        my $target99_a2m_file = $target99_out . ".a2m";

        if (
            !&try_try_again(
                "$SAM -seed $options{fastafile} -out " . $target99_out, 2,
                [$target99_a2m_file], []
            )
          )
        {
            die "sam target99 failed!\n";
        }

        ## run uniqueseq
        my $uniqueseq_a2m_id   = "$options{runid}.uniqueseq";
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
        my $sam_6state = "$options{rundir}/$options{runid}.sam_6state";
        my $sam_ebghtl = "$options{rundir}/$options{runid}.sam_ebghtl";
        my $sam_log    = "$options{rundir}/$options{runid}.sam_log";
        my $samscript_txt_buf =
          qq{ReadAlphabet $SAM_PREDICT_2ND_DIR/std.alphabet
		ReadAlphabet $SAM_PREDICT_2ND_DIR/DSSP.alphabet
		ReadNeuralNet $SAM_PREDICT_2ND_DIR/overrep-3617-IDaa13-7-10-11-10-11-7-7-ebghtl-seeded3-stride-trained.net
		ReadA2M $options{rundir}/$uniqueseq_a2m_file
		PrintRDB $sam_6state
		PrintPredictionFasta $sam_ebghtl
	};
        my $samscript_txt = "$options{rundir}/$options{runid}.samscript.txt";
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
        chdir $options{rundir};

        # condense to 3 state prediction
        my $rows;
        open IN, "<$sam_6state" or die "Cannot open file $sam_6state: $!\n";
        @{$rows} = <IN>;
        close(IN);
        $rows = &condense_rdb6_rdb($rows);
        open OUT, ">$options{runid}.rdb"
          or die "Cannot open file $options{runid}.rdb: $!\n";
        print OUT @{$rows};
        close(OUT);

        # convert rdb to rdb_ss2
        &convert_sam_ss("$options{runid}.rdb");
        if ( !-s "$options{runid}.rdb_ss2" ) {
            die("SAM failed!\n");
        }

        print_debug("sam file ok.");

        push @cleanup_files,
          (
            $samscript_txt, $sam_log, $sam_6state, $sam_ebghtl,
            $target99_a2m_file, $uniqueseq_a2m_file, "$target99_out.cst"
          );
    }
}

# Vall and homolog searches
my %excluded_homologs;
if (-s "$options{runid}.homolog") {
	open(F, "$options{runid}.homolog") or die "ERROR! cannot open $options{runid}.homolog: $!\n";
	while (<F>) {
		if (/^(\S{5})/) {
			$excluded_homologs{$1} = 1;
		} elsif (/^(\S{4})/) {
			$excluded_homologs{$1} = 1;
		}
	}
	close(F);
}

if ( !$options{homs} ) {

    print_debug("excluding homologs.");
    unless ( &nonempty_file_exists("$options{runid}.outn") ) {
        if (
            !&try_try_again(
"$PSIBLAST -t 1 -i $options{fastafile} -j 1 -R $options{runid}.check -o $options{runid}.outn -e 0.05 -d $VALL_BLAST_DB",
                2,
                ["$options{runid}.outn"],
                ["$options{runid}.outn"]
            )
          )
        {
            die("homolog vall blast failed!\n");
        }
    }

    unless ( &nonempty_file_exists("$options{runid}.outn.pdb") ) {
        if (
            !&try_try_again(
"$PSIBLAST -t 1 -i $options{fastafile} -j 1 -R $options{runid}.check -o $options{runid}.outn.pdb -e 0.05 -d $PDB_SEQRES",
                2,
                ["$options{runid}.outn.pdb"],
                ["$options{runid}.outn.pdb"]
            )
          )
        {
            die("homolog pdb blast failed!\n");
        }
    }

    &exclude_blast( $options{runid} );
    &exclude_pdbblast( $options{runid} );
    &exclude_outn( $options{runid} );
}

my %exclude_homologs_by_pdb_date_pdbs;
my %exclude_homologs_by_pdb_date_checked_pdbs;
my $cutoff_date_str;
if ( $options{exclude_homologs_by_pdb_date} ) {
	my $exclude_homologs_by_pdb_date = $options{exclude_homologs_by_pdb_date};
	print_debug("exclude_homologs_by_pdb_date: $exclude_homologs_by_pdb_date");
	my ($cutoff_m,$cutoff_d,$cutoff_y);
	# mm/dd/yy format
	if ( $exclude_homologs_by_pdb_date =~ /^(\d\d)\/(\d\d)\/(\d\d)$/ ) {
		$cutoff_m = $1;
		$cutoff_d = $2;
		$cutoff_y = $3;
	# yyyymmdd format
	} elsif ( $exclude_homologs_by_pdb_date =~ /^(\d\d\d*)(\d\d)(\d\d)$/ ) {
		$cutoff_y = $1;
		$cutoff_m = $2;
		$cutoff_d = $3;
	} else {
		die "ERROR! $exclude_homologs_by_pdb_date date format is incorrect: must be mm/dd/yy or yyyymmdd\n";
	}
	# try to convert 2 digit year to full calendar year - hacky
	$cutoff_date_str = &convert_twodigit_to_full_calendar_year( $cutoff_y, $cutoff_m, $cutoff_d );

	# update pdb_revdat.txt
	print_debug("Getting latest pdb_revdat.txt....");
	system("$Bin/update_revdat.pl $INTERNET_HOST");

	my %month_str_to_num = (
		'JAN' => '01',
		'FEB' => '02',
		'MAR' => '03',
		'APR' => '04',
		'MAY' => '05',
		'JUN' => '06',
		'JUL' => '07',
		'AUG' => '08',
		'SEP' => '09',
		'OCT' => '10',
		'NOV' => '11',
		'DEC' => '12');

	open( F, "$Bin/pdb_revdat.txt" )
		or die "ERROR! cannot open $Bin/pdb_revdat.txt: $!\n";
	foreach my $l (<F>) {
		if ($l =~ /^(\S{4})\s+REVDAT\s+(\d+)\s+(\d+)\-(\w\w\w)\-(\d+)\s+\S{4}\s+(\d+)/) {
			#2w6x REVDAT   1   22-DEC-09 2W6X    0
			my $rcode = lc $1;
			my $rnum = $2;
			my $rday = $3;
			my $rmonth = $month_str_to_num{$4};
			my $ryear = $5;
			my $rstatus = $6;
			$rday = "00$rday";
			$rday =~ s/^.*(\d\d)\s*$/$1/gs;
			# try to convert 2 digit year to full calendar year
			my $rtmpdate = &convert_twodigit_to_full_calendar_year( $ryear, $rmonth, $rday );
			if ($rstatus eq '0' && $rtmpdate > $cutoff_date_str  ) {
				$exclude_homologs_by_pdb_date_pdbs{$rcode} = $rtmpdate;
			}
			$exclude_homologs_by_pdb_date_checked_pdbs{$rcode} = 1;
		}
	}
	close(F);
}

sub convert_twodigit_to_full_calendar_year {
	my ($twodigit, $tmpmonth, $tmpday) = @_;
	# try to convert 2 digit year to full calendar year - hacky

	# current date GMT  (UTC)  PDB gets updated Wed UTC +0  (-8/7 PST depending on daylight savings)
	my ($sec,$min,$hour,$currentmday,$currentmonth,$currentyear,$wday,$yday,$isdst) = gmtime(time);
	$currentyear += 1900;
	$currentmonth++;
	my $mday = "00$currentmday";
	my $mmonth = "00$currentmonth";
	$mday =~ s/^.*(\d\d)\s*$/$1/gs;
	$mmonth =~ s/^.*(\d\d)\s*$/$1/gs;
	my $current_date_str = $currentyear.$mmonth.$mday;

	my $tmpyear = ($twodigit < 100) ? $twodigit + 2000 : $twodigit;
	my $tmpdate = $tmpyear.$tmpmonth.$tmpday;
	if ($tmpdate > $current_date_str) {
		$tmpyear = ($twodigit < 100) ? $twodigit + 1900 : $twodigit;
		$tmpdate = $tmpyear.$tmpmonth.$tmpday;
	}
	return $tmpdate;
}

my @valls = ($VALL);
if ( scalar @use_vall_files ) {
    @valls = @use_vall_files;
}
push( @valls, @add_vall_files ) if ( scalar @add_vall_files );

my $nativeexists = 0;
if ( -s "$options{runid}.pdb" ) {
    print_debug("native pdb exists: $options{runid}.pdb.");
    $nativeexists = 1;
}

# pdbs to vall
if ( scalar @pdbs_to_vall ) {
    print_debug("pdbs_to_vall");
    my $pdbs2valldir = "$options{rundir}/$options{runid}\_pdbs_to_vall";
    ( -d $pdbs2valldir || mkdir($pdbs2valldir) )
      or die "ERROR! cannot mkdir $pdbs2valldir: $!\n";
    chdir($pdbs2valldir);

    if ( !-s "$options{runid}\_pdbs_to_vall.vall" ) {
        if ( !-e $SLAVE_LAUNCHER ) {    # run in series
            foreach my $pdb (@pdbs_to_vall) {
                my $cmd = "$PDB2VALL -p $pdb";
                produce_output_with_cmd( $cmd, "$pdb.vall",
                    $PDB2VALL_IGNORE_ERRORS );
            }
        }
        else {
            my ( @commands, @results );
            foreach my $pdb (@pdbs_to_vall) {
                push( @commands, "$SLAVE_LAUNCHER $PDB2VALL -p $pdb" );
                push( @results,  "$pdb.vall" );
            }
            &run_in_parallel(
                \@commands,              \@results,
                "pdb2vall_parallel_job", $SLAVE_LAUNCHER_MAX_JOBS,
                $SLAVE_MAX_WAIT,         $SLAVE_MAX_ATTEMPTS,
                $PDB2VALL_IGNORE_ERRORS
            );
        }
    }

    # create one template vall - overwrite if exists
    open( F, ">$options{runid}\_pdbs_to_vall.vall" )
      or die "ERROR! cannot open $options{runid}\_pdbs_to_vall.vall: $!\n";
    foreach my $pdb (@pdbs_to_vall) {
        next if ( $PDB2VALL_IGNORE_ERRORS && !-s "$pdb.vall" );
        print_debug("adding $pdb.vall");
        open( NF, "$pdb.vall" ) or die "ERROR! cannot open $pdb.vall: $!\n";
        foreach my $nf (<NF>) {
            print F $nf;
        }
        close(NF);
    }
    close(F);

    if ( -s "$options{runid}\_pdbs_to_vall.vall" ) {
        push( @valls, "$pdbs2valldir/$options{runid}\_pdbs_to_vall.vall" );
    }
    else {
        warn
"WARNING! pdb2vall output $options{runid}\_pdbs_to_vall.vall does not exist.\n";
    }

		chdir($options{rundir});
}


# check valls for date filtering
if ( $options{exclude_homologs_by_pdb_date} ) {
	chdir($options{rundir});
	my %appended;
	open( EXCL,  ">$options{runid}.homolog_by_pdb_date_$cutoff_date_str" ) or die "ERROR! cannot open!\n";
	open( EXCL2, ">>$options{runid}.homolog" );
	print EXCL "# excluded because release date > $cutoff_date_str or missing in pdb_revdat.txt\n";
	foreach my $v (@valls) {
		if (-s "$v.gz") {
			open(F, "gzip -dc $v.gz |") or die "ERROR! cannot open $v.gz: $!\n";
		} else {
			open(F, $v) or die "ERROR! cannot open $v: $!\n";
		}
		while (<F>) {
			next if (/^#/);
			if (/^(\S{4})\S\s+\S/) {
				my $pdbcode = lc $1;
				next if (exists $appended{$pdbcode});
				if (!exists $exclude_homologs_by_pdb_date_checked_pdbs{$pdbcode}) {
					print_debug("excluding $pdbcode because it is missing in pdb_revdat.txt");
					print EXCL "$pdbcode\n";
					print EXCL2 "$pdbcode\n" if (!exists $excluded_homologs{$pdbcode});
					$appended{$pdbcode} = 1;
				} elsif (exists $exclude_homologs_by_pdb_date_pdbs{$pdbcode}) {
					print_debug("excluding $pdbcode released: $exclude_homologs_by_pdb_date_pdbs{$pdbcode}");
					print EXCL "$pdbcode\n";
					print EXCL2 "$pdbcode\n" if (!exists $excluded_homologs{$pdbcode});
					$appended{$pdbcode} = 1;
				}
			}
		}
		close(F);
	}
	close(EXCL);
	close(EXCL2);
}

my $valls_str = join ' ', @valls;

chdir($options{rundir});


if ( !$options{frags} ) {
	cleanup(@cleanup_files);
	exit 0;
}

# picker
foreach my $size (@fragsizes) {
    my $ss_pred_cnt = 0;
    my $score_def;
    my $sparkscmd = "";
    my $nativecmd = "";

    # score file (small fragments)
    if ( $size <= 4 ) {
        $score_def = <<SCORE;
# score name         priority  wght   min_allowed  extras
ProfileScoreL1           700     1.0     -
ProfileScoreStructL1     100     1.4     -
SCORE

        if ($SPARKS) {
            $score_def .= "SolventAccessibility     500     0.5     -\n";
            $score_def .= "Phi                      300     3.9     -\n";
            $score_def .= "Psi                      200     0.9     -\n";
            $sparkscmd =
              "-spine_x                  $options{fastafile}.phipsi\n";
        }

        $ss_pred_cnt = 0;
        foreach my $ss_pred (qw/ psipred sam porter /) {
            if ( $options{$ss_pred} ) {
                $ss_pred_cnt++;
                $score_def .=
                  "SecondarySimilarity      600     1.0     -    $ss_pred\n";
                $score_def .=
                  "RamaScore                400     6.0     -    $ss_pred\n";
            }
        }

        # score file (large fragments)
    }
    else {
        $score_def = <<SCORE;
# score name         priority  wght   min_allowed  extras
ProfileScoreL1           700     1.0     -
ProfileScoreStructL1     100     4.0     -
SCORE

        if ($SPARKS) {
            $score_def .= "SolventAccessibility     500     1.5     -\n";
            $score_def .= "Phi                      300     1.0     -\n";
            $score_def .= "Psi                      200     0.6     -\n";
            $sparkscmd =
              "-spine_x                  $options{fastafile}.phipsi\n";
        }

        $ss_pred_cnt = 0;
        foreach my $ss_pred (qw/ psipred sam porter /) {
            if ( $options{$ss_pred} ) {
                $ss_pred_cnt++;
                $score_def .=
                  "SecondarySimilarity      600     1.0     -    $ss_pred\n";
                $score_def .=
                  "RamaScore                400     0.8     -    $ss_pred\n";
            }
        }
    }

    if ( $options{torsion_bin} && -f $options{torsion_bin} ) {
        $score_def .= "TorsionBinSimilarity    100     1.0     -\n";
    }

    if ($nativeexists) {
        $score_def .= "FragmentCrmsd            30      0.0     -\n";
        $score_def .= "FragmentDME              30      0.0     -\n";
        $nativecmd = "-in:file:s $options{runid}.pdb\n";
    }

    my $score_file_name = "$options{runid}\_scores$size.cfg";
    open FILE, ">$score_file_name";
    print FILE $score_def;
    close FILE or die $!;

    open( PATH_DEFS, ">$options{runid}\_picker_cmd_size$size.txt" );
    my $cmdtxt = <<CMDTXT;
-in:file:fasta          $options{fastafile}
-in:path:database       $ROSETTA_DATABASE
-in:file:vall           $valls_str
-frags:n_candidates     $options{n_candidates}
-frags:n_frags          $options{n_frags}
-frags:frag_sizes       $size
-out:file:frag_prefix   $options{runid}
-frags:scoring:config   $score_file_name
#-out:level 2000
-in:file:checkpoint     $options{runid}.checkpoint
-frags:write_ca_coordinates
-frags:describe_fragments $options{runid}\_frags.$size.score
$sparkscmd
$nativecmd
CMDTXT

    # add fragment secondary structure predictions
    my @frag_ss_tokens;
    if ( -s "$options{runid}.psipred_ss2" ) {
        push @frag_ss_tokens, ( "$options{runid}.psipred_ss2", "psipred" );
    }
    if ( -s "$options{runid}.rdb_ss2" ) {
        push @frag_ss_tokens, ( "$options{runid}.rdb_ss2", "sam" );
    }
    if ( -s "$options{runid}.porter_ss2" ) {
        push @frag_ss_tokens, ( "$options{runid}.porter_ss2", "porter" );
    }
    $cmdtxt .= join ' ', ( '-frags:ss_pred', @frag_ss_tokens );
    $cmdtxt .= "\n";
    if ( $options{torsion_bin} && -f $options{torsion_bin} ) {
        $cmdtxt .= "-in:file:torsion_bin_probs $options{torsion_bin}\n";
    }

    # auto-generate quota votes if more than 1 ss_pred is used
    if ( $ss_pred_cnt > 1 ) {
        my %votes;
        foreach my $ss_pred (qw/ psipred porter sam /) {
            if ( $options{$ss_pred} ) {
                $votes{$ss_pred} = get_ss_quota_vote($ss_pred);
            }
        }

        if ( scalar( keys %votes ) > 0 ) {
            my $quota_fn = "$options{runid}\_quota.def";
            open QUOTA, ">$quota_fn" or die $!;
            print QUOTA "#pool_id        pool_name       fraction\n";
            my $count = 1;

            use List::Util qw/ sum /;
            my $sum = sum( values %votes );
            foreach my $ss_pred ( keys %votes ) {
                my $pct = sprintf( "%4.2f", $votes{$ss_pred} / $sum );
                print QUOTA "$count          $ss_pred        $pct\n";
                $count++;
            }
            $cmdtxt .= "-frags:picking:quota_config_file $quota_fn\n";
        }
    }

    if ( !$options{homs} || $options{exclude_homologs_by_pdb_date} ) {
        $cmdtxt .= "-frags:denied_pdb  $options{runid}.homolog\n";
    }

    print PATH_DEFS $cmdtxt;
    close(PATH_DEFS);

    if ( !-s $SLAVE_LAUNCHER ) {    # run in series
        my $cmd =
"$FRAGMENT_PICKER \@$options{runid}\_picker_cmd_size$size.txt -j $FRAGMENT_PICKER_NUM_CPUS";
        produce_output_with_cmd( $cmd,
            "$options{runid}.$options{n_frags}.$size" . "mers" );
    }
}

if (-e $SLAVE_LAUNCHER) {           # run in parallel
    my ( @commands, @results );
    foreach my $size (@fragsizes) {
        push( @results, "$options{runid}.$options{n_frags}.$size" . "mers" );
        push( @commands,
"$SLAVE_LAUNCHER $FRAGMENT_PICKER \@$options{runid}\_picker_cmd_size$size.txt -j $FRAGMENT_PICKER_NUM_CPUS"
        );
    }
    &run_in_parallel( \@commands, \@results, "picker_parallel_job",
        $SLAVE_LAUNCHER_MAX_JOBS, $SLAVE_MAX_WAIT, $SLAVE_MAX_ATTEMPTS );
}

# Verify that the fragment files have the correct number of positions.
# In rare circumstances, fragment files were truncated at an arbitrary position.
my $sequence_len = length($sequence);
foreach my $fragment_len (@fragsizes) {
	my $fragment_file = "$options{runid}.$options{n_frags}.$fragment_len" . "mers";

	my $expected_positions = $sequence_len - $fragment_len + 1;
	my $actual_positions = 0;
	open FILE, $fragment_file or die $!;
	while (<FILE>) {
		if ($_ =~ "position:") {
			$actual_positions++;
		}
	}
	close(FILE);
	die "Expected $expected_positions, found $actual_positions" if ($expected_positions != $actual_positions);
}
print "Fragment files have correct number of positions\n";

if ( $options{old_name_format} ) {
    foreach my $size (@fragsizes) {
        my $resultfile = "$options{runid}.$options{n_frags}.$size" . "mers";
        my $fragsize = sprintf( "%2.2d", $size );
        my $fragname = "aa$options{runid}$fragsize\_05.$options{n_frags}\_v1_3";

   # just in case scripts use rsync -a make the old_name_format the actual file
   # and then create a link for the original name (for parallel job checkpoints)
        if ( !-s $fragname ) {
            print_debug("mv $resultfile $fragname");
            ( system("mv $resultfile $fragname") == 0 && -s $fragname )
              or die "ERROR! attempt to move $resultfile to $fragname failed\n";
            print_debug("ln -s $fragname $resultfile");
            ( system("ln -s $fragname $resultfile") == 0 && -s $resultfile )
              or die "ERROR! attempt to link $fragname to $resultfile failed\n";
        }
    }
}

# create a success "checkpoint" file if all fragments are generated
my $success = 1;
foreach my $size (@fragsizes) {
    my $resultfile = "$options{runid}.$options{n_frags}.$size" . "mers";
    $success = 0 if ( !-s $resultfile );
}
system("touch $options{runid}.make_fragments.success") if ($success);

# done
print_debug("done!");
exit 0;

# subroutines

# getCommandLineOptions()
# rets: \%opts  pointer to hash of kv pairs of command line options
sub getCommandLineOptions {
    use Getopt::Long;

    my $usage =
      qq{usage: $0  [-rundir <full path to run directory (default = ./)>
		\t-verbose  specify for chatty output
		\t-id  <5 character pdb code/chain id, e.g. 1tum_>
		\t-csbuild_profile Generate sequence profile via csbuild, rather than psiblast.
		\t-nopsipred  don\'t run psipred. (run by default)
		\t-sam  run sam.
		\t-porter  run porter.
		\t-nohoms  specify to omit homologs from search
		\t-exclude_homologs_by_pdb_date <mm/dd/yy or yyyymmdd>
		\t-psipredfile  <path to file containing psipred ss prediction>
		\t-samfile  <path to file containing sam ss prediction>
		\t-porterfile  <path to file containing porter ss prediction>
		\t-nocleanup  specify to keep all the temporary files
		\t-add_vall_files <vall1,vall2,...> add extra Vall files
		\t-use_vall_files <vall1,vall2,...> use only the following Vall files
		\t-add_pdbs_to_vall <codechain,codechain,...> add extra pdb Vall files
		\t-frag_sizes <size1,size2,...n>
		\t-n_frags <number of fragments>
		\t-n_candidates <number of candidates>
		\t-nofrags specify to not make fragments but run everything else
		\t-old_name_format  use old name format e.g. aa1tum_05.200_v1_3
		\t<fasta file>
	};
    $usage = "$usage\n\n" . join ' ', ( 'Version:', VERSION, "\n" );

    # Get args
    my %opts;
    &GetOptions(
        \%opts,             "psipred!", "csbuild_profile!",
        "sam!",             "homs!",
        "frags!",           "porter!",
        "psipredfile=s",    "samfile=s",
        "porterfile=s",     "verbose!",
        "rundir=s",         "id=s",
        "cleanup!",         "frag_sizes=s",
        "n_frags=i",        "n_candidates=i",
        "add_vall_files=s", "use_vall_files=s",
        "torsion_bin=s",    "exclude_homologs_by_pdb_date=s",
        "old_name_format!", "add_pdbs_to_vall=s"
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
        die "Fragment sizes are invalid\n";
    }

    if (wantarray) { return %opts; }
    return \%opts;
}

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

sub nonempty_file_exists {
    my $file = shift;
    return ( -f $file && !-z $file );
}

sub run {
    my ( $cmd, @files ) = @_;
    my $pid;

    print_debug("cmd is: $cmd");

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

            # child process
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
            if ( !&nonempty_file_exists($f) ) {
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
    my ($runid) = @_;
    my @hits;
    my $hit;
    my $hit_pdb;
    my $hit_chain;
    my %uniq_hits;

    open( EXCL,  ">$runid.homolog_vall" );
    open( EXCL2, ">>$runid.homolog" );

    print_debug("exclude_outn: excluding homologs from $runid.outn");
    open( OUTN, "$runid.outn" );
    @hits = grep( /^>/, <OUTN> );
    close(OUTN);

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

        print EXCL "$hit_pdb$hit_chain\n";
        print EXCL2 "$hit_pdb$hit_chain\n" if (!exists $excluded_homologs{$hit_pdb.$hit_chain});
	$excluded_homologs{$hit_pdb.$hit_chain} = 1;
    }

    close(EXCL);
    close(EXCL2);
}

sub exclude_pdbblast {
    my ($runid) = @_;
    my @hits;
    my $hit;
    my $hit_pdb;
    my $hit_chain;
    my %uniq_hits;

    open( EXCL,  ">$runid.homolog_pdb" );
    open( EXCL2, ">>$runid.homolog" );

    print_debug("exclude_pdbblast: excluding homologs from $runid.outn.pdb");
    open( OUTN, "$runid.outn.pdb" );
    @hits = grep( /^>/, <OUTN> );
    close(OUTN);

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

        print EXCL "$hit_pdb$hit_chain\n";
        print EXCL2 "$hit_pdb$hit_chain\n" if (!exists $excluded_homologs{$hit_pdb.$hit_chain});
	$excluded_homologs{$hit_pdb.$hit_chain} = 1;
    }

    close(EXCL);
    close(EXCL2);
}

sub exclude_blast {
    my ($runid) = @_;
    my @hits;
    my $hit;
    my $hit_pdb;
    my $hit_chain;
    my %uniq_hits;

    open( EXCL,  ">$runid.homolog_nr" );
    open( EXCL2, ">>$runid.homolog" );

    print_debug("exclude_blast: excluding homologs from $runid.blast");
    open( BLAST, "$runid.blast" );
    @hits = grep( /pdb\|/, <BLAST> );
    close(BLAST);

    foreach $hit (@hits) {
        $uniq_hits{$hit} = 1;
    }

    foreach $hit ( sort keys %uniq_hits ) {
        while ( $hit =~ s/pdb\|(\w{4})\|(.?)// ) {
            $hit_pdb   = $1;
            $hit_chain = $2;

            $hit_pdb =~ tr/[A-Z]/[a-z]/;
            $hit_chain = '_' if ( $hit_chain =~ /^\s*$/ );

            print EXCL "$hit_pdb$hit_chain\n";
            print EXCL2 "$hit_pdb$hit_chain\n" if (!exists $excluded_homologs{$hit_pdb.$hit_chain});
	    $excluded_homologs{$hit_pdb.$hit_chain} = 1;
	}
    }

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

    open( INPUT, $filename ) or die("Couldn't open $filename for reading.\n");

    read( INPUT, $buf, 4 ) or die("Couldn't read $filename!\n");
    $seqlen = unpack( "i", $buf );

    print_debug("Sequence length: $seqlen");

    read( INPUT, $buf, $seqlen ) or die("Premature end: $filename.\n");
    $seqstr = unpack( "a$seqlen", $buf );

    for ( $i = 0 ; $i < $seqlen ; ++$i ) {
        read( INPUT, $buf, 160 ) or die("Premature end: $filename, line: $i\n");
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
        X => 0
    );

    ( length($s) == scalar(@matrix) )
      or die("Length mismatch between sequence and checkpoint file!\n");

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
    foreach my $blosumline (@blosum62) {
        next if ( $blosumline !~ /^\d/ );
        chomp $blosumline;
        @words = split( /\s/, $blosumline );

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

sub is_sam_6state {
	my $file = shift;
	chomp $file;
	my $returnval = 0;
	open(F, $file) or die "ERROR! cannot open $file: $!\n";
	while (<F>) {
		if (/^[0-9]/) {
			my @parts = split;
			if ($#parts == 7) {
				$returnval = 1;
				last;
			}
		}
	}
	return $returnval;
}


sub convert_sam_ss {
    my $ss_pred = shift;
    my ( $e, $h, $c, @ss );
    my $plusOne = 0;
    open( SAM, $ss_pred ) or die "ERROR! cannot open $ss_pred: $!\n";
    open( SAMCONV, ">$ss_pred\_ss2" )
      or die "ERROR! cannot open $ss_pred\_ss2: $!\n";
    print SAMCONV "# PSIPRED VFORMAT (PSIPRED V2.6 by David Jones)\n\n";
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
              $cols[0] + $plusOne, $c, $h, $e;
        }
        elsif ( $h > $e && $h >= $c ) {
            printf SAMCONV "%4d $cols[1] H   %5.3f  %5.3f  %5.3f\n",
              $cols[0] + $plusOne, $c, $h, $e;
        }
        elsif ( $c > $e && $c > $h ) {
            printf SAMCONV "%4d $cols[1] C   %5.3f  %5.3f  %5.3f\n",
              $cols[0] + $plusOne, $c, $h, $e;
        }
    }
    close(SAMCONV);
    close(SAM);
}

sub print_debug {
    if ( $options{DEBUG} ) {
        foreach my $line (@_) {
            print "$line\n";
        }
    }
}

sub cleanup {
    if ( $options{cleanup} ) {
        foreach my $file (@cleanup_files) {
            print_debug("removing file $file");
            unlink($file);
        }
    }
}

sub options_to_str {
    my $options = shift;

    my $str = '';

    use List::Util qw/ max /;
    my $max_width = max map { length($_) } keys %$options;
    my $format_str = "%" . $max_width . "s";
    foreach my $key ( sort keys %$options ) {
        my $val_str = $options->{$key};
        if ( ref( $options->{$key} ) eq 'ARRAY' ) {
            $val_str = join ' ', @{ $options->{$key} };
        }
        elsif ( ref( $options->{$key} ) eq 'HASH' ) {
            my @keys = keys %{ $options->{$key} };
            $val_str = join ', ',
              map { join ': ', ( $_, $options->{$key}{$_} ) } @keys;
        }
        $str .= join ': ', ( sprintf( $format_str, $key ), $val_str );
        $str .= "\n";
    }
    $str .= '-' x 80 . "\n";
    return $str;
}

sub produce_output_with_cmd {
    my $cmd          = shift;
    my $output_fn    = shift;
    my $no_output_ok = 0;
    $no_output_ok = shift if ( scalar @_ );

    print_debug("Command: $cmd");
    if ( -f $output_fn ) {
        print_debug("Skipping command, $output_fn exists!");
        return;
    }

    my $output = `$cmd`;
    print_debug("Finished running command: $cmd");
    if ( defined $output ) {
        print_debug("Output: $output");
    }

    if ( !$no_output_ok && !-f $output_fn ) {
        die "Error: expected creation of $output_fn after running '$cmd'!\n";
    }
}

sub read_fasta {
    my $fn = shift;

    open( SEQFILE, $fn ) or die "Error opening fasta file $fn!\n";
    my $has_comment = 0;
    my $eof;
    my $seq = '';
    while ( my $line = <SEQFILE> ) {
        $eof = $line;
        $line =~ s/\s//g;
        if ( $line =~ /^>/ ) {
            $has_comment = 1;
        }
        else {
            chomp $line;
            $seq .= $line;
        }
    }

    my $has_eof = ( $eof =~ /\n$/ );
    ( $has_comment && $has_eof )
      or die
      "fasta file (given $fn) must have a comment and end with a new line!\n";
    close SEQFILE or die $!;

    return $seq;
}

sub get_ss_quota_vote {
    my $ss_pred_name = shift;

    if ( $ss_pred_name eq 'psipred' ) {
        return 3;
    }
    elsif ( $ss_pred_name eq 'sam' ) {
        return 1;
    }
    elsif ( $ss_pred_name eq 'porter' ) {
        return 1;
    }
    else {
        die "Error: don't recognize ss_pred_name $ss_pred_name!\n";
    }
}

sub file_overrides_option {
    my $fn  = shift;
    my $tag = shift;

    if ( nonempty_file_exists($fn) ) {
        my $fn_base = basename($fn);
        if ( abs_path("$options{rundir}/$fn_base") ne abs_path($fn) ) {
            copy( $fn, "$options{rundir}/$fn_base" )
              or die "Error copying $fn -> $options{rundir}/$fn_base";
            push( @cleanup_files, $fn_base );
        }

        print "Assuming $options{rundir}/$fn_base is a $tag file.\n";
        return 0;
    }
    print "File for $tag not found! Generating from scratch instead.\n";
    if ( length($fn) ) {
        print "(given fn $fn)\n";
    }
    return 1;
}

sub run_in_parallel {
    my $commands_array_ref    = shift;
    my $resultfiles_array_ref = shift;
    my $runid                 = shift;
    my $maxparalleljobs       = shift;
    my $maxwait               = shift;
    my $maxattempts           = shift;
    my $ignore_errors         = 0;
    $ignore_errors = shift if ( scalar @_ );
    my @commands = @$commands_array_ref;
    my @results  = @$resultfiles_array_ref;
    return if ( !scalar @commands );
    ( scalar @commands == scalar @results )
      or die
      "ERROR! number of commands does not match results in run_in_parallel\n";

    #Only load ForkManager if required.
		eval "use ForkManager";
		use Sys::Hostname;
		my $host = hostname;
    my $logprefix = "$runid-$host";
    eval {
        local $SIG{ALRM} = sub { die "alarm\n" };
        alarm $maxwait;
        foreach my $attempt ( 1 .. $maxattempts ) {
            print_debug("Attempt $attempt");
            my $pm = Parallel::ForkManager->new($maxparalleljobs);
            for ( my $i = 0 ; $i <= $#commands ; ++$i ) {
                my $resultfile = $results[$i];
                next if ( -s $resultfile );
                $pm->start and next;
                my $PID     = $$;
                my $logfile = "$logprefix\_$i.log";
                my $command = "$commands[$i] >& $logfile";
                eval {
                    alarm $maxwait;
                    print_debug("starting process $PID to make $resultfile.");
                    system($command)
                      ; # keep in mind the command can be orphaned if this is timed out
                    print_debug("finished process $PID.");
                    alarm 0;
                };
                if ($@) { alarm 0; }
                $pm->finish;
            }
            $pm->wait_all_children;

            my %checkdone;
            for ( my $i = 0 ; $i <= $#results ; ++$i ) {
                my $resultfile = $results[$i];
                if ( -s $resultfile ) {
                    $checkdone{$resultfile} = 1;
                    print_debug("$resultfile exists!");
                }
                else {
                    print_debug("$resultfile missing!");
                }
            }
            last if ( scalar( keys %checkdone ) == scalar(@results) );
            if ( $attempt >= $maxattempts ) {
                if ($ignore_errors) {
                    warn
"WARNING! all result files do not exist after attempt $attempt.\n";
                    last;
                }
                else {
                    die
                      "ERROR! failed run_in_parallel after $attempt attempts\n";
                }
            }
            warn
"WARNING! all result files do not exist after attempt $attempt, attempting run_in_parallel again\n";
        }
        alarm 0;
    };
    if ($@) {
        alarm 0;
        if ($ignore_errors) {
            warn "WARNING! timed out run_in_parallel!\n";
        }
        else {
            die "ERROR! timed out run_in_parallel\n";
        }
    }
}

sub is_ss2_format {
	my $file = shift;
	chomp $file;
	my $returnval = 0;
	open(F, $file) or do {
		warn "WARNING! cannot open $file: $!\n";
		return 0;
	};
	while (<F>) {
		if (/^# PSIPRED VFORMAT/) {
			$returnval = 1;
			last;
		}
	}
	close(F);
	return $returnval;
}


