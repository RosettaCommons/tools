#!/usr/bin/perl -w
##
##
## Copyright 2000-2008, University of Washington, All rights reserved
## written by Dylan Chivian, Department of Biochemistry.
## use of software governed under the BSD license.
## software can also be obtained via sourceforge (url if available).
##
##
##  Initial Author: Dylan Chivian (dylan@lazy8.com)
##  $Revision: 21627 $
##  $Date: 2008-04-07 22:03:44 +0300 (Mon, 07 Apr 2008) $
##  $Author: yiliu $
##
##
###############################################################################

###############################################################################
package PDButil;
###############################################################################

###############################################################################
# init
###############################################################################

$| = 1;                                                   # don't buffer stdout

local %opts = &getCommandLineOptions ();
local $pdbfile = $opts{pdbfile};
local $chainOI = $opts{chain};

$pdbID = $pdbfile;
$pdbID =~ s!^.*/?p?d?b?(\w\w\w\w)\.[pe][dn][bt]\.?g?z?Z?$!$1!;
$pdbID = lc $pdbID;

###############################################################################
# main
###############################################################################

# read
#
@buf = &fileBufArray ($pdbfile);

for ($i=0; $i <= $#buf; ++$i) {
    if ($buf[$i] =~ /^ATOM/) {
	$chain    = substr ($buf[$i], 21, 1);
	$atomtype = substr ($buf[$i], 13, 2);
	if ($chain eq $chainOI && $atomtype eq 'CA') {
	    $x = substr ($buf[$i], 30, 8);
	    $y = substr ($buf[$i], 38, 8);
	    $z = substr ($buf[$i], 46, 8);
	    printf ("%8.3f %8.3f %8.3f\n", $x, $y, $z);
	}
    }
}

# exit
#
exit 0;

###############################################################################
# subs
###############################################################################

# getCommandLineOptions()
#
#  desc: get the command line options
#
#  args: none
#
#  rets: \%opts  pointer to hash of kv pairs of command line options
#
sub getCommandLineOptions {
    use Getopt::Long;
    local $usage = qq{usage: $0 -pdbfile <pdbfile> -chain <chain>};

    # Get args
    #
    local %opts = ();
    &GetOptions (\%opts, "pdbfile=s", "chain=s");


    # Check for legal invocation
    #
    if (! defined $opts{pdbfile}) {
        print STDERR "$usage\n";
        exit -1;
    }
    &checkExistence ('f', $opts{pdbfile});	


    # Defaults
    #
    $opts{chain} = ' '  if (! defined $opts{chain} || $opts{chain} eq '0' || $opts{chain} eq '_');
    $opts{chain} = uc $opts{chain};
 
    return %opts;
}

###############################################################################
# util
###############################################################################

sub maxInt {
    local ($v1, $v2) = @_;
    return ($v1 > $v2) ? $v1 : $v2;
}

sub tidyDecimals {
    my ($num, $decimal_places) = @_;
    if ($num !~ /\./) {
	$num .= '.' . '0' x $decimal_places;
	$num =~ s/^0+//;
    }
    else {
	if ($num =~ s/(.*\.\d{$decimal_places})(\d).*$/$1/) {
	    my $nextbit = $2;
	    if ($nextbit >= 5) {
		my $flip = '0.' . '0' x ($decimal_places - 1) . '1'; 
		$num += $flip;
	    }
        }
	$num =~ s/^0//;
	my $extra_places = ($decimal_places + 1) - length $num;
	$num .= '0' x $extra_places  if ($extra_places > 0);
    }

    return $num;
}

sub distsq {
    local @dims = @_;
    local $v = 0;
    foreach $dim (@dims) {
	$v += $dim*$dim;
    }
    return $v;
}

sub logMsg {
    local ($msg, $logfile) = @_;

    if ($logfile) {
        open (LOGFILE, ">".$logfile);
        select (LOGFILE);
    }
    print $msg, "\n";
    if ($logfile) {
        close (LOGFILE);
        select (STDOUT);
    }
    return 'true';
}

sub checkExistence {
    local ($type, $path) = @_;
    if ($type eq 'd') {
	if (! -d $path) { 
            print STDERR "$0: dirnotfound: $path\n";
            exit -3;
	}
    }
    elsif ($type eq 'f') {
	if (! -f $path) {
            print STDERR "$0: filenotfound: $path\n";
            exit -3;
	}
    }
}

sub abort {
    local $msg = shift;
    print STDERR "$0: $msg\n";
    exit -2;
}

sub writeBufToFile {
    ($file, $bufptr) = @_;
    if (! open (FILE, '>'.$file)) {
	&abort ("$0: unable to open file $file for writing");
    }
    print FILE join ("\n", @{$bufptr}), "\n";
    close (FILE);
    return;
}

sub fileBufString {
    local $file = shift;
    local $oldsep = $/;
    undef $/;
    if ($file =~ /\.gz|\.Z/) {
	if (! open (FILE, "gzip -dc $file |")) {
	    &abort ("$0: unable to open file $file for gzip -dc");
	}
    }
    elsif (! open (FILE, $file)) {
	&abort ("$0: unable to open file $file for reading");
    }
    local $buf = <FILE>;
    close (FILE);
    $/ = $oldsep;
    return $buf;
}

sub fileBufArray {
    local $file = shift;
    local $oldsep = $/;
    undef $/;
    if ($file =~ /\.gz|\.Z/) {
	if (! open (FILE, "gzip -dc $file |")) {
	    &abort ("$0: unable to open file $file for gzip -dc");
	}
    }
    elsif (! open (FILE, $file)) {
	&abort ("$0: unable to open file $file for reading");
    }
    local $buf = <FILE>;
    close (FILE);
    $/ = $oldsep;
    @buf = split (/$oldsep/, $buf);
    pop (@buf)  if ($buf[$#buf] eq '');
    return @buf;
}

###############################################################################
1; # package end
# end
###############################################################################
