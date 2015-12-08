#!/usr/bin/perl

# Ok, so this script was pretty much copied from the zapChain.pl script, and
# modified so that all hydrogen entries are zapped.

# init
if ($#ARGV < 0) {
    print STDERR "usage: $0 <pdb> \n";
    print STDERR "removes all H (hydrogen) entries\n";
    exit -1;
}
$rosetta_pdb = shift @ARGV;

if (! -f $rosetta_pdb) {
    print STDERR "$0: filenotfound: $rosetta_pdb\n";
    exit -1;
}
@pdb_buf = &fileBufArray ($rosetta_pdb);


# body
for ($i=0; $i <= $#pdb_buf; ++$i) {
    if ($pdb_buf[$i] =~ /^ATOM|^HETATM/) {
#       print substr ($pdb_buf[$i], 12, 1),$chain;
	# if the 13th character is not an H, print the line to output.
        if (! (substr ($pdb_buf[$i], 13, 1) eq 'H')) {
          print $pdb_buf[$i] , "\n";
        }
    }
    else {
        print $pdb_buf[$i] , "\n";
    }
}

exit 0;


###############################################################################
# util
###############################################################################

# fileBufString()
#
sub fileBufString {
    local $file = shift;
    local $oldsep = $/;
    undef $/;
    if (! open (FILE, $file)) {
        &abort ("$0: unable to open file $file for reading");
    }
    local $buf = <FILE>;
    close (FILE);
    $/ = $oldsep;
    return $buf;
}

# fileBufArray()
#
sub fileBufArray {
    local $file = shift;
    local $oldsep = $/;
    undef $/;
    if (! open (FILE, $file)) {
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
# end
###############################################################################
