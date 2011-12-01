#!/usr/bin/perl
##
##Copyright 2000-2008, University of Washington, All rights reserved
## written by Dylan Chivian, Department of Biochemistry.
## Use of software governed under the BSD license.
## Software can also be obtained via sourceforge (url if available)
##  Initial Author: Dylan Chivian (dylan@lazy8.com)
##  $Revision: 21627 $
##  $Date: 2008-04-07 22:03:44 +0300 (Mon, 07 Apr 2008) $
##  $Author: yiliu $
##
###############################################################################
                   

# init
if ($#ARGV < 1) {
    print STDERR "usage: $0 <pdb> <chain>\n";
    exit -1;
}
$rosetta_pdb = shift @ARGV;
$chain = shift @ARGV;
$chain = ($chain eq '_') ? ' ' : uc $chain;

if (! -f $rosetta_pdb) {
    print STDERR "$0: filenotfound: $rosetta_pdb\n";
    exit -1;
}
@pdb_buf = &fileBufArray ($rosetta_pdb);


# body
for ($i=0; $i <= $#pdb_buf; ++$i) {
    if ($pdb_buf[$i] =~ /^ATOM|^HETATM/) {
	substr ($pdb_buf[$i], 21, 1) = $chain;
    }
    print $pdb_buf[$i] , "\n";
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
