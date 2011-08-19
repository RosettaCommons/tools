#!/usr/bin/perl
##
## Copyright 2000-2008, University of Washington, All rights reserved
## written by Dylan Chivian, Department of Biochemistry.
## use of software governed under the BSD license.
## software can also be obtained via sourceforge (url if available).
##
##  Initial Author: Dylan Chivian (dylan@lazy8.com)
##  $Revision: 21627 $
##  $Date: 2008-04-07 22:03:44 +0300 (Mon, 07 Apr 2008) $
##  $Author: yiliu $
##
##
###############################################################################

###############################################################################
# init
###############################################################################

local %opts = &getCommandLineOptions ();

###############################################################################
# paths
###############################################################################

$rcsb_server     = "http://www.rcsb.org";
$rcsb_base_uri   = "/pdb/cgi/export.cgi/";
$pdbobs_server   = "http://pdbobs.sdsc.edu";
$pdbobs_base_uri = "/all_entries";

###############################################################################
# conf
###############################################################################

local $id = lc $opts{'id'};
local ($c1, $c2, $c3, $c4, $underscore, $chain) = split (//, $id);

$pdb_base = lc ($c1.$c2.$c3.$c4);

###############################################################################
# main
###############################################################################

# go to outdir
#
if ($opts{'outdir'}) {
    $outdir = $opts{'outdir'};
    if (! -d $outdir) {
	if (system (qq{mkdir -p $outdir}) != 0) {
	    &abort ("failure creating $outdir");
	}
    }
    chdir $outdir;
}
if (-f "$pdb_base.pdb.Z") {
    print "$0: already have $pdb_base.pdb.Z\nexiting...\n";
    exit 0;
}


# http get it
#
$uc_pdb_base_id = uc $pdb_base;
$url = "$rcsb_server$rcsb_base_uri$uc_pdb_base_id.pdb.Z\?format=PDB\&pdbId=$uc_pdb_base_id\&compression=Z";
$pdb_data = &httpget ($url);


# check the data
#
#if ($pdb_data =~ /^\s*$/) {
if ($pdb_data =~ /^\s*$/ || $pdb_data =~ /404\s+not\s+found/i) {
    print STDERR ("$0: no such pdb $pdb_base.pdb.Z found at $url\n");

    $url = "$pdbobs_server$pdbobs_base_uri/pdb$pdb_base.ent";
    $pdb_data = &httpget ($url);
    if ($pdb_data =~ /^\s*$/ || $pdb_data =~ /404\s+not\s+found/i) {
	&abort ("no such pdb pdb$pdb_base.ent found at $url");
    } else {
	$pdb_src = 'pdbobs';
    }
} else {
    $pdb_src = 'rcsb';
}


# write it
#
if ($pdb_src eq 'rcsb') {
    open (PDB, '>'."$pdb_base.pdb.Z");
    print PDB $pdb_data;
    close (PDB);
} elsif ($pdb_src eq 'pdbobs') {
    open (PDB, '>'."$pdb_base.pdb");
    print PDB $pdb_data;
    close (PDB);
#    system (qq{compress -f $pdb_base.pdb});
    system (qq{gzip -9f -S .Z $pdb_base.pdb});
}


# double check to make sure we made it
#
if (! -s "$pdb_base.pdb.Z") {
    &abort ("failed to create $pdb_base.pdb.Z");

}

# done
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
    local $fail = 'FALSE';
    local $usage = "usage: $0 -id <pdb_id (nochain)> [-outdir <outdir>]";

    # Get args
    #
    local %opts = ();
    &GetOptions (\%opts, "id=s", "outdir=s");

    # Check for legal invocation
    #
    if (!defined $opts{id}) {
        $fail = 'TRUE';
    }
    if ($fail eq 'TRUE') {
        print STDERR "$usage\n";
        exit -1;
    }

    return %opts;
}


# httpget()
#
sub httpget {
    my $showhead = 0;
    my $ret_msg = '';
    ($host, $port, $uri) = &getGetParams (@_);
    ($stat, $head, $page) = &SiteSucker::makeHttpRequest ($host, $port, $uri); 
    
    if ($showhead) {
	foreach $k (sort keys %$head) {
	    $ret_msg .= $head->{$k}."\n";
	}
	$ret_msg .= "\n\n";
    }

    $ret_msg .= "$page\n";
    return $ret_msg;
}


# httppost()
#
sub httppost {
    my $showhead = 0;
    my $ret_msg = '';
    ($host, $port, $uri, $postfile, $flat_or_multipart) = &getPostParams (@_);
    $post = &fileBufString ($postfile);
    ($stat, $head, $page) = &SiteSucker::makeHttpRequestPost ($host, $port, $uri, $post, $flat_or_multipart); 
    
    if ($showhead) {
	foreach $k (sort keys %$head) {
	    $ret_msg .= $head->{$k}."\n";
	}
	$ret_msg .= "\n\n";
    }
    $ret_msg .= "$page\n";
    return $ret_msg;
}


# getGetParams()
#
sub getGetParams {
    my @argv = @_;
    my ($host, $port, $uri);

    if ($#argv < 0) {
        print "usage: $0 <url>\n";
        exit 0;
    }
    $url = $argv[0];;

    $url =~ s!^http://!!;
    if ($url =~ m!^([^:/]+):?(\d*)(.*)!) {
        $host = $1;
        $port = $2;
        $uri = $3;
    }
    else {
        print STDERR "$0: malformed url '$url'\n";
        exit -2;
    }

    $port = 80  if (! $port);
    $uri = '/'  if (! $uri);

    #$uri .= '/'  if ($uri !~ /\/$/ && $uri !~ /\./ && $uri !~ /\?/);

    return ($host, $port, $uri);
}


# getPostParams
#
sub getPostParams {
    my @argv = @_;
    my ($host, $port, $uri);

    if ($#argv < 0) {
	print "usage: $0 <url> <postfile> <flat_or_multipart>\n";
	exit 0;
    }
    $url               = $argv[0];
    $postfile          = $argv[1];
    $flat_or_multipart = $argv[2];

    $url =~ s!^http://!!;
    if ($url =~ m!^([^:/]+):?(\d*)(.*)!) {
	$host = $1;
	$port = $2;
	$uri  = $3;
    }
    else {
	print STDERR "$0: malformed url '$url'\n";
	exit -2;
    }

    $port = 80  if (! $port);
    $uri = '/'  if (! $uri);

    #$uri .= '/'  if ($uri !~ /\/$/ && $uri !~ /\./ && $uri !~ /\?/);

    return ($host, $port, $uri, $postfile, $flat_or_multipart);
}


###############################################################################
# util
###############################################################################

# maxInt()
#
sub maxInt {
    local ($v1, $v2) = @_;
    return ($v1 > $v2) ? $v1 : $v2;
}
# end maxInt()


# makeGetArgs (\%a)
#
#
#   desc:  make a string suitable for GET args out of associative array
#
#   args:  \%a    array with key->val pairs
#
#   rets:  $ret   string for GET args
#
#
sub makeGetArgs {
    local ($a) = @_;
    my ($k, $s, $ret);

    foreach $k (keys %$a) {
        $s = $a->{$k};
        $s = &escapeGetArg ($s);
        $ret .= "$k=$s&";
    }
    chop $ret;                                      # Remove trailing ampersand
    return $ret;
}
# end makeGetArgs ()


# escapeGetArg ()
#
sub escapeGetArg {
    local $str = shift;

    $str =~ s/ /+/go;
    $str =~ s/\%/\%25/go;
    $str =~ s/([^0-9a-zA-Z\%+])/"\%".&charToHex($1)/ge;

    return $str;
}
# end escapeGetArg ()


# charToHex ()
#
sub charToHex {
    my $ascii = ord($_[0]);
    my %hexMap = (  0 => '0',
                    1 => '1',
                    2 => '2',
                    3 => '3',
                    4 => '4',
                    5 => '5',
                    6 => '6',
                    7 => '7',
                    8 => '8',
                    9 => '9',
                   10 => 'a',
                   11 => 'b',
                   12 => 'c',
                   13 => 'd',
                   14 => 'e',
                   15 => 'f'
		    );

    return $hexMap{(($ascii & 0xf0) >> 4)} . $hexMap{($ascii & 0x0f)};
}
# end charToHex ()


#  hexToChar ()
#
sub hexToChar {
    my $ascii = hex($_[0]);
    return chr $ascii;
}
# end hexToChar ()


# checkExist()
#
sub checkExist {
    local ($type, $path) = @_;
    if ($type eq 'f' && ! -f $path) {
        &abort ("filenotfound $path");

    }
    if ($type eq 'd' && ! -d $path) {
        &abort ("dirnotfound $path");
    }
    return 'true';
}


# checkExistAndCreate()
#
sub checkExistAndCreate {
    local ($type, $path) = @_;
    if ($type eq 'f' && ! -f $path) {
        print "creating $path...\n";
        open (FILE, '>'.$path);
        close (FILE);
    }
    if ($type eq 'd' && ! -d $path) {
        print "creating $path...\n";
        $mode = 0777  if (! $mode);
        mkdir ($path, $mode);
    }
    return 'true';
}

# abort()
#
sub abort {
    local $msg = shift;
    print STDERR "$0: $msg\n";
#    print "$0: $msg\n";
    exit -2;
}


# writeBuf ()
#
sub writeBuf {
    my ($buf, $outfile) = @_;
    if ($outfile) {
        #print "creating $outfile\n";
        open (OUTFILE, '>'.$outfile);
        select (OUTFILE);
    }
    print $buf;
    if ($outfile) {
        close (OUTFILE);
        select (STDOUT);
    }
    return;
}


# fileBufString ()
#
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


# fileBufArray ()
#
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
package SiteSucker;
###############################################################################

use Socket;
use FileHandle;


$accept_cookies = 0;
$keep_sid = 0;
$cookie = '';
$debug = 0;

$| = 1  if ($debug);
%uriVisited = ();


sub crawlURL {
    my ($host, $port, $uri, $refuri, $sid) = @_;
    my ($req_uri, $clean_uri, $cwd_uri) = ('','','');
    my %snd_hdrs = ();
    my $head = undef;
    my $page = '';
    my @links = ();
    my ($new_host, $new_port, $link, $path_elt) = ('','','','');
    my @new_uri = ();
    my $nxt_uri = '';
    my $head_end = -1;
    my $set_cookie = '';
    

    if ($uri !~ m!/cgi-bin/!) {
	$uri =~ s!(/[^\./]+)$!$1/!;
    }
    $clean_uri = $uri;
    $clean_uri =~ s!\#[^/]+$!!;
    $clean_uri =~ s!\?[^/]+$!!;
    return  if ($uriVisited{$clean_uri});
    $uriVisited{$clean_uri} = 1;

    $cwd_uri = $uri;
    $cwd_uri =~ s![^/]*$!!;
    $cwd_uri = '/'  if (! $cwd_uri);


    # Read page
    #
    $req_uri = $uri;
    if ($cookie) {
	$snd_hdrs{Cookie} = $cookie;
    }
    else {
	if (! $accept_cookies && $keep_sid && $sid) {
	    $req_uri =  &insertSession ($uri, $sid) ;
	}
    }
    ($stat, $head, $page) = &makeHttpRequest($host,$port,$req_uri,%snd_hdrs);

    if ($stat != 1) {
	print STDERR "$0: failure: unable to makeHttpRequest()\n";
	exit 2;
    }


    # Play with mime and page
    #
    if (defined &main::examineMime) {
	&main::examineMime ($host, $port, $uri, $refuri, $head);
    }
    if (defined &main::examinePage) {    
	&main::examinePage ($host, $port, $uri, $refuri, $page);
    }


    # Crawl links
    #
#    if (($head_end = index ($page, "</head>")) == -1) {
#	$head_end = index ($page, "</HEAD>");
#    }
#    if ($head_end != -1) {
#	$head_end += length ("</head>");
#	$page = substr ($page, $head_end, length ($page) - $head_end);
#    }
    while ($page =~ s/<a[^<>]+href\s*=\s*\"?([^\"\'\s>]+)[^<>]*>//im) {
	push (@links, $1);
    }
    while ($page =~ s/<frame[^<>]+src\s*=\s*\"?([^\"\'\s>]+)[^<>]*>//im) {
	push (@links, $1);
    }
    undef $page;                                                  # save memory
    undef $head;


    foreach $link (@links) {
	next  if ($link =~ /^mailto:/i);
	next  if ($link =~ /^ftp:/i);
	next  if ($link =~ /^javascript:/i);

	next  if ($link =~ /^\#/);
	$link =~ s/\#.*$//;
	if (! $accept_cookies && $keep_sid) {
	    $sid = &getSession ($link);
	}
	$link = &removeSession ($link);

	if ($link =~ s!^http://([^\:/]+)!!i) {
	    $new_host = $1;
	    next  if ($new_host ne $host);
	}

	if ($link =~ s!^:(\d+)!!) {
	    $new_port = $1;
	    next  if ($new_port ne $port);
	}

	if ($link =~ m!^/!) {
	    $nxt_uri = $link;
	}
	else {
	    $cwd_uri =~ s!^/|/$!!g;
	    @new_uri = split (/\/{1,2}/, $cwd_uri);
	    
	    foreach $path_elt (split (/\/{1,2}/, $link)) {
		next if ($path_elt eq '.' || $path_elt eq '');
		if ($path_elt eq '..') {
		    pop (@new_uri);
		}
		else {
		    push (@new_uri, $path_elt);
		}
	    }
	    $nxt_uri = '/' . join ('/', @new_uri);
	    $nxt_uri .= '/'  if ($new_uri[$#new_uri] !~ /\./);
	}
	&crawlURL ($host, $port, $nxt_uri, $uri, $sid);
    }

    return;
}


sub makeHttpRequest {
    my ($host, $port, $uri, %snd_hdrs) = @_;
    my ($nic);
    my ($name, $aliases, $proto, $type, $len, $addr);
    my $head = undef;
    my $data = '';
    my $status = 1;
    my $linesep = $/;

    
    ($name, $aliases, $proto) = getprotobyname('tcp');
    ($name, $aliases, $type, $len, $addr) = gethostbyname($host);
    $nic = pack('S n N x8', AF_INET, $port, unpack("N", $addr));
    
    if (!socket(S, PF_INET, SOCK_STREAM, $proto)) {
	$status = -1;
	goto done;
    }
    if (!connect(S, $nic)) {
	$status = -2;
	goto done;
    }

    S->autoflush(1);
    $request  = "GET $uri HTTP/1.0\r\n";
    foreach $i (keys %snd_hdrs) {
	next if (! $i);
	$request .= "$i: " . $snd_hdrs{$i} . "\r\n";
    }
    $request .= "\r\n";
    print S $request;


    # Get mime header
    #
    $head->{status} = <S>;
    while (<S>) {
	last  if ($_ =~ /^\s*$/);
	$_ =~ /^([^:]+):\s*(.*)/;
	$head->{$1} = $2;
	if ($accept_cookies && $_ =~ /^set-cookie: (.*)/i) {
	    $set_cookie = $1;
	    $set_cookie =~ s/expires=[^\;]+//;
	    $set_cookie =~ s/path=[^\;]+//;
	    $set_cookie =~ s/domain=[^\;]+//;
	    $set_cookie =~ s/secure//;
	    $set_cookie =~ s/\s*\;\s*//g;
	    $cookie .= " $set_cookie\;";
	}
    }
    $cookie =~ s/^\s+|\s*\;$//g;


    # Get data
    #
    undef $/;
    $data = <S>;
    close (S);
    $/ = $linesep;
    
done:
    return ($status, $head, $data);
}


sub makeHttpRequestPost {
    my ($host, $port, $uri, $post, $flat_or_multipart, %snd_hdrs) = @_;
    my ($nic);
    my ($name, $aliases, $proto, $type, $len, $addr);
    my $head = undef;
    my $data = '';
    my $status = 1;
    my $linesep = $/;

    my $boundary;

    if ($flat_or_multipart !~ /^f/i) { 
	($post, $boundary) = &makeMultiPartPost ($post);
    }

    ($name, $aliases, $proto) = getprotobyname('tcp');
    ($name, $aliases, $type, $len, $addr) = gethostbyname($host);
    $nic = pack('S n N x8', AF_INET, $port, unpack("N", $addr));
    
    if (!socket(S, PF_INET, SOCK_STREAM, $proto)) {
	$status = -1;
	goto done;
    }
    if (!connect(S, $nic)) {
	$status = -2;
	goto done;
    }

    S->autoflush(1);
    #$request  = "POST $uri HTTP/1.0\r\n";
    $request  = "POST $uri HTTP/1.1\r\n";
    foreach $i (keys %snd_hdrs) {
	next if (! $i);
	$request .= "$i: " . $snd_hdrs{$i} . "\r\n";
    }
    my $content_length = length ($post);
    $request .= "Connection: " . "close"          . "\r\n";
    #$request .= "Connection: " . "Keep-Alive"     . "\r\n";
    $request .= "User-Agent: " . "HTTPposter/1.0" . "\r\n";
    $request .= "Host: "       . "$host:$port"    . "\r\n";
    $request .= ($flat_or_multipart !~ /^f/i) 
	        ? "Content-Type: " . "multipart/form-data; boundary=$boundary" 
		: "Content-Type: " . "application/x-www-form-urlencoded";
    $request .= "\r\n"; 
    $request .= "Content-Length: " . $content_length . "\r\n";
    $request .= "\r\n";
    $request .= $post . "\r\n";
    print S $request;

    # debug
    #print STDOUT $request;

    # Get mime header
    #
    $head->{status} = <S>;
    while (<S>) {
	last  if ($_ =~ /^\s*$/);
	$_ =~ /^([^:]+):\s*(.*)/;
	$head->{$1} = $2;
	if ($accept_cookies && $_ =~ /^set-cookie: (.*)/i) {
	    $set_cookie = $1;
	    $set_cookie =~ s/expires=[^\;]+//;
	    $set_cookie =~ s/path=[^\;]+//;
	    $set_cookie =~ s/domain=[^\;]+//;
	    $set_cookie =~ s/secure//;
	    $set_cookie =~ s/\s*\;\s*//g;
	    $cookie .= " $set_cookie\;";
	}
    }
    $cookie =~ s/^\s+|\s*\;$//g;


    # Get data
    #
    undef $/;
    $data = <S>;
    close (S);
    $/ = $linesep;
    
done:
    return ($status, $head, $data);
}


sub makeMultiPartPost {
    my $post = shift;
    my ($ret_post, $boundary);
    my ($pair, $k, $v);
    my $rand;

    srand();
    $boundary = '-'x29;
    for (my $i=0; $i < 28; ++$i) {
	$rand = rand 10;
	$rand = int ($rand);
	$boundary .= $rand;
    }

    my @args = split (/\&/, $post);
    foreach $pair (@args) {
	next if ($pair =~ /^\s*$/);
	($k, $v) = split (/=/, $pair);
	$k = &unEscapeGetArg ($k);
	$v = &unEscapeGetArg ($v);
	$ret_post .= qq{$boundary\r\n};
	$ret_post .= qq{Content-Disposition: form-data; name="$k"\r\n\r\n$v\r\n};
    }
    $ret_post .= "$boundary--";

    return ($ret_post, $boundary);
}


sub getSession {
    my $uri = shift;
    $uri =~ /\@SK\@([^\@]+)\@\@/;
    return $1;
}


sub removeSession {
    my $uri = shift;

    $uri =~ s/\@SK\@[^\@]+\@\@//;
    return $uri;
}


sub insertSession {
    my ($uri, $sid) = @_;
    my $new_uri = '';
    my $sid_path = "\@SK\@$sid\@\@";
    my $uri_len = length ($uri);
    my $part1 = '';
    my $part1_len = 0;
    my $marker = 0;
    my ($octothorpe_marker, $questionmark_marker) = (-1,-1);


    if (($octothorpe_marker = index ($uri, "\#")) != -1) {
	$marker = $octothorpe_marker;
    }
    if (($questionmark_marker = index ($uri, "?")) != -1
	&& $questionmark_marker < $marker) {
	$marker = $questionmark_marker;
    }
    if ($marker > 0) {
	--$marker;
    }
    else {
	$marker = $uri_len - 1;
    }


    if (substr ($uri, $marker, 1) eq '/') {
	$new_uri = substr ($uri, 0, $marker + 1);
	$new_uri .= $sid_path;
	++$marker;
	$new_uri .= substr ($uri, $marker)  if ($marker < $uri_len);
    }
    else {
	$part1 = substr ($uri, 0, $marker + 1);
	$part1 =~ s/[^\.\/]*$//;

	if ($part1 =~ /\.$/) {
	    chop $part1;
	}
	# else ($part1 =~ /\/$/)
	
	$new_uri = $part1;
	$new_uri .= $sid_path;
	$part1_len = length ($part1);
	$new_uri .= substr ($uri, $part1_len)  if ($part1_len < $uri_len);
    }

    return $new_uri;
}


# makeGetArgs (\%a)
#
#
#   desc:  make a string suitable for GET args out of associative array
#
#   args:  \%a    array with key->val pairs
#
#   rets:  $ret   string for GET args
#
#
sub makeGetArgs {
    local ($a) = @_;
    my ($k, $s, $ret);
    foreach $k (keys %$a) {
        $s = $a->{$k};
        $s = &escapeGetArg ($s);
        $ret .= "$k=$s&";
    }
    chop $ret;                                      # Remove trailing ampersand
    return $ret;
}
# end makeGetArgs ()


# escapeGetArg ()
#
sub escapeGetArg {
    local $str = shift;
    $str =~ s/ /+/go;
    $str =~ s/\%/\%25/go;
    $str =~ s/([^0-9a-zA-Z\%+])/"\%".&charToHex($1)/ge;
    return $str;
}
# end escapeGetArg ()


# unEscapeGetArg ()
#
sub unEscapeGetArg {
    local $str = shift;
    $str =~ s/\%25/PRESERVE_PERCENT_SYMBOL/go;
    $str =~ s/\%(\d[0-9a-fA-F])/&hexToChar($1)/ge;
    $str =~ s/PRESERVE_PERCENT_SYMBOL/\%/go;
    $str =~ s/\+/ /go;
    return $str;
}
# end escapeGetArg ()


# charToHex ()
#
sub charToHex {
    my $ascii = ord($_[0]);
    my %hexMap = (  0 => '0',
                    1 => '1',
                    2 => '2',
                    3 => '3',
                    4 => '4',
                    5 => '5',
                    6 => '6',
                    7 => '7',
                    8 => '8',
                    9 => '9',
                   10 => 'a',
                   11 => 'b',
                   12 => 'c',
                   13 => 'd',
                   14 => 'e',
                   15 => 'f'
                 );

    return $hexMap{(($ascii & 0xf0) >> 4)} . $hexMap{($ascii & 0x0f)};
}
# end charToHex ()


#  hexToChar ()
#
sub hexToChar {
    my $ascii = hex($_[0]);
    return chr $ascii;
}
# end hexToChar ()


1;

###############################################################################
# END
###############################################################################
