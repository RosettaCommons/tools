#!/usr/bin/perl -w

use strict;
use File::Find;

my $base_path   = '/Users/tex/src/mini';
my @source_dirs = qw| src demo test |;
my @extensions  = qw/ cc hh /;

for my $source_dir (@source_dirs) {
  find( 
    { wanted => \&spaces2tabs },
    "$base_path/$source_dir"
  );
}

# Takes a file as an argument and removes all trailing whitespace from its
# lines.

sub spaces2tabs {
  my $file = $_;

  if ( $file =~ /\.svn/ ) {
    #print "ignoring $file\n";
    return 0;
  }

  if ( ! -f $file ) {
    #print "$file isn't a file!\n";
    return 0;
  }

  my $success = 0;
  for my $extension ( @extensions ) {
    if ( $file =~ /\.$extension$/ ) {
      $success = 1;
      last;
    }
  }
  if ( ! $success ) {
    return 0;
  }

  # warning: changing the following code will void your warranty!
  print "converting spaces to tabs in $file ... ";
  open FILE, "<$file" or die $!;
  my @file = <FILE>;
  close FILE or die $!;

  open FILE, ">$file" or die $!;
  foreach my $line (@file) {
    chomp $line;
		#$line =~ s/^\t/  /g;
		if ( $line =~ /^(\s+)(.*)$/ ) {
			my $space = $1;
			my $rest  = $2;
			$space =~ s/  /\t/g;
			$line = $space . $rest;
		}
		#$line =~ s/^  /\t/g;
    print FILE $line, "\n";
  }
  close FILE or die $!;
  print "done.\n";
  # end warning.
}
