#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my $app = 'minirosetta';
my $boincrev = 15926;
my $output;

my $mini_url = "https://svn.rosettacommons.org/source/trunk/mini";
&GetOptions(
	"mini_url=s" => \$mini_url
);

my $version = $ARGV[0];
if ( !$version || $version !~ /(\d+)\.(\d+)/ ) {
	die "Error: don't recognize version ($version)!\n";
}

my $platform = lc `uname`; # should be one of linux, darwin or cygwin
$platform =~ s/[\s\n]+//g; # remove newlines and whitespace
if ( $platform =~ /cygwin/i ) {
	$platform = "cygwin";
}

print STDERR "building minirosetta on platform $platform.\n";

my $dir = "$app\_$version";
print $dir, "\n";

#if ( -d "./$dir" ) {
	#die "Error: directory $dir already exists! Build a new version or remove the old directory!\n";
#}

system( "mkdir -p $dir" );

# mini
if ( ! -d "$dir/mini" ) {
	my $checkout_cmd = "cd $dir; svn co $mini_url";
	print "checking out mini from $mini_url ";
	$output .= `$checkout_cmd`;
	print "done.\n";
}

# BOINC
if ( ! -d "$dir/mini/external/boinc" ) {
	my $boinc_checkout_cmd
		= "cd $dir/mini/external; svn co http://boinc.berkeley.edu/svn/trunk/boinc -r $boincrev";
	print "checking out BOINC libraries ... ";
	$output .= `$boinc_checkout_cmd`;
	print "done.\n";
}

if ( $platform ne 'cygwin' ) {
	my $boinc_build_cmd = <<BOINC_BUILD;
	cd $dir/mini/external/boinc
	./_autosetup
	./configure --disable-client --disable-server
	cd lib; make; cd ..
	cd api; make; cd ..
	cd zip; make; cd ..
BOINC_BUILD

	if ( $platform eq 'darwin' ) {
		system( "cp -r ~/jpeg6b $dir/mini/external/boinc" );
		$boinc_build_cmd .= "cd mac_build\n";
		$boinc_build_cmd .= "source BuildMacBOINC.sh -lib\n";
	}

	print "building boinc libraries ... ";
	$output .= `$boinc_build_cmd`;
	print "done.\n";
} else { # platform eq 'cygwin'
	my $cmd = "cd $dir/mini; python svn_version.py";
	system( $cmd );
}

# jpeg-6b
if ( $platform eq 'darwin' ) {
	my $jpeg_cmd = "cp -r $ENV{HOME}/jpeg-6b $dir/mini/external";
	$output .= `$jpeg_cmd`;
}

# glut
if ( $platform eq "cygwin" ) {
	my $glut_cmd = "cp -r ~/boinc_build/minirosetta_dependencies/glut $dir/mini/external";
	$output .= "$glut_cmd\n";
	$output .= `$glut_cmd`;
}

# Visual Studio setup
if ( $platform eq "cygwin" ) {
	my $vs_cmd = "cd \'$dir/mini/Visual Studio\'; source makeFiles.csh";
	print $vs_cmd, "\n";
	$output .= `$vs_cmd`;
}

# build applications
if ( $platform eq 'darwin' || $platform eq 'linux' ) {
	print "building minirosetta applications ... ";
	my $suffix;
	if ( $platform eq 'linux' ) {
		$suffix = "linuxgccrelease";
	} elsif ( $platform eq 'darwin' ){
		$suffix = "macosgccrelease";
	}

	my $build_cmd1 = "cd $dir/mini; scons -j4 mode=release extras=boinc,static bin/minirosetta.$suffix";
	my $build_cmd2 = "cd $dir/mini; scons -j4 mode=release extras=boinc,static bin/minirosetta_graphics.$suffix";

	# make a little file with the build.cmd in it.
	open FILE, ">$dir/mini/build.cmd" or die $!;
	print FILE $build_cmd1, "\n";;

	$output .= `$build_cmd1`;
	if ( $platform eq 'darwin' ) {
		$output .= `$build_cmd2`;
		print FILE $build_cmd2, "\n";;
	}
	close FILE or die $!;

}
print "done.\n";

# print the logfile
open LOG, ">$dir/autobuild.log" or die $!;
print LOG $output, "\n";
print LOG "finished building minirosetta on $platform.\n";
close LOG or die $!;

sub create_make_target_lines {
	my $target = shift;	
	my $lines  = shift;

}
