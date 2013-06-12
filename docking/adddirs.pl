#!/usr/bin/perl

if (@ARGV < 1) {
    
    while(<STDIN>) {
	
	print substr($_,0,2),"/",$_;
	
    };
    
}
else {
    
    open (IN,$ARGV[0]);
    while(<IN>) {
	
	print substr($_,0,2),"/",$_;
	
    };
    
};
