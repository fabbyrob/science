#!/usr/bin/perl
$in = $ARGV[0];
$out = $ARGV[1];
@codes{@ARGV[2..$#ARGV]}=();

open(INFILE, "$in");
open(OUTFILE, ">$out");

print ("FILES OPEN\n\n");
while(<INFILE>){
	chomp;
    if(/(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)/){
    	if (exists $codes{$6}){
    		#print OUTFILE "$1\t$2\n";
    		print OUTFILE "$_\n";
    	}
    }
}

close(INFILE);
close(OUTFILE);
