#!/usr/bin/perl
#take an annotated reference, and prints out the list of sites that fit into a given category
print "INFILE OUTFILE [NUMBERS YOU WANT]\n";
$in = $ARGV[0];
$out = $ARGV[1];
open(INFILE, "$in");
open(OUTFILE, ">$out");

$category = $ARGV[2];

print "FILES OPEN\n";

while(<INFILE>){
	chomp;
	if(/(\S*)\s*(\d*)\s*\w*\s*\w*\s*\S*\s*(\d*)/){
		if ($3 == $category){
			print OUTFILE "$1\t$2\n";
		}
	}
}
close(INFILE);
close(OUTFILE);
print "DONE\n";
