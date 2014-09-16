#!/usr/bin/perl
$in = $ARGV[0];
$out = $ARGV[1];

open(INFILE, "$in");
open(OUTFILE, ">$out");

while(<INFILE>){
	if(/\D*(\d*)\s*(\d*)\s*(\d*)(.*)/){
		$i = $2;
		for($i = $2; $i <= $3; $i++){
			print OUTFILE "scaffold_$1\t$i$4\n";
		}
	}
}

print ("\nDONE");
close(INFILE);
close(OUTFILE);
