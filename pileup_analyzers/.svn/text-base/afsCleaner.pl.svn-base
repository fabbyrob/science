#!/usr/bin/perl
$in = $ARGV[0];
$out = $ARGV[1];
$n = $ARGV[2];

open(INFILE, "$in");
open(OUTFILE, ">$out");

while(<INFILE>){
	if (/\[afs\]/){
		@data = split(/\s/, $_);
		foreach (@data){
			@items = split(/:/, $_);
			if ($items[0] > $n){
				last;
			}
			print OUTFILE "$items[1] "
		}
		print OUTFILE "\n"
	}
}