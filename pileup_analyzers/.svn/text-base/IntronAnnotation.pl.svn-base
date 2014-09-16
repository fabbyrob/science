#!/usr/bin/perl
$in = $ARGV[0];
$out = $ARGV[1];

open(INFILE, "$in");
open(OUTFILE, ">$out");

$previous = "";
$pend = 0;
while (<INFILE>){
	chomp;
	#print "$_";
	if (/(\S*)\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)/){
		#print "$previous\t$3\n";
		if ($previous eq $4 or $previous eq ""){#intron!
		    $s = $pend+1;
		    $e = $2-1;
		    $pend = $3;
		    if ($previous ne ""){
			 print OUTFILE "$1\t$s\t$e\t$4\t$5\t$6\n";
		    } else{
		    	$previous = $4;
		    }
		} else{#intergene!
			$previous = $4;
		}
	}
}

close(INFILE);
close(OUTFILE);