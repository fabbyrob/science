#!/usr/bin/perl

#for getting CDS locations from a gff3

$in = $ARGV[0];
$out = $ARGV[1];

open(INFILE, "$in");
open(OUTFILE1, ">$out");
open(ERROR, ">error.txt");
print ("FILES OPEN\n\n");
while(<INFILE>){
    chomp;
    if(/(\S*)\s*\S*\s*(\S*)\s*(\S*)\s*(\S*)\s*\S*\s*(\S*)\s*\S*\s*transcript_id \"(\S*)\";./){
        if ($2 eq "CDS"){
                print OUTFILE1 "$1\t$3\t$4\t$6\t$5\n";
        }else{
            print ERROR "Line was not gene or CDS:\n\t$_\n";
        }
    }elsif((/(\S*)\s*\S*\s*(\S*)\s*(\S*)\s*(\S*)\s*\S*\s*(\S*)\s*\S*\s*Parent=([^.]*)./)){
        if ($2 eq "CDS"){
                print OUTFILE1 "$1\t$3\t$4\t$6\t$5\n";
        }else{
            print ERROR "Line was not gene or CDS:\n\t$_\n";
        }
    }else{
        print ERROR "Line did not match filter:\n\t$_\n";
    }
}
print ("DONE\n");
close(INFILE);
close(OUTFILE);
