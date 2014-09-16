#!/usr/bin/perl

#Forgetting 5UTR and 3UTR locations fom a gff3

$in = $ARGV[0];
$out = $ARGV[1];
$out2 = $ARGV[2];

open(INFILE, "$in");
open(OUTFILE1, ">$out");
open(OUTFILE2, ">$out2");
open(ERROR, ">error.txt");
print ("FILES OPEN\n\n");
while(<INFILE>){
    chomp;
    if(/(\S*)\s*\S*\s*(\S*)\s*(\S*)\s*(\S*)\s*\S*\s*(\S*)\s*\S*\s*ID=([^.]*)./){
        if ($2 eq "five_prime_UTR" or $2 eq "three_prime_UTR"){
            if ($2 eq "five_prime_UTR"){
                print OUTFILE1 "$1\t$3\t$4\t$6\t$5\n";
            } else{
                print OUTFILE2 "$1\t$3\t$4\t$6\t$5\n";
            }
        }else{
            print ERROR "Line was not gene or CDS:\n\t$_\n";
        }
    }elsif((/(\S*)\s*\S*\s*(\S*)\s*(\S*)\s*(\S*)\s*\S*\s*(\S*)\s*\S*\s*Parent=([^.]*)./)){
        if ($2 eq "five_prime_UTR" or $2 eq "three_prime_UTR"){
            if ($2 eq "five_prime_UTR"){
                print OUTFILE1 "$1\t$3\t$4\t$6\t$5\n";
            } else{
                print OUTFILE2 "$1\t$3\t$4\t$6\t$5\n";
            }
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
