#!/usr/bin/perl
use POSIX;
#print "INFILE outfile sampSize\n";
$in = $ARGV[0];
$out = $ARGV[1];
$num = $ARGV[2];

open(INFILE, "$in");
open(OUTFILE, ">$out");

#print ("FILES OPEN\n\n");
while(<INFILE>){
    #/(1.chromosome)/t(2.site)/t(3.ReferenceBase)/t(4.InferedBase)/t(5.ConsensusScore)/t(6.SNPScore)/t(7.MapScore)/t(8.Depth)/t(9.Reads)/t(10.ReadScores)/
    if(/^(\S+)\s+(\d+).+AF=(\S+);/){
        $af1 = floor($3*$num+0.5);

        print OUTFILE "$1\t$2\t$af1\n";
        
    }
}

#print ("\nDONE");
close(INFILE);
close(OUTFILE);
