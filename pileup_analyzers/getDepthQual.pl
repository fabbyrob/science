#!/usr/bin/perl
use POSIX;
print "INFILE DEPTH QUAL SNPQUAL AF1 [last 4 are filenames where data is output]\n";
$in = $ARGV[0];
$out = $ARGV[1];
$out2 = $ARGV[2];
$out3 = $ARGV[3];
$out4 = $ARGV[4];

open(INFILE, "$in");
open(OUTFILE, ">$out");
open(OUTFILE2, ">$out2");
open(OUTFILE3, ">$out3");
open(OUTFILE4, ">$out4");

$first = 1;

print ("FILES OPEN\n\n");
while(<INFILE>){
	if($first){
		$num = $#samps-8;
	}
    #/(1.chromosome)/t(2.site)/t(3.ReferenceBase)/t(4.InferedBase)/t(5.ConsensusScore)/t(6.SNPScore)/t(7.MapScore)/t(8.Depth)/t(9.Reads)/t(10.ReadScores)/
    if(/\S*\s*\S*\s*\S*\s*\S*\s*\S*\s*(\S*)\s*(\S*)\s*\S*DP=(\d*);AF1=(.*);\S*\s*\S*\s*(.*)/){
    	if($first){
            @samps = split(/\s/, $5);
            $num = ($#samps+1)*2;
    		$first = 0;
    		
    		for($i = 0; $i <= $num; $i++){
    		  $afs{$i} = 0;	
    		}
    	}
    	
    	
        if(not exists $depths{$3}){
            $depths{$3} = 1;
        } else{
            $depths{$3} = $depths{$3}+1;
        }
        
        $q = int($1);
        if(not exists $quals{$q}){
            $quals{$q} = 1;
        } else{
            $quals{$q} = $quals{$q}+1;
        }
        
        if(not exists $snpquals{$2}){
            $snpquals{$2} = 1;
        } else{
            $snpquals{$2} = $snpquals{$2}+1;
        }
        
        $af1 = floor($4*$num+0.5);

        if(not exists $afs{$af1}){
            $afs{$af1} = 1;
        } else{
            $afs{$af1} = $afs{$af1}+1;
        }
#        for($i = 1; $i < 26; $i++){
#        	if($4 < ($i+0.01)/26){
#        		if(not exists $afs{$i}){
#		            $afs{$i} = 1;
#		        } else{
#		            $afs{$i} = $afs{$i}+1;
#		        }
#		        last;
#        	}
#        }
        
    }
}

print OUTFILE "#DEPTH\tCOUNT\n";
foreach $key (sort { $a <=> $b } keys %depths){
    print OUTFILE "$key $depths{$key}\n";
}

print OUTFILE2 "#QUAL\tCOUNT\n";
foreach $key (sort { $a <=> $b } keys %quals){
    print OUTFILE2 "$key $quals{$key}\n";
}

print OUTFILE3 "#SNPQUAL\tCOUNT\n";
foreach $key (sort { $a <=> $b } keys %snpquals){
    print OUTFILE3 "$key $snpquals{$key}\n";
}

print OUTFILE4 "#AF1\tCOUNT\n";
foreach $key (sort { $a <=> $b } keys %afs){
    print OUTFILE4 "$key $afs{$key}\n";
}

print ("\nDONE");
close(INFILE);
close(OUTFILE);
close(OUTFILE2);
close(OUTFILE3);
close(OUTFILE4);
