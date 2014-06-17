#!/usr/bin/perl

$out = $ARGV[0];
open(OUTFILE, ">$out");
$i = 0;
foreach (@ARGV){
    $i = $i +1;
	
	if ($i == 1){
		next;
	}
	
	$in = $_;
	open(INFILE, "$in");
	@vals = ();
	
	while(<INFILE>){
	    if (/#/){
	    #donothing
	    }else{
	        @data = split(/\s/, $_);
	        push(@vals, $data[1]);
	    }
	}
	
	#print "array - @vals\n";
	@folded = ();
	#print "len - $#vals\n";
	for ($c = 0; $c < $#vals/2; $c++){
		#print "$vals[$c]\t$vals[$#vals-$c]\n";
		push(@folded, $vals[$c]+$vals[$#vals-$c]);
	}
	push(@folded, $vals[$c]);
	
	foreach (@folded){
		print OUTFILE "$_ ";
	}
	print OUTFILE "\n";
}