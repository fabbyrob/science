#!/usr/bin/perl
#take an annotated reference, and prints out the list of sites that fit into a given category
print "INFILE OUTFILE \n";
$in = $ARGV[0];
$out = $ARGV[1];
open(INFILE, "$in");
open(OUTFILE, ">$out");

print "FILES OPEN\n";

%other = ();
%samps1 = ();
%samps2 = ();

$samp = "";
$n = 0;
while(<INFILE>){
	chomp;
    if(/>\s?(\S*)(.*)/){
        $samp = $1;
        $other{$samp} = $2;
        #print "$2\n";
        if (exists $samps1{$samp}){
        	$n = 1;
        	$samps2{$samp} = ""
        }
        else{
        	$n = 0;
            $samps1{$samp} = ""
        }
    } else{
    	if ($n){#already have one sample by this name
            $tmp = $samps2{$samp};
            $samps2{$samp} = "$tmp$_";
    	}
    	else{
    		$tmp = $samps1{$samp};
    		$samps1{$samp} = "$tmp$_";
    	}
    }
}

while (($key, $value) = each(%samps1)){
	#print ("$key - \n\t$value\n");
	if (exists $samps2{$key}){#have both strands
		@one = split(//, $value);
		@two = split(//, $samps2{$key});
		@res = ();
		for ($i = 0; $i <= $#one; $i++){
			if (@one[$i] eq @two[$i]){
				@res[$i] = @one[$i];
			} else{
				@res[$i] = "N";
			}
		}
		$samps1{$key} = join('', @res);
	}
}

$others = $other{'Rubella'};
print OUTFILE ">Rubella $others\n";
$value = $samps1{'Rubella'};
print OUTFILE "$value\n\n";

while (($key, $value) = each(%samps1)){
	if ($key ne 'Rubella'){
        print OUTFILE ">$key$other{$key}\n";
        print OUTFILE "$value\n\n";
	}
}

close(INFILE);
close(OUTFILE);
print "DONE\n";
