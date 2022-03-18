use strict;
use warnings;
use List::Util qw(sum);

my%dropgenes = ("gene27686","1",
"_ambiguous","1",
"_unmapped","1",
"_no_feature","1",
"gene37698","1",
"gene37699","1"
);

open(IN,"../readcounts_EXP0.txt");
my@samples;
my%data;
my%genecount;
while(<IN>){
	chomp;
	if($. == 1){
		@samples = split("\t");
	} else {
		my@l = split("\t");
		
		if(exists($dropgenes{$l[0]})){next};
		
		for my$i (1 .. $#l){
			if($samples[$i] =~ m/^(H\d+)_/){
				$data{$1}{$l[0]} += $l[$i];
				$genecount{$l[0]} += $l[$i];
			}
		}
	}
}
close IN;

foreach my$gene (keys %genecount){
	if($genecount{$gene} < 100){
		foreach my$volunteer (keys %data){
			delete $data{$volunteer}{$gene};
			delete $genecount{$gene};
		}
	}
}
my@sortedgenes = sort {$a cmp $b} keys %genecount;

my%medcorr;
foreach my$volunteer (keys %data){
	my@sortedarray = sort {$a <=> $b} values(%{$data{$volunteer}});
	my$median = $sortedarray[int($#sortedarray/2)];
	print "Median ".$volunteer."\t".$median."\t".sum(@sortedarray)."\n";
	$medcorr{$volunteer} = $median;
}
my@volsort = sort {$a cmp $b} keys %data;

open(RESP,'../response.csv');
my%resp;
while(<RESP>){
	if(m/"(H\d+)",(\d)/){
		$resp{$1} = $2;
	}
}
close RESP;

open(OUT,">day0_normcounts.txt");
print OUT "Patient\t".join("\t",@sortedgenes)."\tResponse\n";
foreach my$vol (@volsort){
	if(exists($resp{$vol})){
	print OUT $vol;
	foreach my$gene (@sortedgenes){
		print OUT "\t".int($data{$vol}{$gene}*10000/$medcorr{$vol})/10000
	}
	print OUT "\t".$resp{$vol}."\n";
	}
}
close OUT;