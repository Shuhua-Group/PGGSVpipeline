#!/usr/bin/env perl
use File::Basename;
use Env;
die("Argument: indi_Dir SampleID\n") if (@ARGV !=2);
my $indi_dir = shift @ARGV;
my $sample = shift @ARGV;
if(! -e $indi_dir){
	die("Please give a valid directory of individual MetaSV results!\n");
}

my @annotations;
my $header;
my @SVs;

for $c (1..24,'X','Y'){
	next if ((! -e "${indi_dir}/chr${c}") and (! -e "${indi_dir}/Chr${c}"));
	my $cur_dir;
	if(-e "${indi_dir}/chr${c}"){
		$cur_dir = "${indi_dir}/chr${c}";
	}elsif( -e "${indi_dir}/Chr${c}"){
		$cur_dir = "${indi_dir}/Chr${c}";
	}else{
		die("Impossible!\n");
	}
	my $vcf = "${cur_dir}/variants.vcf.gz";
	open VCF, "zcat -c $vcf |" or die("Cannot open $vcf!\n");
	if(!@annotations){						# the basic annotation lines (,as the first chromosome/vcf)
		while($line=<VCF>){
			if($line =~ /^##/){
				push @annotations, $line;
			}else{
				$header = $line;
				@SVs = <VCF>;
			}
		}
	}else{
		while($line=<VCF>){					# annotation lines of rest files
			if($line =~ /^##/){
				my $item = (split /=/,$line)[0];
				my $pos = -1;
				for $i (0..$#annotations){
					my $anno_item = (split /=/, $annotations[$i])[0];	# I try to put the annotation lines saying the same information together.
					next if ($item ne $anno_item);
					$pos = $i;
					if($line eq $annotations[$i]){	# already existed line
						$pos = -2;
						last;
					}
				}
				if($pos == -2){
					## ready existed annotation line
				}elsif(($pos == -1) or ($pos == $#annotation)){
					## not existed, add to the last
					push @annotations, $line;
				}else{
					## with same information, but different discreption, put together
					@annotations = (@annotations[0..$pos],$line,@annotations[($pos+1)..$#annotations]);
				}
			}else{
				die("Unmatched header for $vcf!\n") if ($header ne $line);
				@SVs = (@SVs, <VCF>);
			}
		}
	}
	close VCF;
}

my @genome = (@annotations, $header, @SVs);
my $genome = "${indi_dir}/${sample}.MetaSV.vcf";
open VCF,">$genome" or die("Cannot open $genome!\n");
print VCF @genome;
close VCF;
print "Finish $sample!\n";
