#!/usr/bin/env perl
## specially for convert into CNVnator genotype input format

use File::Basename;
die("Argument: .genotype.vcf[.gz]\n") if (@ARGV != 1);
my $vcf = shift @ARGV;
my ($fn, $dir, undef) = fileparse($vcf, qw/.vcf.gz .vcf/);
my ($chr_col, $pos_col, $id_col, $ref_col, $alt_col, $info_col, $format_col, $sample_col) = (1-1, 2-1, 3-1, 4-1, 5-1, 8-1, 9-1, 10-1);

if($vcf =~ /\.vcf\.gz$/){
	open VCF, "zcat -c $vcf|" or die("Cannot open $vcf!\n");
}elsif($vcf =~ /\.vcf$/){
	open VCF, "<$vcf" or die("Cannot open $vcf!\n");
}else{
	die("Please give .vcf(.gz) as first argument!\n");
}

my $header;
while($line=<VCF>){
	next if ($line=~/^##/);
	$header = $line;
	last;
}
my @header = split /\s/, $header;

while($line=<VCF>){
	my $cn_num_flag = 0;
	my @line = split /\s+/, $line;
        my $chr = $line[$chr_col];
	my $pos = $line[$pos_col];
	my $id = $line[$id_col];
	my $ref = $line[$ref_col];
	if ($ref =~ /CN(\d+)/){
		$cn_num_flag = 1;
		$ref = $1;
	}
	my @alt = split /,/, $line[$alt_col];
	for $i (0..$#alt){
		if($alt[$i] =~ /CN(\d+)/){
			$cn_num_flag = 1;
			$alt[$i] = $1;
		}else{
			die("Check $line") if ($cn_num_flag);
		}
	}
	my @alleles = ($ref,@alt);
	my @info = split /;/, $line[$info_col];
	my ($svtype, $end, $svlen);
	for $item (@info){
		if($item =~ /^SVTYPE=/){
			$svtype = $';
		}elsif($item =~ /^END=/){
			$end = $';
		}elsif($item =~ /^SVLEN=/){
			$svlen = $';
		}else{
		}
	}
	die("1.Check $line") if (length($svtype)==0);
        if(length($end)*length($svlen)!=0){
		#               print "2. Check $line" if (($pos - $end + 1)!=$svlen)
	}
	my @formats = split /:/, $line[$format_col];
	my $GT_flag;
	my $genotype_index;
	for $i (0..$#formats){
		if($formats[$i] eq 'GT'){			# for vcf genotype format ./. or .|.
			$GT_flag = 1;
			$genotype_index = $i;
		}
		if($formats[$i] eq 'CN'){		# for single copynumber format .
			$GT_flag = 0;
			$genotype_index = $i;
			last;				# 'CN' has higher priority
		}
	}
	die("Unknown genotype format in $vcf!\n") if (length($GT_flag)<=0);
	if(length($end)!=0){
#		push @output, "$id\t$chr\t$pos\t$end\t$ref\t$line[$alt_col]\t$svtype";
		$chr =~ s/^23$/X/;
		$chr =~ s/^24$/Y/;
		push @output, "${chr}:${pos}-${end}\n";
	}elsif(length($svlen)!=0){
		my $sv_end = $pos + $svlen - 1;
#		push @output, "$id\t$chr\t$pos\t$sv_end\t$ref\t$line[$alt_col]\t$svtype";
		$chr =~ s/^23$/X/;
		$chr =~ s/^24$/Y/;
		push @output, "${chr}:${pos}-${sv_end}\n";
	}else{
		die("Check $line!");
	}
#	for $i ($sample_col..$#line){
#		my $gt = (split /:/, $line[$i])[$genotype_index];
#		if ($gt =~ /\./){
#		## Missing data
#			$output[-1] .= "\tNA";
#			next;
#		}
#		if($GT_flag){
#			## SNP-like genotype format: ./. or .|.
#			## convert into copy number
#			my @gt = (split /\/|\|/, $gt);
#			my $cn = 0;
#			foreach $a (@gt){
#				if($line[$alt_col] =~ /del/i){
#				## DEL, GenomeSTRiP deletion vcf
#					if($a == 0){
#						$cn += 1;
#					}elsif($a == 1){
#						$cn += 0;
#					}else{
#						die("Unknown genotype in $line");
#					}
#				}elsif($cn_num_flag){
#				## Copy Numbers, like <CN#>
#					if($a == 0){
#						$cn += $ref;
#					}elsif($a == 1){
#						$cn += $alt[0];
#					}else{
#						die("Unknown genotype in $line");
#					}
#				}else{
#					die("Check $line in $vcf!\n");
#				}
#			}
#			$output[-1] .= "\t$cn";
#		}else{
#			## CNV, GenomeSTRiP CNV vcf
#			$output[-1] .= "\t$gt";
#		}
#	}
#	$output[-1] .= "\n";
}
close VCF;

$calls_file = $dir . $fn . '.region.txt';
open CALL, ">$calls_file" or die("Cannot open $calls_file!\n");
#print CALL (join "\t",('ID','Chr','Pos','End','Ref','Alt','SVtype',@header[$sample_col..$#header]));
#print CALL "\n";
print CALL @output;
close CALL;
print "Finish $vcf!\n";
