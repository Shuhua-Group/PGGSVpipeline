#!/usr/bin/env perl
## specially for convert into CNVnator genotype input format

use File::Basename;
die("Argument: .genotype.vcf[.gz] region.Regenotype\n") if (@ARGV != 2);
my $vcf = shift @ARGV;
my $geno = shift @ARGV;
my ($fn, $dir, undef) = fileparse($vcf, qw/.vcf.gz .vcf/);
my ($chr_col, $pos_col, $id_col, $ref_col, $alt_col, $info_col, $format_col, $sample_col) = (1-1, 2-1, 3-1, 4-1, 5-1, 8-1, 9-1, 10-1);
my $cnvnator_genotype_col = -2;			# genotype coloumn of cnvnator genotype file (last but one)
my $cnvnator_region_col = 2-1;
my $cutoff = 2;

my $region2geno_exclude={};
open GENO,"<$geno" or die("Cannot open $geno!\n");
while($line=<GENO>){
	next if $line =~ /^Assuming/;
	my @line = split /\s+/, $line;
	my ($region, $genotype) = ($line[$cnvnator_region_col], $line[$cnvnator_genotype_col]);
	$region =~ s/^chr//i; $region =~ s/^X:/23:/i; $region =~ s/^Y:/24:/i;
	if (! exists $region2geno_exclude{$region}){
		$region2geno_exclude{$region}=$genotype;
	}else{
		die("Duplicated region in $geno!\n");
	}
}
close GENO;

if($vcf =~ /\.vcf\.gz$/){
	open VCF, "zcat -c $vcf|" or die("Cannot open $vcf!\n");
}elsif($vcf =~ /\.vcf$/){
	open VCF, "<$vcf" or die("Cannot open $vcf!\n");
}else{
	die("Please give .vcf(.gz) as first argument!\n");
}

my @output;

while($line=<VCF>){
	if ($line =~ /^#/){
		push @output, $line;
		next;
	}
	my $cn_num_flag = 0;
	my @line = split /\s+/, $line;
        my $chr = $line[$chr_col];
	$chr =~ s/^chr//i; $chr =~ s/^X:/23:/i; $chr =~ s/^Y:/24:/i;
	my $pos = $line[$pos_col];
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
	if(length($end)!=0){
		if(exists $region2geno_exclude{"${chr}:${pos}-${end}"}){
			next if($svtype=~/^del$/i and $region2geno_exclude{"${chr}:${pos}-${end}"} > $cutoff);
			next if($svtype=~/^dup$/i and $region2geno_exclude{"${chr}:${pos}-${end}"} < $cutoff);
		}
	}elsif(length($svlen)!=0){
		my $sv_end = $pos + $svlen - 1;
		if(exists $region2geno_exclude{"${chr}:${pos}-${sv_end}"}){
			next if($svtype=~/^del$/i and $region2geno_exclude{"${chr}:${pos}-${sv_end}"} > $cutoff);
			next if($svtype=~/^dup$/i and $region2geno_exclude{"${chr}:${pos}-${sv_end}"} < $cutoff);
		}
	}else{
		die("Check $line!");
	}
	push @output, $line;
}
close VCF;

$pass_file = $dir . $fn . '.reGenoPass.vcf';
open PASS, ">$pass_file" or die("Cannot open $pass_file!\n");
print PASS @output;
close PASS;
print "Finish $vcf!\n";
