#!/usr/bin/perl -w
use strict;
use warnings;
###
#set file input from command line argument
my $sample = $ARGV[0];
my $dir = $ARGV[1];
my $genesfile = $ARGV[2];
my $batch = $ARGV[3];

my $file = "$dir/star_2pass/ReadsPerGene.out.tab";

my $strandedoutput = "$dir/genecounts/$sample.$batch.raw.genecounts.txt";
my $rpkmoutput = "$dir/genecounts/$sample.genes.rpkm.gct";

open my $countsfile, '<', $file or die "Can't read $file : $!";
my @counts = <$countsfile>;
close $countsfile;
#my $genesfile = "~/data/target_coordinates/nmd274genes_gtf106_ensg.txt";
open my $geneslist, '<', $genesfile or die "Can't read $genesfile : $!";
my @genes = <$geneslist>;
close $geneslist;

my @strandedcounts;
my $strandedtotal = 0;

open my $strandedfile, '>', $strandedoutput or die "Can't read $strandedoutput : $!";

COUNTS: foreach my $countline (@counts){
	chomp $countline;
	my ($fullgene, $unstranded, $na, $stranded) = split('\t', $countline);
	my ($gene, $vers) = split('\.', $fullgene);
	GENE: foreach my $currentgene (@genes){
		chomp $currentgene;
		my ($geneid, $ensgid, $genelength) = split('\t', $currentgene);
		if ( $gene eq $geneid || $gene eq $ensgid ){
			my $strandedline = join("\t", $geneid, $stranded);
			print $strandedfile "$strandedline\n";
			my $strandedcalc = ((10 ** 9) * ($stranded)) / $genelength;
			my $rpkmline = join("\t", $ensgid, $geneid, $strandedcalc);
			push @strandedcounts, $rpkmline;
			$strandedtotal += $stranded;
			next COUNTS;
		} else {
			next GENE;
		}
	}
}
close $strandedfile;
open my $rpkmfile, '>', $rpkmoutput or die "Can't read $rpkmoutput : $!";
print $rpkmfile "##\n";
print $rpkmfile "###\n";
print $rpkmfile "Name\tDescription\t$sample\n";
foreach my $line (@strandedcounts){
	chomp $line;
	my ($ensg, $hgnc, $count) = split('\t', $line);
	my $rpkm = $count / $strandedtotal;
	print $rpkmfile "$ensg\t$hgnc\t$rpkm\n";
}
close $rpkmfile;