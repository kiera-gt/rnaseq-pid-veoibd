#!/usr/bin/perl -w
use strict;
use warnings;
my $basedir = $ARGV[0];
my $type = $ARGV[1];

my $sampleroster = "$basedir/$type.splicefiles.txt";
open my $samplelist, '<', $sampleroster or die "Can't read $sampleroster : $!";
my @samples = <$samplelist>;
close $samplelist;

my $intronfile = "$basedir/$type.sjs.uniq.txt";
open my $introns_from_samples, '<', $intronfile or die "Can't read $intronfile : $!";
my @introns = <$introns_from_samples>;
close $introns_from_samples;

my $outputfile = "$basedir/$type.splices.txt";
my $outputfile2 = "$basedir/$type.splices.sorted.txt";
###
#separate annotated and unannotated
my @annotated;
my @unannotated;
foreach my $intron (@introns) {
	chomp $intron;
	if ($intron =~ /UNANNOTATED/){
		push @unannotated, $intron;
	} else {
		push @annotated, $intron;
	}
}
###
#put back in the correct order
my $innerindex2 = 0;
my $outer2;
my $inner2;
my $nmd2;
if ($#annotated >= $#unannotated){
	$outer2 = \@annotated;
	$inner2 = \@unannotated;
	$nmd2 = 1;
} else {
	$outer2 = \@unannotated;
	$inner2 = \@annotated;
	$nmd2 = 2;
}
my @finalgtf;
OUTER: foreach my $entry1 ( @{$outer2} ){
	chomp $entry1;
	my ($gene1, $transcript1, $strand1, $chr1, $start1, $stop1) = split('\t', $entry1);
	INNER:
	{
		if (! defined ${$inner2}[$innerindex2]) {
			my $gtfline = join( "\t", $gene1, $transcript1, $strand1, $chr1, $start1, $stop1 );
			push @finalgtf, $gtfline;
			next OUTER;
		}
		chomp ${$inner2}[$innerindex2];
		my ($gene2, $transcript2, $strand2, $chr2, $start2, $stop2) = split('\t', ${$inner2}[$innerindex2]);
		if ($chr1 !~ /^$chr2$/){
			if ($chr1 eq 'Y' && $chr2 eq 'X'){
				my $gtfline = join( "\t", $gene2, $transcript2, $strand2, $chr2, $start2, $stop2 );
				push @finalgtf, $gtfline;
				$innerindex2++;
				redo INNER;
			} elsif ($chr1 eq 'X' && $chr2 eq 'Y'){
				my $gtfline = join( "\t", $gene1, $transcript1, $strand1, $chr1, $start1, $stop1 );
				push @finalgtf, $gtfline;
				next OUTER;
			} elsif ($chr1 == 22 && $chr2 eq 'X'){
				my $gtfline = join( "\t", $gene1, $transcript1, $strand1, $chr1, $start1, $stop1 );
				push @finalgtf, $gtfline;
				next OUTER;
			} elsif ($chr1 eq 'X' && $chr2 == 22){
				my $gtfline = join( "\t", $gene2, $transcript2, $strand2, $chr2, $start2, $stop2 );
				push @finalgtf, $gtfline;
				$innerindex2++;
				redo INNER;
			} elsif ($chr1 > $chr2){
				my $gtfline = join( "\t", $gene2, $transcript2, $strand2, $chr2, $start2, $stop2 );
				push @finalgtf, $gtfline;
				$innerindex2++;
				redo INNER;
			} else {
				my $gtfline = join( "\t", $gene1, $transcript1, $strand1, $chr1, $start1, $stop1 );
				push @finalgtf, $gtfline;
				next OUTER;
			}
		} else {
			if ($start1 != $start2){
				if ($start2 > $start1){
					#inner is further along than outer
					my $gtfline = join( "\t", $gene1, $transcript1, $strand1, $chr1, $start1, $stop1 );
					push @finalgtf, $gtfline;
					next OUTER;
				} else {
					#outer is further along than inner
					my $gtfline = join( "\t", $gene2, $transcript2, $strand2, $chr2, $start2, $stop2 );
					push @finalgtf, $gtfline;
					$innerindex2++;
					redo INNER;
				}
			} else {
				if ($stop1 != $stop2){
					if ($stop2 > $stop1){
						#inner is further along than outer
						my $gtfline = join( "\t", $gene1, $transcript1, $strand1, $chr1, $start1, $stop1 );
						push @finalgtf, $gtfline;
						next OUTER;
					} else {
						#outer is further along than inner
						my $gtfline = join( "\t", $gene2, $transcript2, $strand2, $chr2, $start2, $stop2 );
						push @finalgtf, $gtfline;
						$innerindex2++;
						redo INNER;
					}
				} elsif ($strand1 ne $strand2) {
					my $gtfline = join( "\t", $gene1, $transcript1, $strand1, $chr1, $start1, $stop1 );
					my $intronline = join( "\t", $gene2, $transcript2, $strand2, $chr2, $start2, $stop2 );
					if ($nmd2 == 1){
						push @finalgtf, $gtfline;
						push @finalgtf, $intronline;
					} else {
						push @finalgtf, $intronline;
						push @finalgtf, $gtfline;
					}
					$innerindex2++;
					next OUTER;
				} else {
					my $gtfline = join( "\t", $gene1, $transcript1, $strand1, $chr1, $start1, $stop1 );
					my $intronline = join( "\t", $gene2, $transcript2, $strand2, $chr2, $start2, $stop2 );
					if ($nmd2 == 1){
						push @finalgtf, $gtfline;
						push @finalgtf, $intronline;
					} else {
						push @finalgtf, $intronline;
						push @finalgtf, $gtfline;
					}
					$innerindex2++;
					next OUTER;
				}
			}
		}
	}
}
#This section goes through each sample's SJ output file and finds a count for each gtf line
my $sampleindex = 0;
my $finalgtfindex = 0;
my %sjcounts;
my $columninfo = join(":", "gene_id", "transcripts", "strand", "chr", "start", "end:");
my @outputheader;
push @outputheader, $columninfo;
foreach my $sample (@samples) {
	chomp $sample;
	my ($sampleid, $sampleheader, $path) = split('\t', $sample);
	my $sample_file = "$path/splicing/$sampleid.nodups.raw.splices.txt";
	open my $sjfile, '<', $sample_file or die "Can't read $sample_file : $!";
	my @samplesjs = <$sjfile>;
	close $sjfile;
	push @outputheader, $sampleheader;
	#$sjcounts{'sample_id'}[$sampleindex] = $sampleid;
	my $sample_sjindex = 0;
	FINALGTF: foreach my $finalgtfline (@finalgtf) {
		chomp $finalgtfline;
		my ($finalgene, $finaltranscript, $finalstrand, $finalchr, $finalstart, $finalstop) = split('\t', $finalgtfline);
		SAMPLEFILE:
		{
			if (! defined $samplesjs[$sample_sjindex]) {
				my $key = join(":", $finalgene, $finaltranscript, $finalstrand, "chr$finalchr", $finalstart, $finalstop);
				$sjcounts{$key}[$sampleindex] = 0;
				$sample_sjindex++;
				next FINALGTF;
			}
			chomp $samplesjs[$sample_sjindex];
			my ($samplegene, $sampletrans, $samplestrand, $samplechr, $samplestart, $samplestop, $samplecount) = split('\t', $samplesjs[$sample_sjindex]); #split each line by the tabs
			$samplechr =~ s/chr//;
			if ($samplechr ne $finalchr){
				if ($samplechr eq 'M'){
					$sample_sjindex++;
					redo SAMPLEFILE;
				} elsif ($finalchr eq 'Y' and $samplechr eq 'X'){
					$sample_sjindex++;
					redo SAMPLEFILE;
				} elsif ($finalchr eq 'X' and $samplechr eq 'Y'){
					my $key = join(":", $finalgene, $finaltranscript, $finalstrand, "chr$finalchr", $finalstart, $finalstop);
					$sjcounts{$key}[$sampleindex] = 0;
					next FINALGTF; #go to next gtf entry but keep same SJ entry
				} elsif ($finalchr == 22 and $samplechr eq 'X'){
					my $key = join(":", $finalgene, $finaltranscript, $finalstrand, "chr$finalchr", $finalstart, $finalstop);
					$sjcounts{$key}[$sampleindex] = 0;
					next FINALGTF; #go to next gtf entry but keep same SJ entry
				} elsif ($finalchr eq 'X' and $samplechr == 22){
					$sample_sjindex++;
					redo SAMPLEFILE;
				} elsif ($finalchr > $samplechr){
					$sample_sjindex++;
					redo SAMPLEFILE;
				} else {
					my $key = join(":", $finalgene, $finaltranscript, $finalstrand, "chr$finalchr", $finalstart, $finalstop);
					$sjcounts{$key}[$sampleindex] = 0;
					next FINALGTF; #go to next gtf entry but keep same SJ entry
				}
			} elsif ($samplestart != $finalstart){
				if ($samplestart < $finalstart){
					$sample_sjindex++;
					redo SAMPLEFILE;
				} else {
					my $key = join(":", $finalgene, $finaltranscript, $finalstrand, "chr$finalchr", $finalstart, $finalstop);
					$sjcounts{$key}[$sampleindex] = 0;
					next FINALGTF; #go to next gtf entry but keep same SJ entry
				}
			} elsif ($samplestop != $finalstop){
				if ($samplestop < $finalstop){
					$sample_sjindex++;
					redo SAMPLEFILE;
				} else {
					my $key = join(":", $finalgene, $finaltranscript, $finalstrand, "chr$finalchr", $finalstart, $finalstop);
					$sjcounts{$key}[$sampleindex] = 0;
					next FINALGTF; #go to next gtf entry but keep same SJ entry
				}
			} elsif ($samplestrand ne $finalstrand){
				$sample_sjindex++;
				redo SAMPLEFILE;
			} else {
				my $key = join(":", $finalgene, $finaltranscript, $finalstrand, "chr$finalchr", $finalstart, $finalstop);
				$sjcounts{$key}[$sampleindex] = $samplecount;
				$sample_sjindex++;
				next FINALGTF;
			}
		}
	}
	$sampleindex++;
}

open my $sortedoutput, '>', $outputfile2 or die "Can't read $outputfile2 : $!";
print $sortedoutput join("\t",@outputheader),"\n";
close $sortedoutput;
open my $output, '>', $outputfile or die "Can't read $outputfile : $!";
foreach my $entry (keys %sjcounts) { 
	my $keycounts = join("\t", @{$sjcounts{$entry}});
	print $output "$entry\t$keycounts\n";
}
close $output;