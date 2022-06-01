#!/usr/bin/perl -w
use strict;
use warnings;
#
my $path = $ARGV[0];
chomp $path;
my $sample = $ARGV[1];
chomp $sample;
my $geneintervalsfile = $ARGV[2];
chomp $geneintervalsfile;
my $introngtf = $ARGV[3];
chomp $introngtf;
my $batch = $ARGV[4];
chomp $batch;
###
my $outputfile = "$path/splicing/$sample.nodups.splices.txt";
my $outputfile2 = "$path/splicing/$sample.nodups.raw.splices.txt";
###
#This section takes the list of SJs among all samples and uses the gene intervals file to assign gene name to SJ entries (if they are within one of the 274 genes of interest)
#my $geneintervalsfile = '/nv/hp10/kberger9/data/target_coordinates/nmd274genes_intervals_vsort.txt';
open my $geneintervals, '<', $geneintervalsfile or die "Can't read $geneintervalsfile : $!";
my @genes = <$geneintervals>;
close $geneintervals;
###
my $intronfile = "$path/splicing/$sample.$batch.nodups.sj.out.tab";
open my $introns_from_sample, '<', $intronfile or die "Can't read $intronfile : $!";
my @introns = <$introns_from_sample>;
close $introns_from_sample;
###
#my $introngtf = '/nv/hp10/kberger9/data/target_coordinates/condensed_introns_nooverlap.txt';
open my $gtf, '<', $introngtf or die "Can't read $introngtf : $!";
my @gtfintrons = <$gtf>;
close $gtf;
###
#this section finds all uniq SJs from sample within the panel genes
my @labeledintrons;
my $sjindex = 0;
###
GENES: foreach my $currentgene (@genes) {
	chomp $currentgene;
	my ($geneid, $genestrand, $genechr, $genestart, $genestop) = split('\t', $currentgene);
	$genechr =~ s/chr//;
	INTRONFILE:
	{
		if (! defined $introns[$sjindex]){
			next GENES;
		}
		#iterate over file but use array slice to skip over entries no longer useful
		chomp $introns[$sjindex];
		my @splitsj = split('\t', $introns[$sjindex]); #split each line by the tabs
		my $sjchr = $splitsj[0]; #add column variable
		my $sjstart = $splitsj[1];
		my $sjstop = $splitsj[2];
		my $samplecount = $splitsj[6];
		$sjchr =~ s/chr//;
		my $sjstrand;
		if ($splitsj[3] == 0) {
			$sjstrand = 0;
		} elsif ($splitsj[3] == 1) {
			$sjstrand = '+';
		} elsif ($splitsj[3] == 2) {
			$sjstrand = '-';
		}
		#chomp $introns[$sjindex];
		#my ($sjchr, $sjstart, $sjstop, $sjstrand) = split('\t', $introns[$sjindex]);
		#$sjchr =~ s/chr//;
		my $transcript = "UNANNOTATED";
		my $type = "intron";
		my $sjgene = "UNKNOWN";
		if ($genechr !~ /^$sjchr$/) {
			if ($sjchr !~ /^([1-9]|1[0-9]|2[0-2]|[XY])$/) {
				#sjchr is not a number and is not X or Y, cannot be relevant to our genes
				$sjindex++; #increase sjindex by one
				redo INTRONFILE; #go to next intron line but keep same gene interval
			} elsif ($genechr =~ /^([1-9]|1[0-9]|2[0-2])$/ && $sjchr =~ /^[X]$/){
				#genechr is 1-22 and sjchr is X
				next GENES;
			} elsif ($genechr =~ /^[X]$/){
				if ($sjchr =~ /^([1-9]|1[0-9]|2[0-2])$/){
					#genechr is further along than sjchr
					$sjindex++; #increase sjindex by one
					redo INTRONFILE; #go to next intron line but keep same gene interval
				} else {
					#sjchr is Y and genechr is X
					next GENES;
				}
			} elsif ($genechr =~ /^[Y]$/ && $sjchr =~ /^([1-9]|1[0-9]|2[0-2]|[X])$/){
				#genechr is further along than sjchr
				$sjindex++; #increase sjindex by one
				redo INTRONFILE;
			} elsif ($genechr > $sjchr){
				$sjindex++; #increase sjindex by one
				redo INTRONFILE; #go to next intron line but keep same gene interval
			} else {
				next GENES;
			}
			
		} else {
			if ($sjstart < $genestart) {
				#sj is earlier in genome than current gene and therefor not fully contained in one of the 274 genes
				#my $intronline = join( "\t", $sjgene, $sjchr, $sjstart, $sjstop, $sjstrand );
				#push @labeledintrons, $intronline; #add entry to labeled introns array for use in next section
				$sjindex++; #increase sjindex by one
				redo INTRONFILE; #go to next intron line but keep same gene interval
			} elsif ($sjstart > $genestop) {
				#sj starts beyond the end of current gene, so we need to go to next gene
				next GENES;
			} else {
				if ($sjstop <= $genestop) {
					if (($sjstrand ne $genestrand) and ($sjstrand ne '0')) {
						#strand doesn't match so need to move on to next SJ
						my $intronline = join( "\t", $sjgene, $transcript, $type, $sjstrand, $sjchr, $sjstart, $sjstop, $samplecount );
						push @labeledintrons, $intronline; #add entry to labeled introns array for use in next section
						$sjindex++; #increase sjindex by one
						redo INTRONFILE; #go to next intron line but keep same gene interval
					} else {
						#everything matches so we label this SJ with current gene name
						my $intronline = join( "\t", $geneid, $transcript, $type, $sjstrand, $sjchr, $sjstart, $sjstop, $samplecount );
						push @labeledintrons, $intronline; #add entry to labeled introns array for use later
						$sjindex++; #increase sjindex by one
						redo INTRONFILE; #go to next intron line but keep same gene interval
					}
				} else {
					#sj ends beyond current gene and therefor not fully contained in one of the 274
					#my $intronline = join( "\t", $sjgene, $sjchr, $sjstart, $sjstop, $sjstrand );
					#push @labeledintrons, $intronline; #add entry to labeled introns array for use in next section
					$sjindex++; #increase sjindex by one
					redo INTRONFILE; #go to next intron line but keep same gene interval
				}
			}
		}
		
	}
}
###
#this section goes through the condensed introns file and puts the introns within any of the panel genes into a new array
my @nmdgtf;
my $intronindex = 0;
###
NMDGENES: foreach my $currentgene (@genes){
	chomp $currentgene;
	my ($geneid, $genestrand, $genechr, $genestart, $genestop) = split('\t', $currentgene);
	$genechr =~ s/chr//;
	INTRONGTF:
	{
		if (! defined $gtfintrons[$intronindex]){
			next NMDGENES;
		}
		chomp $gtfintrons[$intronindex];
		my ($introngene, $introntranscript, $introntype, $intronstrand, $intronchr, $intronstart, $intronstop) = split('\t', $gtfintrons[$intronindex]);
		$intronchr =~ s/chr//;
		if ($intronchr !~ /^$genechr$/) {
			if ($genechr eq 'Y' && $intronchr eq 'X'){
				$intronindex++;
				redo INTRONGTF;
			} elsif ($genechr eq 'X' && $intronchr eq 'Y'){
				next NMDGENES;
			} elsif ($genechr == 22 && $intronchr eq 'X'){
				next NMDGENES;
			} elsif ($genechr eq 'X' && $intronchr == 22){
				$intronindex++;
				redo INTRONGTF;
			} elsif ($genechr > $intronchr){
				$intronindex++;
				redo INTRONGTF;
			} else {
				next NMDGENES;
			}
		} else {
			if ($intronstart < $genestart) {
				#intron entry occurs before current gene
				$intronindex++;
				redo INTRONGTF;
			} else {
				if ($intronstart > $genestop){
					#intron entry occurs after current gene
					next NMDGENES;
				} else {
					if ($intronstop <= $genestop){
						#intron entry is within gene
						my $count = 0;
						my $nmdgtfline = join("\t", $introngene, $introntranscript, $introntype, $intronstrand, $intronchr, $intronstart, $intronstop, $count);
						push @nmdgtf, $nmdgtfline;
						$intronindex++;
						redo INTRONGTF;
					} else {
						$intronindex++;
						redo INTRONGTF;
					}
				}
			}
		}
	}
}
###
#This section creates an array of all possible splice sites among all the samples and the annotation file
###
my @annotated;
my @unannotated;
my $innerindex = 0;
my $outer;
my $inner;
my $nmd;
if ($#nmdgtf >= $#labeledintrons){
	$outer = \@nmdgtf;
	$inner = \@labeledintrons;
	$nmd = 1;
} else {
	$outer = \@labeledintrons;
	$inner = \@nmdgtf;
	$nmd = 2;
}
open my $output2, '>', $outputfile2 or die "Can't read $outputfile2 : $!";
###
OUTER: foreach my $intron1 ( @{$outer} ){
	chomp $intron1;
	my ($gene1, $transcript1, $type1, $strand1, $chr1, $start1, $stop1, $count1) = split('\t', $intron1);
	INNER:
	{
		if (! defined ${$inner}[$innerindex]) {
			my $gtfline = join( "\t", $gene1, $transcript1, $strand1, $chr1, $start1, $stop1, $count1 );
			print $output2 "$gtfline\n";
			if ($transcript1 eq 'UNANNOTATED') {
				push @unannotated, $gtfline;
			} else {
				push @annotated, $gtfline;
			}
			next OUTER;
		}
		chomp ${$inner}[$innerindex];
		my ($gene2, $transcript2, $type2, $strand2, $chr2, $start2, $stop2, $count2) = split('\t', ${$inner}[$innerindex]);
		if ($chr1 !~ /^$chr2$/){
			if ($chr1 eq 'Y' && $chr2 eq 'X'){
				my $gtfline = join( "\t", $gene2, $transcript2, $strand2, $chr2, $start2, $stop2, $count2 );
				print $output2 "$gtfline\n";
				if ($transcript2 eq 'UNANNOTATED') {
					push @unannotated, $gtfline;
				} else {
					push @annotated, $gtfline;
				}
				$innerindex++;
				redo INNER;
			} elsif ($chr1 eq 'X' && $chr2 eq 'Y'){
				my $gtfline = join( "\t", $gene1, $transcript1, $strand1, $chr1, $start1, $stop1, $count1 );
				print $output2 "$gtfline\n";
				if ($transcript1 eq 'UNANNOTATED') {
					push @unannotated, $gtfline;
				} else {
					push @annotated, $gtfline;
				}
				next OUTER;
			} elsif ($chr1 == 22 && $chr2 eq 'X'){
				my $gtfline = join( "\t", $gene1, $transcript1, $strand1, $chr1, $start1, $stop1, $count1 );
				print $output2 "$gtfline\n";
				if ($transcript1 eq 'UNANNOTATED') {
					push @unannotated, $gtfline;
				} else {
					push @annotated, $gtfline;
				}
				next OUTER;
			} elsif ($chr1 eq 'X' && $chr2 == 22){
				my $gtfline = join( "\t", $gene2, $transcript2, $strand2, $chr2, $start2, $stop2, $count2 );
				print $output2 "$gtfline\n";
				if ($transcript2 eq 'UNANNOTATED') {
					push @unannotated, $gtfline;
				} else {
					push @annotated, $gtfline;
				}
				$innerindex++;
				redo INNER;
			} elsif ($chr1 > $chr2){
				my $gtfline = join( "\t", $gene2, $transcript2, $strand2, $chr2, $start2, $stop2, $count2 );
				print $output2 "$gtfline\n";
				if ($transcript2 eq 'UNANNOTATED') {
					push @unannotated, $gtfline;
				} else {
					push @annotated, $gtfline;
				}
				$innerindex++;
				redo INNER;
			} else {
				my $gtfline = join( "\t", $gene1, $transcript1, $strand1, $chr1, $start1, $stop1, $count1 );
				print $output2 "$gtfline\n";
				if ($transcript1 eq 'UNANNOTATED') {
					push @unannotated, $gtfline;
				} else {
					push @annotated, $gtfline;
				}
				next OUTER;
			}
		} else {
			if ($start1 != $start2){
				if ($start2 > $start1){
					#inner is further along than outer
					my $gtfline = join( "\t", $gene1, $transcript1, $strand1, $chr1, $start1, $stop1, $count1 );
					print $output2 "$gtfline\n";
					if ($transcript1 eq 'UNANNOTATED') {
						push @unannotated, $gtfline;
					} else {
						push @annotated, $gtfline;
					}
					next OUTER;
				} else {
					#outer is further along than inner
					my $gtfline = join( "\t", $gene2, $transcript2, $strand2, $chr2, $start2, $stop2, $count2 );
					print $output2 "$gtfline\n";
					if ($transcript2 eq 'UNANNOTATED') {
						push @unannotated, $gtfline;
					} else {
						push @annotated, $gtfline;
					}
					$innerindex++;
					redo INNER;
				}
			} else {
				if ($stop1 != $stop2){
					if ($stop2 > $stop1){
						#inner is further along than outer
						my $gtfline = join( "\t", $gene1, $transcript1, $strand1, $chr1, $start1, $stop1, $count1 );
						print $output2 "$gtfline\n";
						if ($transcript1 eq 'UNANNOTATED') {
							push @unannotated, $gtfline;
						} else {
							push @annotated, $gtfline;
						}
						next OUTER;
					} else {
						#outer is further along than inner
						my $gtfline = join( "\t", $gene2, $transcript2, $strand2, $chr2, $start2, $stop2, $count2 );
						print $output2 "$gtfline\n";
						if ($transcript2 eq 'UNANNOTATED') {
							push @unannotated, $gtfline;
						} else {
							push @annotated, $gtfline;
						}
						$innerindex++;
						redo INNER;
					}
				} elsif ($strand1 ne $strand2) {
					my $gtfline = join( "\t", $gene1, $transcript1, $strand1, $chr1, $start1, $stop1, $count1 );
					my $intronline = join( "\t", $gene2, $transcript2, $strand2, $chr2, $start2, $stop2, $count2 );
					if ($transcript1 eq 'UNANNOTATED') {
						print $output2 "$intronline\n";
						print $output2 "$gtfline\n";
						push @unannotated, $gtfline;
						push @annotated, $intronline;
					} else {
						print $output2 "$gtfline\n";
						print $output2 "$intronline\n";
						push @unannotated, $intronline;
						push @annotated, $gtfline;
					}
					$innerindex++;
					next OUTER;
				} else {
					#everything matches so get info from nmdgtf but counts from labeledintrons
					if ($nmd == 1){
						my $gtfline = join( "\t", $gene1, $transcript1, $strand1, $chr1, $start1, $stop1, $count2 );
						print $output2 "$gtfline\n";
						push @annotated, $gtfline;
						$innerindex++;
						next OUTER;
					} else {
						my $gtfline = join( "\t", $gene2, $transcript2, $strand2, $chr2, $start2, $stop2, $count1 );
						print $output2 "$gtfline\n";
						push @annotated, $gtfline;
						$innerindex++;
						next OUTER;
					}
				}
			}
		}
	}
}
close $output2;
#this section checks each unannotated splice junction for a shared splice start or stop in the annotated junctions and creates an array of unannotated splice junctions that gives a ratio in comparison to the highest shared annotated splice. 
my @percents;
UNANNOTATED: foreach my $currentunannotated (@unannotated){
	chomp $currentunannotated;
	my ($geneid, $transcript, $strand, $chr, $start, $stop, $rawcount) = split('\t', $currentunannotated);
	if ($geneid eq 'UNKNOWN'){
		my $entry = join( "\t", $geneid, $transcript, $strand, $chr, $start, $stop, $rawcount);
		push @percents, $entry;
		next UNANNOTATED;
	}
	my $annotatedindex = 0;
	my $topcount = 0;
	#$chr =~ s/chr//;
	ANNOTATED:
	{
		if (! defined $annotated[$annotatedindex]){
			if ($topcount == 0){
				my $entry = join( "\t", $geneid, $transcript, $strand, $chr, $start, $stop, $rawcount);
				push @percents, $entry;	
			} else {
				my $ratio = $rawcount / $topcount;
				my $entry = join( "\t", $geneid, $transcript, $strand, $chr, $start, $stop, $ratio);
				push @percents, $entry;	
			}
			next UNANNOTATED;
		}
		chomp $annotated[$annotatedindex];
		my ($geneidA, $transcriptA, $strandA, $chrA, $startA, $stopA, $countA) = split('\t', $annotated[$annotatedindex]);
		if ($geneidA ne $geneid) {
			$annotatedindex++;
			redo ANNOTATED;
		} elsif ($strandA ne $strand) {
			$annotatedindex++;
			redo ANNOTATED;
		} elsif (($startA == $start) || ($stopA == $stop)) {
			if ($countA > $topcount) {
				$topcount = $countA;
				$annotatedindex++;
				redo ANNOTATED;
			} else {
				$annotatedindex++;
				redo ANNOTATED;
			}
		} else {
			$annotatedindex++;
			redo ANNOTATED;
		}
	}
}
#This section combines the annotated and unannotated (ratios) junctions arrays in the correct order and prints to the output file
###
my $innerindex2 = 0;
my $outer2;
my $inner2;
my $nmd2;
if ($#annotated >= $#percents){
	$outer2 = \@annotated;
	$inner2 = \@percents;
	$nmd2 = 1;
} else {
	$outer2 = \@percents;
	$inner2 = \@annotated;
	$nmd2 = 2;
}
open my $output, '>', $outputfile or die "Can't read $outputfile : $!";
###
OUTER: foreach my $entry1 ( @{$outer2} ){
	chomp $entry1;
	my ($gene1, $transcript1, $strand1, $chr1, $start1, $stop1, $count1) = split('\t', $entry1);
	INNER:
	{
		if (! defined ${$inner2}[$innerindex2]) {
			my $gtfline = join( "\t", $gene1, $transcript1, $strand1, $chr1, $start1, $stop1, $count1 );
			print $output "$gtfline\n";
			next OUTER;
		}
		chomp ${$inner2}[$innerindex2];
		my ($gene2, $transcript2, $strand2, $chr2, $start2, $stop2, $count2) = split('\t', ${$inner2}[$innerindex2]);
		if ($chr1 !~ /^$chr2$/){
			if ($chr1 eq 'Y' && $chr2 eq 'X'){
				my $gtfline = join( "\t", $gene2, $transcript2, $strand2, $chr2, $start2, $stop2, $count2 );
				print $output "$gtfline\n";
				$innerindex2++;
				redo INNER;
			} elsif ($chr1 eq 'X' && $chr2 eq 'Y'){
				my $gtfline = join( "\t", $gene1, $transcript1, $strand1, $chr1, $start1, $stop1, $count1 );
				print $output "$gtfline\n";
				next OUTER;
			} elsif ($chr1 == 22 && $chr2 eq 'X'){
				my $gtfline = join( "\t", $gene1, $transcript1, $strand1, $chr1, $start1, $stop1, $count1 );
				print $output "$gtfline\n";
				next OUTER;
			} elsif ($chr1 eq 'X' && $chr2 == 22){
				my $gtfline = join( "\t", $gene2, $transcript2, $strand2, $chr2, $start2, $stop2, $count2 );
				print $output "$gtfline\n";
				$innerindex2++;
				redo INNER;
			} elsif ($chr1 > $chr2){
				my $gtfline = join( "\t", $gene2, $transcript2, $strand2, $chr2, $start2, $stop2, $count2 );
				print $output "$gtfline\n";
				$innerindex2++;
				redo INNER;
			} else {
				my $gtfline = join( "\t", $gene1, $transcript1, $strand1, $chr1, $start1, $stop1, $count1 );
				print $output "$gtfline\n";
				next OUTER;
			}
		} else {
			if ($start1 != $start2){
				if ($start2 > $start1){
					#inner is further along than outer
					my $gtfline = join( "\t", $gene1, $transcript1, $strand1, $chr1, $start1, $stop1, $count1 );
					print $output "$gtfline\n";
					next OUTER;
				} else {
					#outer is further along than inner
					my $gtfline = join( "\t", $gene2, $transcript2, $strand2, $chr2, $start2, $stop2, $count2 );
					print $output "$gtfline\n";
					$innerindex2++;
					redo INNER;
				}
			} else {
				if ($stop1 != $stop2){
					if ($stop2 > $stop1){
						#inner is further along than outer
						my $gtfline = join( "\t", $gene1, $transcript1, $strand1, $chr1, $start1, $stop1, $count1 );
						print $output "$gtfline\n";
						next OUTER;
					} else {
						#outer is further along than inner
						my $gtfline = join( "\t", $gene2, $transcript2, $strand2, $chr2, $start2, $stop2, $count2 );
						print $output "$gtfline\n";
						$innerindex2++;
						redo INNER;
					}
				} elsif ($strand1 ne $strand2) {
					my $gtfline = join( "\t", $gene1, $transcript1, $strand1, $chr1, $start1, $stop1, $count1 );
					my $intronline = join( "\t", $gene2, $transcript2, $strand2, $chr2, $start2, $stop2, $count2 );
					if ($nmd2 == 1){
						print $output "$gtfline\n";
						print $output "$intronline\n";
					} else {
						print $output "$intronline\n";
						print $output "$gtfline\n";
					}
					$innerindex2++;
					next OUTER;
				} else {
					my $gtfline = join( "\t", $gene1, $transcript1, $strand1, $chr1, $start1, $stop1, $count1 );
					my $intronline = join( "\t", $gene2, $transcript2, $strand2, $chr2, $start2, $stop2, $count2 );
					if ($nmd2 == 1){
						print $output "$gtfline\n";
						print $output "$intronline\n";
					} else {
						print $output "$intronline\n";
						print $output "$gtfline\n";
					}
					$innerindex2++;
					next OUTER;
				}
			}
		}
	}
}
###
close $output;