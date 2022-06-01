#!/usr/bin/perl -w
use strict;
use warnings;
###
my $annotatedvcf = "$ARGV[0]";
open my $vcfinput, '<', $annotatedvcf or die "Can't read $annotatedvcf : $!";
my @vcf = <$vcfinput>;
close $vcfinput;
###
my $trapfile = "$ARGV[1]";
open my $trapscores, '<', $trapfile or die "Can't read $trapfile : $!";
my @trap= <$trapscores>;
close $trapscores;
###
my $outputfile = "$ARGV[2]";
chomp $outputfile;
open my $output, '>', $outputfile or die "Can't read $outputfile : $!";
###
my @finalvcf;
#remove header from annotated vcf before processing, insert trap score column after avsnp column (13th), put header into final vcf array
my $header = shift @vcf;
chomp $header;
my @headsplit = split('\t', $header);
my $newhead = join("\t", @headsplit[0..12], "TraP", @headsplit[13..($#headsplit)]);
print $output $newhead."\n";
#print STDERR "$newhead\n";
#push @finalvcf, $newhead;
###
my $noscore = ".";
my $trapindex = 0;
my $vartype = "HET";
VCF: foreach my $variant (@vcf) {
	#process each column, if column is ".." then insert the tab, recognize the GT:AD (113 with old annovar, now 122,  - don't count TRAP yet) column as last before samples. check output from yesterday to see if this is the issue
	next if ($variant =~ m/^Chr/); 
	chomp $variant;
	my @splitline = split('\t', $variant);
	my $final = join("\t", @splitline[0..12]);
	my $chr = $splitline[0]; #add column variable
	$chr =~ s/chr//;
	my $start = $splitline[1];
	my $stop = $splitline[2];
	my $ref = $splitline[3];
	my $alt = $splitline[4];
	my $gene = $splitline[6];
	my $finalsamples = $splitline[121]; #was 112 in old version of annovar
	foreach my $sample (@splitline[122..($#splitline)]) {
		chomp $sample;
		#print STDERR $sample;
		if ($sample =~ m/HET/) {
			#sample is heterozygous
			my @sampleinfo = split(':', $sample);
			my @alleles = split(',', $sampleinfo[1]);
			if ($#alleles > 1) {
				$finalsamples = join("\t", $finalsamples, "$sample:MULTI");
			} elsif ($alleles[1] < 10) {
				$finalsamples = join("\t", $finalsamples, "$sample:LowAD");
			} else {
				$finalsamples = join("\t", $finalsamples, $sample);
			}
			#print STDERR $finalsamples;
			#exit;
		} else {
			$finalsamples = join("\t", $finalsamples, $sample);
		}	
	}
	if ($start != $stop) {
		#variant is not a SNP, no trap score will be possible
		$final = join("\t", $final, $noscore, @splitline[13..120], $finalsamples);
		print $output "$final\n";
		#push @finalvcf, $final;
		next VCF;
	}
	if ($ref eq '-' || $alt eq '-' || length($ref) > 1 || length($alt) > 1) {
		#variant is not a SNP, no trap score will be possible
		$final = join("\t", $final, $noscore, @splitline[13..120], $finalsamples);
		print $output "$final\n";
		#push @finalvcf, $final;
		next VCF;
	}
	TRAP:
	{
		if (! defined $trap[$trapindex]){
			$final = join("\t", $final, $noscore, @splitline[13..120], $finalsamples);
			print $output "$final\n";
			#push @finalvcf, $final;
			next VCF;
		}
		chomp $trap[$trapindex];
		my ($trapchr, $trapsite, $trapref, $trapalt, $trapgene, $score) = split('\s+', $trap[$trapindex]);
		if ($trapchr eq "chr8_KI270821v1_alt"){
			$trapindex++; #increase TRAPindex by one
			redo TRAP; #go to trap line but keep same variant 
		}
		$trapchr =~ s/chr//;
		if ($chr ne $trapchr) {
			if ($chr eq 'M'){
				$final = join("\t", $final, $noscore, @splitline[13..120], $finalsamples);
				print $output "$final\n";
				#push @finalvcf, $final;
				next VCF;
			} elsif ($chr eq 'Y'  && $trapchr eq 'X'){
				$trapindex++; #increase TRAPindex by one
				redo TRAP; #go to trap line but keep same variant 
			} elsif ($chr eq 'X' && $trapchr eq 'Y'){
				#site does not have a matching trap score
				$final = join("\t", $final, $noscore, @splitline[13..120], $finalsamples);
				print $output "$final\n";
				#push @finalvcf, $final;
				next VCF;
			} elsif ($chr == 22 && $trapchr eq 'X'){
				#site does not have a matching trap score
				$final = join("\t", $final, $noscore, @splitline[13..120], $finalsamples);
				print $output "$final\n";
				#push @finalvcf, $final;
				next VCF;
			} elsif ($chr eq 'X' && $trapchr == 22){
				$trapindex++;
				redo TRAP; 
			} elsif ($chr > $trapchr){
				$trapindex++; 
				redo TRAP; 
			} else {
				#site does not have a matching trap score
				$final = join("\t", $final, $noscore, @splitline[13..120], $finalsamples);
				print $output "$final\n";
				#push @finalvcf, $final;
				next VCF;
			}
		} else {
			if ($stop < $trapsite) {
				#site does not have a matching trap score
				$final = join("\t", $final, $noscore, @splitline[13..120], $finalsamples);
				print $output "$final\n";
				#push @finalvcf, $final;
				next VCF;
			} elsif ($stop > $trapsite) {
				$trapindex++; 
				redo TRAP; 
			} else {
				if ($ref ne $trapref) {
					#site does not have a matching trap score
					$final = join("\t", $final, $noscore, @splitline[13..120], $finalsamples);
					print $output "$final\n";
					#push @finalvcf, $final;
					next VCF;
				} else {
					if ($alt ne $trapalt || $gene ne $trapgene) {
						$trapindex++; 
						redo TRAP;
					} else {
						#match!
						$final = join("\t", $final, $score, @splitline[13..120], $finalsamples);
						print $output "$final\n";
						#push @finalvcf, $final;
						$trapindex = $trapindex - 3; #back up a few spots so that if the next variant is the same position if can still get matched
						next VCF;
					}
				}
			}
		}
	}
}

#foreach (@finalvcf) {
 # print $output "$_\n";
 #}
close $output;