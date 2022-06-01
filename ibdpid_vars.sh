#!/bin/bash
#variables specific to ibdpid panel
humangenome=~/scratch/gencode.v29 #base folder name
genomeDir=$humangenome/starindex #STAR index
wholeGenomeFasta=$humangenome/GRCh38.primary_assembly.genome.fa
annotations=$humangenome/gencode.v29.primary_assembly.annotation.gtf
altsplices=~/scratch/Homo_sapiens_jun12/NCBI/GRCh38/Annotation/Genes/aberrant_splices.txt
exonparts=$humangenome/ibd.pid.genes/ibdpid.exonic_parts.vsort.forPSI.gff
###
genes_for_counts=$humangenome/ibd.pid.genes/ibdpid.hgnc.ensg.txt
intervals=$humangenome/ibd.pid.genes/ibdpidgenes.intervals
sortedintervals=$humangenome/ibd.pid.genes/ibdpid.collapsed.genes.intervals.hgnc.txt
introns=$humangenome/condensed_introns_nooverlap_hgnconly.txt
dexseq_gff=$humangenome/ibd.pid.genes/ibdpid.genes.gencode.v29.pri.annot.collapsed.gff
trapscores=$humangenome/ibd.pid.genes/ibdpid.hg38.trapscores.txt

