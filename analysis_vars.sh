#!/bin/bash
user=kberger9

basedir=~/p-ggibson3-0/rich_project_pb1/targetedRNA


tools=~/p-ggibson3-0/rich_project_pb1/tools
scripts=~/p-ggibson3-0/rich_project_pb1/scripts

ibdpid_vars=$scripts/ibdpid_vars.sh
nmd_vars=$scripts/nmd_vars.sh
mousenmd_vars=$scripts/mousenmd_vars.sh
krabbe_vars=$scripts/krabbe_vars.sh

fastqc=$tools/fastqc/fastqc
trimmomatic=$tools/trimmomatic/trimmomatic-0.36.jar
adapters=$tools/trimmomatic/adapters/TruSeq3-PE-2.fa
STAR=$tools/star/v2.6.1d/bin/Linux_x86_64/STAR
bedtools=$tools/bedtools2/bin
qorts=$tools/qorts/QoRTs-STABLE.jar
picard=$tools/picard/picard.jar
dexseq_count=/storage/home/hcoda1/7/kberger9/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py
gatk=$tools/gatk/GenomeAnalysisTK.jar
phase1snps=~/p-ggibson3-0/rich_project_pb1/homo_sapiens/gatk_resource_bundle/1000G_phase1.snps.high_confidence.hg38.vcf
millsIndels=~/p-ggibson3-0/rich_project_pb1/homo_sapiens/gatk_resource_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf
dbsnp=~/p-ggibson3-0/rich_project_pb1/homo_sapiens/gatk_resource_bundle/dbsnp_146.hg38.vcf
qortsscript=/storage/home/hcoda1/7/kberger9/R/x86_64-pc-linux-gnu-library/3.6/QoRTs/extdata/scripts/qortsGenMultiQC.R
annovar=$tools/annovar/table_annovar.pl

psi_script=$scripts/psi.sh
genecounts_script=$scripts/genecounts.pl
genecounts_nodupsscript=$scripts/genecounts.nodups.pl
ind_splicing=$scripts/persample_splicing.pl

joint_splicing=$scripts/joint_splicing.pl
process_vars=$scripts/process_variants.pl
zscores_script=$scripts/zscores.R

pbshead_big="#PBS -A GT-ggibson3-biocluster\n#PBS -l nodes=1:ppn=12\n#PBS -l mem=40gb\n#PBS -l walltime=12:00:00\n#PBS -j oe\n\ncd \$PBS_O_WORKDIR"
pbshead_little="#PBS -A GT-ggibson3-biocluster\n#PBS -l nodes=1:ppn=4\n#PBS -l mem=8gb\n#PBS -l walltime=12:00:00\n#PBS -j oe\n\ncd \$PBS_O_WORKDIR"

####old
#STAR=$tools/star/bin/Linux_x86_64/STAR
homedir="/Users/kieraberger/Dropbox (GaTech)/nmd_analysis"
#user=kberger9
#basedir=~/scratch/nmd_samples_jun12
#tools=~/data/tools
#scripts=~/data/scripts
#humangenome=~/scratch/Homo_sapiens_jun12

#genomeDir=~/scratch/genomeDir_jun12/illumina_grch38_aug2017
#wholeGenomeFasta=$humangenome/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa
#annotations=$humangenome/NCBI/GRCh38/Annotation/Genes/genes.gtf
#altsplices=$humangenome/NCBI/GRCh38/Annotation/Genes/aberrant_splices.txt
#trapscores=$humangenome/NCBI/GRCh38/Annotation/Genes/nmd274genes.hg38.trapscores.txt
#exonparts=$humangenome/NCBI/GRCh38/Annotation/Genes/ncbi274.exonic_parts.forPSI.gff
leafexons=$humangenome/NCBI/GRCh38/Annotation/Genes/ncbi_grch38.leafexons.txt.gz

#fastqc=$tools/fastqc/fastqc


#trimmomatic=$tools/trimmomatic/trimmomatic-0.36.jar
#adapters=$tools/trimmomatic/adapters/TruSeq3-PE-2.fa
#qorts=$tools/qorts/QoRTs-STABLE.jar
#picard=$tools/picard/picard.jar
#bedtools=$tools/bedtools2/bin
leafcutter=$tools/leafcutter/scripts

#gatk=$tools/gatk/GenomeAnalysisTK.jar
#phase1snps=$humangenome/gatk_resource_bundle/1000G_phase1.snps.high_confidence.hg38.vcf
#millsIndels=$humangenome/gatk_resource_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf
#dbsnp=$humangenome/gatk_resource_bundle/dbsnp_146.hg38.vcf
#annovar=$tools/annovar/table_annovar.pl

#qortsscript=/nv/hp10/kberger9/R/x86_64-pc-linux-gnu-library/3.3/QoRTs/extdata/scripts/qortsGenMultiQC.R
tissuecheck_script=$tools/macarthur/MendelianRNA-seq-master/QC/MuscleCheck.R

#genecounts_script=$scripts/genecounts_strandedonly.pl
#genecounts_nodupsscript=$scripts/genecounts.nodups.pl
#dexseq_count=/nv/hp10/kberger9/R/x86_64-pc-linux-gnu-library/3.3/DEXSeq/python_scripts/dexseq_count.py
#dexseq_gff=$humangenome/NCBI/GRCh38/Annotation/Genes/genes.onlynormchr.noxm.gff
#psi_script=$scripts/psi.sh
#ind_splicing=$scripts/persample_splicing.pl
#joint_splicing=$scripts/joint_splicing.pl
#process_vars=$scripts/process_variants.pl

peutils=$tools/myenv/bin/pe_utils
constexons=$humangenome/NCBI/GRCh38/Annotation/Genes/exons/genes.ncbi274.min_1000.const_exons.gff

#intervals=~/data/target_coordinates/nmd274genes.intervals
#sortedintervals=~/data/target_coordinates/nmd274genes_intervals_vsort.txt
#geneinfo=~/data/target_coordinates/nmd274genes_gtf106_ensg.txt
musclegvcfs=~/p-ggibson3-0/rich_project_pb1/variants/muscle_gvcfs.list
bloodgvcfs=~/p-ggibson3-0/rich_project_pb1/variants/wholeblood_gvcfs.list
#introns=~/data/target_coordinates/condensed_introns_nooverlap.txt


#pbshead_big="#PBS -l nodes=1:ppn=12\n#PBS -l mem=40gb\n#PBS -l walltime=12:00:00\n#PBS -q iw-shared-6\n#PBS -j oe\n\ncd \$PBS_O_WORKDIR"
#pbshead_little="#PBS -l nodes=1:ppn=4\n#PBS -l mem=8gb\n#PBS -l walltime=12:00:00\n#PBS -q iw-shared-6\n#PBS -j oe\n\ncd \$PBS_O_WORKDIR"
