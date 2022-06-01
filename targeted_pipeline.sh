#!/bin/bash
scriptname=$0   # $0 is the name of the program
#make list of tools/paths/versions to echo onto a roster of files for each sample? Also list of required files? & modules?
panel=ibdpid #either nmd or ibdpid or krabbe
batch=""
trim=no
kmersize=35
readLength=151
vars=~/p-ggibson3-0/rich_project_pb1/scripts/analysis_vars.sh
###
#help function
HELP () {
	echo -e "\n\tUSAGE: $scriptname -b BATCH [-p (ibdpid|nmd|krabbe)] [-t (yes|no)] [-v VARIABLES_FILE]\n"
	echo -e "\tMANDATORY OPTIONS:\n"
	echo -e "\t\t-b BATCH\t\tName of the sequencing run, usually in the form of GG##\n\n"
	echo -e "\tADDITIONAL OPTIONS:\n"
	echo -e "\t\t-v VARIABLES_FILE\tLocation of file containing variable information for analysis scripts"
	echo -e "\t\t\t\t\tDEFAULT: ~/data/scripts/analysis_vars.sh\n"	
	echo -e "\t\t-p PANEL\t\tPanel used for Targeted RNA sequencing"  
	echo -e "\t\t\t\t\tDEFAULT: ibdpid\n"						
	echo -e "\t\t-t TRIM\t\t\tWhether the FASTQ files to be aligned have been trimmed"
	echo -e "\t\t\t\t\tDEFAULT: no\n"
	echo -e "\t\t-k KMER_SIZE\t\tkmer size to use in GATK HaplotypeCaller. Must be an integer."
	echo -e "\t\t\t\t\tDEFAULT: 35\n"
	echo -e "\t\t-l READ_LENGTH\t\tRead Length. Must be an integer."
	echo -e "\t\t\t\t\tDEFAULT: 151\n\n"
	echo -e "\tNOTE: this script expects PE reads with fastq filenames containing \"_L00#_R#_001.fastq.gz\""
	exit 0
}
###
#Check the number of arguments. If none are passed, print help and exit.
NUMARGS=$#
if [ $NUMARGS -eq 0 ]; then
  HELP
fi
###
#Get Options from command line
while getopts :k:t:l:b:p:v:h opt
do
	case $opt in
	t)trim=$OPTARG;;
	b)batch=$OPTARG;;
	p)panel=$OPTARG;;
	k)kmersize=$OPTARG;;
	l)readLength=$OPTARG;;
	v)vars=$OPTARG;;
    h) HELP; exit;;
    \?) echo "ERROR: Invalid option: -$OPTARG" >&2; echo -e "Run $scriptname -h to see usage and help" >&2; exit 1;;
    :) echo "ERROR: Option -$OPTARG requires an argument" >&2; echo -e "Run $scriptname -h to see usage and help" >&2; exit 1;;
    *) echo "ERROR: Unimplemented option: -$OPTARG" >&2; echo -e "Run $scriptname -h to see usage and help" >&2; exit 1;;
  esac
done
shift $(($OPTIND -1))
###
#check to make sure variables file exists
if [ ! -f $vars ] 
then
    echo -e "ERROR: File $vars DOES NOT exist" >&2
	echo -e "Edit this scripts default or use the -v option to point the script to your file that contains variables for analysis scripts" >&2
    exit 1
fi
source $vars
###
#check if panel is nmd or ibdpid, set panelid to appropriate name
panelid=""
if [ $panel == 'ibdpid' ]
then 
	panelid=ibdpid
elif [ $panel == 'nmd' ]
then
	panelid=nmd
elif [ $panel == 'krabbe' ]
then
	panelid=krabbe
else
	echo -e "ERROR: \"$panel\" is not a valid argument for [-p] option. Valid arguments are \"idbpid\" or \"nmd\" or \"krabbe\""
	exit 1
fi
#source panel variables
panelvars="${panelid}_vars"
if [ ! -f ${!panelvars} ] 
then
    echo -e "ERROR: File ${!panelvars} DOES NOT exist" >&2
	echo -e "Edit your file $panelvars to point to the correct ${panelid}_vars file" >&2
    exit 1
fi
source ${!panelvars}
###
#make sure kmer size is a positive integer
if ! [[ "$kmersize" =~ ^[0-9]+$ ]]
    then
        echo -e "ERROR: \"$kmersize\" is not a valid argument for [-k] option. kmer size must be an integer."
		exit 1
fi

if ! [[ "$readLength" =~ ^[0-9]+$ ]]
    then
        echo -e "ERROR: \"$readLength\" is not a valid argument for [-l] option. Read Length must be an integer."
		exit 1
fi
###
#check whether or not to trim, set fastq_dir appropriately
fastq_type=""
if [ $trim == 'yes' ]
then 
	fastq_type=fastq_files/trimmed
	echo "Will align TRIMMED fastq files"
elif [ $trim == 'no' ]
then
	fastq_type=fastq_files
	echo "Will align UNTRIMMED fastq files"
else
	echo -e "ERROR: \"$trim\" is not a valid argument for [-t] option. Valid arguments are \"yes\" or \"no\""
	exit 1
fi
###
#ensure directory containing samples folder and samples.txt exists
if [ ! -d $basedir/$panelid/$batch ] 
then
    echo -e "ERROR: Directory $basedir/$panelid/$batch DOES NOT exist" >&2
	echo -e "Check that option arguments -b (& -p if used) is correct" >&2
    exit 1
fi
###
#ensure samples.txt exists
if [ ! -f $basedir/$panelid/$batch/samples.txt ] 
then
    echo -e "ERROR: File $basedir/$panelid/$batch/samples.txt DOES NOT exist" >&2
	echo -e "Create the standard samples.txt file for this batch" >&2
    exit 1
fi
###
star_2pass () { 
	local sample=$1
	local dir=$2
	local readgroup=$3
	if [ ! -d $dir/star_2pass ]
	then
		mkdir $dir/star_2pass #STAR needs outFileNamePrefix to already exist before running
	fi
	local step=star
	local pbsfile=$dir/star_2pass/$step.pbs
	fastq_dir=$dir/$fastq_type
	rgline=""
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile #start pbs file for star run
	#get read 1s

	read1=""
	counter=1
	for file in $fastq_dir/*_R1_001.fastq.gz
	do
		if [ "$counter" -eq 1 ]
		then
			read1="${file}"
			counter=$((counter+1))
		else
			read1="${read1},${file}"
		fi
	done
	#reset and get R2s
	counter=1
	read2=""
	for file in $fastq_dir/*_R2_001.fastq.gz
	do
		lane=`echo "${file: -17}" | sed 's/_R2_001.fastq.gz//g'`
		if [ "$counter" -eq 1 ]
		then
			read2="${file}"
			rgline="ID:${readgroup}.${lane} PL:illumina PU:${readgroup}.${lane} LB:${panelid} SM:${sample}.${batch}"
			counter=$((counter+1))
		else
			read2="${read2},${file}"
			rgline="${rgline} , ID:${readgroup}.${lane} PL:illumina PU:${readgroup}.${lane} LB:${panelid} SM:${sample}.${batch}"
			counter=$((counter+1))
		fi
	done
	echo "$STAR --runThreadN 16 --genomeDir $genomeDir --sjdbGTFfile $annotations --sjdbFileChrStartEnd $altsplices --readFilesIn $read1 $read2 --readFilesCommand zcat --outFileNamePrefix $dir/star_2pass/ --sjdbOverhang 100 --outSAMtype BAM SortedByCoordinate --outSAMattrRGline $rgline --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMultimapNmax 1 --outSJfilterCountUniqueMin 5 5 5 5 --outSJfilterCountTotalMin 5 5 5 5 --outFilterType BySJout --outSJfilterReads Unique --alignInsertionFlush Right --quantMode GeneCounts --twopassMode Basic" >> $pbsfile #add commands for star run
	echo -e "mv $dir/star_2pass/Aligned.sortedByCoord.out.bam $dir/star_2pass/$sample.bam" >> $pbsfile
	echo -e "mv $dir/star_2pass/Log.final.out $dir/star_2pass/$sample.$batch.log.final.out" >> $pbsfile
	echo -e "mv $dir/star_2pass/SJ.out.tab $dir/star_2pass/$sample.$batch.sj.out.tab" >> $pbsfile
	echo -e "samtools index $dir/star_2pass/$sample.bam" >> $pbsfile
	echo -e "qsub $dir/star_2pass/counts_tissuecheck.pbs\n" >> $pbsfile
	echo -e "qsub $dir/star_2pass/qorts.pbs\n" >> $pbsfile
	echo -e "cd $dir/psi/" >> $pbsfile 
	echo -e "qsub $dir/psi/psi.pbs" >> $pbsfile
	echo -e "cd $dir/splicing/" >> $pbsfile 
	echo -e "qsub $dir/splicing/removedups.pbs" >> $pbsfile
	echo -e "cd $dir/variant_calling/" >> $pbsfile #move into variant_calling dir so pbs outfile goes to here
	echo -e "qsub $dir/variant_calling/markdups.pbs\n" >> $pbsfile #submit next pbs file in chain
}

psi () {
	local sample=$1
	local dir=$2
	local step=psi
	if [ ! -d $dir/psi ]
	then
		mkdir $dir/psi
	fi
	local pbsfile=$dir/psi/$step.pbs
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile #start pbs file for step
	echo -e "$psi_script -s $sample -b $batch -d $dir -l $readLength -v $vars -p ${!panelvars}" >> $pbsfile
	echo -e "sed -i \"1i #\t#\t$sample\t$sample\t$sample.$batch\" $basedir/$panelid/$batch/samples/$sampledir/psi/$sample.$batch.exonic_parts.psi" >> $pbsfile
}

readcounts_tissuecheck () {
	local sample=$1
	local dir=$2
	local step=counts_tissuecheck
	if [ ! -d $dir/genecounts ]
	then
		mkdir $dir/genecounts
	fi
	local pbsfile=$dir/star_2pass/$step.pbs
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_little" >> $pbsfile #start pbs file for step
	echo -e "$genecounts_script $sample $dir $genes_for_counts $batch" >> $pbsfile
	echo -e "ontarget=\`awk '{sum += \$2} END {print sum}' $dir/genecounts/$sample.$batch.raw.genecounts.txt\`; total=\`head -n 9 $dir/star_2pass/$sample.$batch.log.final.out | tail -n 1 | awk '{print \$6}'\`; awk -v tot=\"\$total\" -v ont=\"\$ontarget\" 'BEGIN{print \"Reads on target: \"ont/tot}'" >> $pbsfile
	#echo -e "Rscript $tissuecheck_script -patient_rpkm $dir/genecounts/$sample.genes.rpkm.gct -out_file $basedir/$panelid/$batch/metrics/tissuechecks/$sample.tissuecheck" >> $pbsfile
}

qorts () {
	local sample=$1
	local dir=$2
	local readgroup=$3
	local step=qorts
	local pbsfile=$dir/star_2pass/$step.pbs
	fastq_dir=$dir/$fastq_type
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile #start pbs file for step

	for file in $fastq_dir/*_R2_001.fastq.gz
	do
		lane=`echo "${file: -17}" | sed 's/_R2_001.fastq.gz//g'`
		if [ ! -d $basedir/$panelid/$batch/metrics/qorts/${sample}.${readgroup}.${lane} ]
		then
			mkdir $basedir/$panelid/$batch/metrics/qorts/${sample}.${readgroup}.${lane}
		fi
		echo -e "java -jar $qorts QC --verbose --stranded --readGroup ${readgroup}.${lane} $dir/star_2pass/${sample}.bam $annotations $basedir/$panelid/$batch/metrics/qorts/${sample}.${readgroup}.${lane}/" >> $pbsfile
	done
}

remove_dups () {
	local sample=$1
	local dir=$2
	local step=removedups
	local pbsfile=$dir/splicing/$step.pbs
	if [ ! -d $dir/splicing ]
	then
		mkdir $dir/splicing
	fi
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile #start pbs file for step

	if [[ -f $dir/splicing/$sample.nodups.R*.fastq.gz ]]
	then
		rm $dir/splicing/$sample.nodups.R*.fastq.gz
	fi
	echo -e "java -jar $picard MarkDuplicates I=$dir/star_2pass/$sample.bam O=$dir/splicing/$sample.nodups.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true REMOVE_SEQUENCING_DUPLICATES=true M=$dir/splicing/removedups_metrics.txt" >> $pbsfile
	echo -e "samtools view -h -u -f 0x2 $dir/splicing/$sample.nodups.bam | samtools sort -n - | samtools fastq -1 $dir/splicing/$sample.nodups.R1.fastq -2 $dir/splicing/$sample.nodups.R2.fastq -" >> $pbsfile
	echo -e "gzip $dir/splicing/$sample.nodups.R1.fastq" >> $pbsfile
	echo -e "gzip $dir/splicing/$sample.nodups.R2.fastq" >> $pbsfile
	echo -e "cd $dir/splicing/" >> $pbsfile
	echo -e "qsub $dir/splicing/star_remap.pbs" >> $pbsfile
}

splicing_remap () {
	local sample=$1
	local dir=$2
	local step=star_remap
	local pbsfile=$dir/splicing/$step.pbs
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile #start pbs file for star run
	echo "$STAR --runThreadN 16 --genomeDir $genomeDir --sjdbGTFfile $annotations --sjdbFileChrStartEnd $altsplices --readFilesIn $dir/splicing/$sample.nodups.R1.fastq.gz $dir/splicing/$sample.nodups.R2.fastq.gz --readFilesCommand zcat --outFileNamePrefix $dir/splicing/ --sjdbOverhang 100 --outSAMtype BAM SortedByCoordinate --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMultimapNmax 1 --outSJfilterCountUniqueMin 5 5 5 5 --outSJfilterCountTotalMin 5 5 5 5 --outFilterType BySJout --outSJfilterReads Unique --alignInsertionFlush Right --quantMode GeneCounts --twopassMode Basic" >> $pbsfile #add commands for star run
	echo -e "mv $dir/splicing/Aligned.sortedByCoord.out.bam $dir/splicing/$sample.remap.nodups.bam" >> $pbsfile
	echo -e "samtools index $dir/splicing/$sample.remap.nodups.bam" >> $pbsfile
	echo -e "mv $dir/splicing/Log.final.out $dir/splicing/$sample.$batch.nodups.log.final.out" >> $pbsfile
	echo -e "mv $dir/splicing/SJ.out.tab $dir/splicing/$sample.$batch.nodups.sj.out.tab" >> $pbsfile
	echo -e "cd $dir/splicing/" >> $pbsfile
	echo -e "qsub $dir/splicing/splice_counts.pbs" >> $pbsfile
	echo -e "qsub $dir/splicing/counts_nodups.pbs" >> $pbsfile
	echo -e "qsub $dir/splicing/qortsremap.pbs" >> $pbsfile
	echo -e "cd $dir/exoncounts/" >> $pbsfile
	echo -e "qsub $dir/exoncounts/exoncounts.nodups.pbs" >> $pbsfile
}

qorts_remap () {
	local sample=$1
	local dir=$2
	local readgroup=$3
	local step=qortsremap
	local pbsfile=$dir/splicing/$step.pbs
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile #start pbs file for step

	if [ ! -d $basedir/$panelid/$batch/metrics/qorts/${sample}.${readgroup}.dedup ]
	then
		mkdir $basedir/$panelid/$batch/metrics/qorts/${sample}.${readgroup}.dedup
	fi
	echo -e "java -jar $qorts QC --verbose --stranded $dir/splicing/${sample}.remap.nodups.bam $annotations $basedir/$panelid/$batch/metrics/qorts/${sample}.${readgroup}.dedup/" >> $pbsfile
}

readcounts_nodups () {
	local sample=$1
	local dir=$2
	local step=counts_nodups
	local pbsfile=$dir/splicing/$step.pbs
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_little" >> $pbsfile #start pbs file for step
	echo -e "$genecounts_nodupsscript $sample $dir $genes_for_counts $batch" >> $pbsfile
	echo -e "ontarget=\`awk '{sum += \$2} END {print sum}' $dir/genecounts/$sample.$batch.nodups.raw.genecounts.txt\`; total=\`head -n 9 $dir/splicing/$sample.$batch.nodups.log.final.out | tail -n 1 | awk '{print \$6}'\`; awk -v tot=\"\$total\" -v ont=\"\$ontarget\" 'BEGIN{print \"Reads on target: \"ont/tot}'" >> $pbsfile
}

norm_label_splices () {
	local sample=$1
	local dir=$2
	local step=splice_counts
	local pbsfile=$dir/splicing/$step.pbs
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_little" >> $pbsfile #start pbs file for star run
	echo -e "$ind_splicing $dir $sample $sortedintervals $introns $batch" >> $pbsfile
}

mark_duplicates () {
	local sample=$1
	local dir=$2
	local step=markdups
	if [ ! -d $dir/variant_calling ]
	then
		mkdir $dir/variant_calling
	fi
	if [ ! -d $dir/variant_calling/temp ]
	then
		mkdir $dir/variant_calling/temp
	fi
	local pbsfile=$dir/variant_calling/$step.pbs
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile #start pbs file for step
	echo -e "java -jar $picard MarkDuplicates I=$dir/star_2pass/$sample.bam O=$dir/variant_calling/temp/$sample.$batch.md.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=false REMOVE_SEQUENCING_DUPLICATES=true M=$dir/variant_calling/markdups_metrics.txt" >> $pbsfile
	if [ $panelid == 'nmd' ]
	then
		echo -e "samtools view -h -u $dir/variant_calling/temp/$sample.$batch.md.bam chr2:71450000-71690000 | samtools sort -o $dir/variant_calling/$sample.$batch.md.DYSF.bam" >> $pbsfile
		echo -e "samtools index -b $dir/variant_calling/$sample.$batch.md.DYSF.bam" >> $pbsfile
	fi
	if [ $panelid == 'krabbe' ]
	then
		echo -e "samtools view -h -u $dir/variant_calling/temp/$sample.$batch.md.bam chr14:87830000-87995000 | samtools sort -o $dir/variant_calling/$sample.$batch.md.GALC.bam" >> $pbsfile
		echo -e "samtools index -b $dir/variant_calling/$sample.$batch.md.GALC.bam" >> $pbsfile
	fi
	echo -e "cd $dir/variant_calling/" >> $pbsfile #move into variant_calling dir so pbs outfile goes to here
	echo -e "qsub $dir/variant_calling/splitn.pbs\n" >> $pbsfile
	echo -e "cd $dir/exoncounts/" >> $pbsfile
	echo -e "qsub $dir/exoncounts/exoncounts.pbs" >> $pbsfile
}

exon_counts () {
	local sample=$1
	local dir=$2
	if [ ! -d $dir/exoncounts ]
	then
		mkdir $dir/exoncounts
	fi
	local step=exoncounts
	local pbsfile=$dir/exoncounts/$step.pbs
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile
	echo -e "python $dexseq_count -p yes -r pos -f bam -s reverse $dexseq_gff $dir/variant_calling/temp/$sample.$batch.md.bam $dir/exoncounts/$sample.$batch.dexseqcounts.txt" >> $pbsfile
}

exoncounts_nodups () {
	local sample=$1
	local dir=$2
	local step=exoncounts.nodups
	local pbsfile=$dir/exoncounts/$step.pbs
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile
	echo -e "python $dexseq_count -p yes -r pos -f bam -s reverse $dexseq_gff $dir/splicing/$sample.remap.nodups.bam $dir/exoncounts/$sample.$batch.nodups.dexseqcounts.txt" >> $pbsfile
}

split_reads () {
	local sample=$1
	local dir=$2
	local step=splitn
	local pbsfile=$dir/variant_calling/$step.pbs
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile #start pbs file for step
	echo -e "java -jar $gatk -T SplitNCigarReads -R $wholeGenomeFasta -I $dir/variant_calling/temp/$sample.$batch.md.bam -o $dir/variant_calling/temp/$sample.$batch.md.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS" >> $pbsfile
	echo -e "cd $dir/variant_calling/" >> $pbsfile #move into variant_calling dir so pbs outfile goes to here
	echo -e "qsub $dir/variant_calling/gvcf.pbs\n" >> $pbsfile
	echo -e "qsub $dir/variant_calling/vcf.pbs\n" >> $pbsfile
}

#base_recalibration () {
#	local sample=$1
#	local dir=$2
#	local step=baserecal
#	local pbsfile=$dir/variant_calling/$step.pbs
#	echo -e "#PBS -N $step" > $pbsfile
#	echo -e "$pbshead_big" >> $pbsfile #start pbs file for step
#	echo -e "java -jar $gatk -T BaseRecalibrator -R $wholeGenomeFasta -I $dir/variant_calling/temp/$sample.dedup.split.bam -knownSites $phase1snps -knownSites $millsIndels -knownSites $dbsnp -o $dir/variant_calling/recal_data.table" >> $pbsfile
#	echo -e "java -jar $gatk -T PrintReads -R $wholeGenomeFasta -I $dir/variant_calling/temp/$sample.dedup.split.bam -BQSR $dir/variant_calling/recal_data.table -o $dir/variant_calling/temp/$sample.dedup.split.recal.bam" >> $pbsfile
#	echo -e "cd $dir/variant_calling/" >> $pbsfile #move into variant_calling dir so pbs outfile goes to here
#	echo -e "qsub $dir/variant_calling/gvcf.pbs\n" >> $pbsfile
#	echo -e "qsub $dir/variant_calling/vcf.pbs\n" >> $pbsfile
#}

gvcf () {
	local sample=$1
	local dir=$2
	local step=gvcf
	local pbsfile=$dir/variant_calling/$step.pbs
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile #start pbs file for step
	echo -e "java -jar $gatk -T HaplotypeCaller -R $wholeGenomeFasta -I $dir/variant_calling/temp/$sample.$batch.md.split.bam --dbsnp $dbsnp -kmerSize 10 -kmerSize 25 -kmerSize 35 -kmerSize 45 -dontUseSoftClippedBases -minDanglingBranchLength 1 -minPruning 3 -stand_call_conf 20.0 --emitRefConfidence GVCF -L $intervals -o $dir/variant_calling/$sample.$batch.nodups.raw.snps.indels.g.vcf" >> $pbsfile
}

vcf () {
	local sample=$1
	local dir=$2
	local step=vcf
	local pbsfile=$dir/variant_calling/$step.pbs
	echo -e "#PBS -N $step" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile #start pbs file for step
	echo -e "java -jar $gatk -T HaplotypeCaller -R $wholeGenomeFasta -I $dir/variant_calling/temp/$sample.$batch.md.split.bam --dbsnp $dbsnp -kmerSize 10 -kmerSize 25 -kmerSize 35 -kmerSize 45 -dontUseSoftClippedBases -minDanglingBranchLength 1 -minPruning 3 -stand_call_conf 20.0 -L $intervals -o $dir/variant_calling/$sample.$batch.nodups.raw.snps.indels.vcf" >> $pbsfile
	echo -e "java -jar $gatk -T VariantFiltration -R $wholeGenomeFasta -V $dir/variant_calling/$sample.$batch.nodups.raw.snps.indels.vcf -window 35 -cluster 3 -filterName FS -filter \"FS > 30.0\" -filterName QD -filter \"QD < 2.0\" -filterName DP -filter \"DP < 11.0\" -G_filterName GQ -G_filter \"GQ < 5.0\" -G_filterName HET -G_filter \"isHet == 1 && DP > 10\" -o $dir/variant_calling/$sample.$batch.nodups.filtered.vcf" >> $pbsfile
	#echo -e "perl $annovar $dir/variant_calling/$sample.wdups.filtered.vcf $tools/annovar/humandb/ -buildver hg38 -out $dir/variant_calling/$sample.wdups.filtered.annotated.vcf -remove -protocol refGene,cytoBand,genomicSuperDups,avsnp150,dbnsfp33a,clinvar_20170905,exac03,gnomad_exome -operation g,r,r,f,f,f,f,f -argument '-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs' -nastring . -vcfinput" >> $pbsfile
	#echo -e "rm $dir/variant_calling/$sample.wdups.raw.snps.indels.vcf\nrm $dir/variant_calling/$sample.wdups.filtered.annotated.vcf.avinput\nrm $dir/variant_calling/$sample.wdups.raw.snps.indels.vcf.idx" >> $pbsfile
}
###
if [ ! -d $basedir/$panelid/$batch/metrics/qorts ]
then
	mkdir $basedir/$panelid/$batch/metrics/qorts
fi
if [ ! -d $basedir/$panelid/$batch/metrics/tissuechecks ]
then
	mkdir $basedir/$panelid/$batch/metrics/tissuechecks
fi
while read sample sampledir readgroup fastqname
do
	dir=$basedir/$panelid/$batch/samples/$sampledir
	
	star_2pass $sample $dir $readgroup
	remove_dups $sample $dir
	splicing_remap $sample $dir
	readcounts_nodups $sample $dir
	norm_label_splices $sample $dir
	psi $sample $dir
	readcounts_tissuecheck $sample $dir
	qorts $sample $dir $readgroup
	qorts_remap $sample $dir $readgroup
	mark_duplicates $sample $dir
	exon_counts $sample $dir
	exoncounts_nodups $sample $dir
	split_reads $sample $dir
	#base_recalibration $sample $dir
	gvcf $sample $dir
	vcf $sample $dir
    #cd $dir/psi/
    #qsub $dir/psi/psi.pbs
    #cd $dir/splicing/
    #qsub $dir/splicing/splice_counts.pbs
    #qsub $dir/splicing/counts_nodups.pbs
	cd $dir/star_2pass/
	#qsub $dir/star_2pass/counts_tissuecheck.pbs
	#qsub $dir/star_2pass/star.pbs
done <$basedir/$panelid/$batch/samples.txt

echo "The following one liners can be used to copy the sample roster, STAR log files, STAR splice junction files, Aligned bams,  from the cluster to your own computer after pbs job has completed:"
echo -e "scp $user@iw-dm-4.pace.gatech.edu:$basedir/$panelid/$batch/samples/*/variant_calling/temp/*.dedup.ba* ./results/$panelid/bamfiles/$batch/"
echo -e "scp $user@iw-dm-4.pace.gatech.edu:$basedir/$panelid/$batch/samples/*/star_2pass/*.sj.out.tab ./results/$panelid/metrics/$batch/star_sjout/"
echo -e "scp $user@iw-dm-4.pace.gatech.edu:$basedir/$panelid/$batch/samples/*/star_2pass/*.log.final.out ./results/$panelid/metrics/$batch/star_logs/"
echo -e "scp $user@iw-dm-4.pace.gatech.edu:$basedir/$panelid/$batch/samples/*/splicing/*.nodups.sj.out.tab ./results/$panelid/metrics/$batch/star_sjout/"
echo -e "scp $user@iw-dm-4.pace.gatech.edu:$basedir/$panelid/$batch/samples/*/splicing/*.nodups.log.final.out ./results/$panelid/metrics/$batch/star_logs/"
echo -e "scp $user@iw-dm-4.pace.gatech.edu:$basedir/$panelid/$batch/samples/*/genecounts/*.raw.genecounts.txt ./results/$panelid/genecounts/$batch/"
echo -e "scp $user@iw-dm-4.pace.gatech.edu:$basedir/$panelid/$batch/samples/*/genecounts/*.genes.rpkm.gct ./results/$panelid/genecounts/$batch/"
echo -e "scp $user@iw-dm-4.pace.gatech.edu:$basedir/$panelid/$batch/samples/*/genecounts/*.nodups.raw.genecounts.txt ./results/$panelid/genecounts/$batch/"
echo -e "scp $user@iw-dm-4.pace.gatech.edu:$basedir/$panelid/$batch/samples/*/exoncounts/*.nodups.dexseqcounts.txt ./results/$panelid/exoncounts/$batch/"
echo -e "scp $user@iw-dm-4.pace.gatech.edu:$basedir/$panelid/$batch/samples/*/exoncounts/*.dexseqcounts.txt ./results/$panelid/exoncounts/$batch/"
echo -e "scp $user@iw-dm-4.pace.gatech.edu:$basedir/$panelid/$batch/metrics/tissuechecks/* ./results/$panelid/metrics/$batch/tissuechecks/"
