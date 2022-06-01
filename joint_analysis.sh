#!/bin/bash
scriptname=$0   # $0 is the name of the program
#make list of tools/paths/versions to echo onto a roster of files for each sample? Also list of required files? & modules?
panel=ibdpid #either ibdpid or nmd or krabbe
group=all #group id from sampleroster.txt that can subset the analysis to fewer samples. for group: either blood or muscle - other group types will attempt to run but will produce a warning
vars=~/p-ggibson3-0/rich_project_pb1/scripts/analysis_vars.sh
###
#help function
HELP () {
	echo -e "\n\tThis script will execute joint analysis scripts on samples from the sample roster to produce an annotated joint vcf, a matrix of non-duplicate splice counts, and a matrix of PSI calculations.\n\n"
	echo -e "\n\tUSAGE: $scriptname [-g GROUP] [-p (ibdpid|nmd|krabbe)] [-v VARIABLES_FILE]\n"
	echo -e "\tOPTIONS:\n"
	echo -e "\t\t-a GROUP\t\tgroup ID of samples to perform joint analyses for, default is all samples\n\n"
	echo -e "\t\t-v VARIABLES_FILE\tLocation of file containing variable information for analysis scripts"
	echo -e "\t\t\t\t\tDEFAULT: ~/p-ggibson3-0/rich_project_pb1/scripts/analysis_vars.sh\n"	
	echo -e "\t\t-p PANEL\t\tPanel used for Targeted RNA sequencing"  
	echo -e "\t\t\t\t\tDEFAULT: ibdpid\n"						
	exit 0
}
#Check the number of arguments. If none are passed, print help and exit.
NUMARGS=$#
if [ $NUMARGS -eq 0 ]; then
  HELP
fi
#Get Options from command line
while getopts :g:p:v:h opt
do
	case $opt in
	g)group=$OPTARG;;
	p)panel=$OPTARG;;
	v)vars=$OPTARG;;
    h) HELP; exit;;
    \?) echo "ERROR: Invalid option: -$OPTARG" >&2; echo -e "Run $scriptname -h to see usage and help" >&2; exit 1;;
    :) echo "ERROR: Option -$OPTARG requires an argument" >&2; echo -e "Run $scriptname -h to see usage and help" >&2; exit 1;;
    *) echo "ERROR: Unimplemented option: -$OPTARG" >&2; echo -e "Run $scriptname -h to see usage and help" >&2; exit 1;;
  esac
done
shift $(($OPTIND -1))
#check to make sure variables file exists
if [ ! -f $vars ] 
then
    echo -e "ERROR: File $vars DOES NOT exist" >&2
	echo -e "Edit this scripts default or use the -v option to point the script to your file that contains variables for analysis scripts" >&2
    exit 1
fi
###
source $vars
###
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
#ensure samples.txt exists
if [ ! -f $basedir/$panel/sampleroster.txt ] 
then
    echo -e "ERROR: File $basedir/$panel/sampleroster.txt DOES NOT exist" >&2
	echo -e "Create the sample roster for $panel containing these five columns for each sample: ID, group, batch, directory, id2" >&2
    exit 1
fi

if [ ! -d $basedir/$panel/joint_analysis ] 
then
    echo -e "Directory $basedir/$panel/joint_analysis DOES NOT exist - creating directory" >&2
	mkdir $basedir/$panel/joint_analysis
fi
###
if [ $group == 'all' ]
then
	echo -e "All samples will be analyzed jointly:\n"
	cat $basedir/$panel/sampleroster.txt
else
	echo -e "Samples that match group ID \"$group\":\n"
	awk -v var="$group" '$2 ~ var' $basedir/$panel/sampleroster.txt
fi

###
variant_call () {
	local pbsfile=$basedir/$panel/joint_analysis/$group.jointvarcall.pbs
	incl=""
	excl=""
	echo -e "#PBS -N $group.jointvarcall" > $pbsfile
	echo -e "$pbshead_big" >> $pbsfile #start pbs file for step
	#prepare gvcf list based on group type
	if [ -f $basedir/$panel/joint_analysis/$group.gvcfs.list ]
	then
		rm $basedir/$panel/joint_analysis/$group.gvcfs.list
	fi
	if [ -f $basedir/$panel/joint_analysis/$group.$panel.filtered.annotated.processed.txt ]
	then
		rm $basedir/$panel/joint_analysis/$group.$panel.filtered.annotated.processed.txt
	fi
	while read sample gtype batch sampledir sid
	do
		if [ $group == 'all' ]
		then
			if [ ! -f $basedir/$panel/$batch/samples/$sampledir/variant_calling/$sample.$batch.nodups.raw.snps.indels.g.vcf ]
			then
				excl="${excl}\n${sample}"
			else
				incl="${incl}\n${sample}"
				echo -e "$basedir/$panel/$batch/samples/$sampledir/variant_calling/$sample.$batch.nodups.raw.snps.indels.g.vcf" >> $basedir/$panel/joint_analysis/$group.gvcfs.list
			fi
		else
			if [[ $gtype =~ $group ]]
			then
				if [ ! -f $basedir/$panel/$batch/samples/$sampledir/variant_calling/$sample.$batch.nodups.raw.snps.indels.g.vcf ]
				then
					excl="${excl}\n${sample}"
				else
					incl="${incl}\n${sample}"
					echo -e "$basedir/$panel/$batch/samples/$sampledir/variant_calling/$sample.$batch.nodups.raw.snps.indels.g.vcf" >> $basedir/$panel/joint_analysis/$group.gvcfs.list
				fi
			fi
		fi
	done <$basedir/$panel/sampleroster.txt
	echo -e "Samples that will be included in joint variant call:\n$incl\n"
	echo -e "Samples that do not have a gvcf file and will be excluded from joint variant call:\n$excl\n"
	#create joint variant call pbs file
	echo -e "java -jar $gatk -T GenotypeGVCFs -nt 2 -L $intervals -R $wholeGenomeFasta -stand_call_conf 20.0 -V $basedir/$panel/joint_analysis/$group.gvcfs.list -o $basedir/$panel/joint_analysis/$group.$panel.raw.snps.indels.vcf" >> $pbsfile
	echo -e "java -jar $gatk -T VariantFiltration -R $wholeGenomeFasta -V $basedir/$panel/joint_analysis/$group.$panel.raw.snps.indels.vcf -window 35 -cluster 3 -filterName FS -filter \"FS > 30.0\" -filterName QD -filter \"QD < 2.0\" -G_filterName DP -G_filter \"DP < 10\" -G_filterName GQ -G_filter \"GQ < 5.0\" -G_filterName HET -G_filter \"isHet == 1 && DP > 9\" -G_filterName REF -G_filter \"isHomRef == 1 && DP > 9\" -G_filterName HOM -G_filter \"isHomVar == 1 && DP > 9\" -o $basedir/$panel/joint_analysis/$group.$panel.filtered.vcf" >> $pbsfile
	echo -e "perl $annovar $basedir/$panel/joint_analysis/$group.$panel.filtered.vcf $tools/annovar/humandb/ -buildver hg38 -out $basedir/$panel/joint_analysis/$group.$panel.filtered.annotated.vcf -remove -protocol refGene,cytoBand,genomicSuperDups,avsnp150,dbnsfp35a,clinvar_20190305,exac03,gnomad211_exome -operation g,r,r,f,f,f,f,f -argument '-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs' -intronhgvs 50 -nastring . -otherinfo -vcfinput" >> $pbsfile
	echo -e "$process_vars $basedir/$panel/joint_analysis/$group.$panel.filtered.annotated.vcf.hg38_multianno.txt $trapscores $basedir/$panel/joint_analysis/$group.$panel.filtered.annotated.processed.txt" >> $pbsfile
	cd $basedir/$panel/joint_analysis/
	qsub $pbsfile
}
###
splicing () {
	local pbsfile=$basedir/$panel/joint_analysis/$group.jointsplicing.pbs
	incl=""
	excl=""
	echo -e "#PBS -N $group.jointsplicing" > $pbsfile
	echo -e "$pbshead_little" >> $pbsfile #start pbs file for step
	#prepare splicing file list based on group type
	#remove old files from joint analysis

	if [ -f $basedir/$panel/joint_analysis/$group.splicefiles.txt ]
	then
		rm $basedir/$panel/joint_analysis/$group.splicefiles.txt
	fi
	if [ -f $basedir/$panel/joint_analysis/$group.sjs.txt ]
	then
		rm $basedir/$panel/joint_analysis/$group.sjs.txt
	fi
	if [ -f $basedir/$panel/joint_analysis/$group.sjs.uniq.txt ]
	then
		rm $basedir/$panel/joint_analysis/$group.sjs.uniq.txt
	fi
	if [ -f $basedir/$panel/joint_analysis/$group.sampleroster.subset.txt ]
	then
		rm $basedir/$panel/joint_analysis/$group.sampleroster.subset.txt
	fi
	if [ -f $basedir/$panel/joint_analysis/$group.splices.sorted.txt ]
	then
		rm $basedir/$panel/joint_analysis/$group.splices.sorted.txt
	fi
	while read sample ttype batch sampledir sid
	do
		if [ $group == 'all' ]
		then
			if [ ! -f $basedir/$panel/$batch/samples/$sampledir/splicing/$sample.nodups.raw.splices.txt ]
			then
				excl="${excl}\n${sample}"
			else
				incl="${incl}\n${sample}"
				echo -e "$sample\t$ttype\t$batch\t$sampledir" >> $basedir/$panel/joint_analysis/$group.sampleroster.subset.txt
				echo -e "$sample\t$basedir/$panel/$batch/samples/$sampledir" >> $basedir/$panel/joint_analysis/$group.splicefiles.txt
				awk '{OFS="\t"}{print $1, $2, $3, $4, $5, $6}' $basedir/$panel/$batch/samples/$sampledir/splicing/$sample.nodups.raw.splices.txt >> $basedir/$panel/joint_analysis/$group.sjs.txt
			fi
		else
			if [[ $ttype =~ $group ]]
			then
				if [ ! -f $basedir/$panel/$batch/samples/$sampledir/splicing/$sample.nodups.raw.splices.txt ]
				then
					excl="${excl}\n${sample}"
				else
					incl="${incl}\n${sample}"
					echo -e "$sample\t$ttype\t$batch\t$sampledir" >> $basedir/$panel/joint_analysis/$group.sampleroster.subset.txt
					echo -e "$sample\t${sample}.${batch}\t$basedir/$panel/$batch/samples/$sampledir" >> $basedir/$panel/joint_analysis/$group.splicefiles.txt
					awk '{OFS="\t"}{print $1, $2, $3, $4, $5, $6}' $basedir/$panel/$batch/samples/$sampledir/splicing/$sample.nodups.raw.splices.txt >> $basedir/$panel/joint_analysis/$group.sjs.txt
				fi
			fi
		fi
	done <$basedir/$panel/sampleroster.txt
	echo -e "Samples that will be included in joint splicing:\n$incl\n"
	echo -e "Samples that do not have a splicing file and will be excluded from joint splicing:\n$excl\n"
	#sort by chr, then start, then stop, then gene (for A-T genes, order will be correct, beyond that the "unknown" will come before the actual gene which is backwards)
	sort -k4V,4 -k5n,5 -k6n,6 -k1,1 $basedir/$panel/joint_analysis/$group.sjs.txt | uniq > $basedir/$panel/joint_analysis/$group.sjs.uniq.txt
	echo -e "$joint_splicing $basedir/$panel/joint_analysis $group" >> $pbsfile
	echo -e "sort -t \":\" -k4V,4 -k5n,5 -k6n,6 $basedir/$panel/joint_analysis/$group.splices.txt >> $basedir/$panel/joint_analysis/$group.splices.sorted.txt" >> $pbsfile
	#echo -e "Rscript $zscores_script -sampleinfo $basedir/$panel/joint_analysis/$group.sampleroster.subset.txt -panel $panel -splices $basedir/$panel/joint_analysis/$group.splices.sorted.txt" >> $pbsfile
	cd $basedir/$panel/joint_analysis/
	qsub $pbsfile
}
###
psi () {
	incl=""
	excl=""
	counter=1
	first=""
	second=""
	others=""
	while read sample ttype batch sampledir sid
	do
		if [ $group == 'all' ]
		then
			if [ ! -f $basedir/$panel/$batch/samples/$sampledir/psi/$sample.$batch.exonic_parts.psi ]
			then
				excl="${excl}\n${sample}"
			else
				incl="${incl}\n${sample}"
				if [ "$counter" -eq 1 ]
				then
					first=$basedir/$panel/$batch/samples/$sampledir/psi/$sample.$batch.exonic_parts.psi
					counter=2
					continue
				fi
				if [ "$counter" -eq 2 ]
				then
					join $first $basedir/$panel/$batch/samples/$sampledir/psi/$sample.$batch.exonic_parts.psi > $basedir/$panel/joint_analysis/$group.$panel.psi.txt
					#second=$basedir/$panel/$batch/samples/$sampledir/psi/$sample.exonic_parts.psi
					counter=3
					continue
				fi
				if [ "$counter" -eq 3 ]
				then
					join $basedir/$panel/joint_analysis/$group.$panel.psi.txt $basedir/$panel/$batch/samples/$sampledir/psi/$sample.$batch.exonic_parts.psi > $basedir/$panel/joint_analysis/$group.$panel.temp.psi.txt
					mv $basedir/$panel/joint_analysis/$group.$panel.temp.psi.txt $basedir/$panel/joint_analysis/$group.$panel.psi.txt
					continue
				fi
			fi
		else
			if [[ $ttype =~ $group ]]
			then
				if [ ! -f $basedir/$panel/$batch/samples/$sampledir/psi/$sample.$batch.exonic_parts.psi ]
				then
					excl="${excl}\n${sample}"
				else
					incl="${incl}\n${sample}"
					#sed -i "1i #\t#\t$sample\t$sample\t$sample" $basedir/$panel/$batch/samples/$sampledir/psi/$sample.exonic_parts.psi
					if [ "$counter" -eq 1 ]
					then
						first=$basedir/$panel/$batch/samples/$sampledir/psi/$sample.$batch.exonic_parts.psi
						counter=2
						continue
					fi
					if [ "$counter" -eq 2 ]
					then
						join $first $basedir/$panel/$batch/samples/$sampledir/psi/$sample.$batch.exonic_parts.psi > $basedir/$panel/joint_analysis/$group.$panel.psi.txt
						counter=3
						continue
					fi
					if [ "$counter" -eq 3 ]
					then
						join $basedir/$panel/joint_analysis/$group.$panel.psi.txt $basedir/$panel/$batch/samples/$sampledir/psi/$sample.$batch.exonic_parts.psi > $basedir/$panel/joint_analysis/$group.$panel.temp.psi.txt
						mv $basedir/$panel/joint_analysis/$group.$panel.temp.psi.txt $basedir/$panel/joint_analysis/$group.$panel.psi.txt
						continue
					fi
				fi
			fi
		fi
	done <$basedir/$panel/sampleroster.txt
	echo -e "Samples that were included in psi matrix:\n$incl\n"
	echo -e "Samples that do not have a psi file and were excluded from psi matrix:\n$excl\n"
}
###
variant_call
splicing
psi
	