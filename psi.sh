#!/bin/bash
scriptname=$0   # $0 is the name of the program
sampledir=""
sample=""
readLength=""
vars=""
panelvars=""
batch=""
###
#help function
HELP () {
	echo -e "\n\tUSAGE: $scriptname -d DIR -s SAMPLEID -b BATCH -l READ_LENGTH -v GENERAL_VARIABLES -p PANEL_VARIABLES\n"
	echo -e "\tMANDATORY OPTIONS:\n"
	echo -e "\t\t-d DIR\t\tAbsolute Path to the sample directory\n"
	echo -e "\t\t-s SAMPLEID\tSample ID used in file names\n\n"
	echo -e "\t\t-v GENERAL_VARIABLES\tLocation of file containing variable information for analysis scripts"
	echo -e "\t\t-p PANEL_VARIABLES\tLocation of file containing panel specific variable information for analysis scripts"
	echo -e "\t\t-l READ_LENGTH\tRead Length for sample reads. Must be a positive integer."
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
while getopts :s:b:d:l:v:p:h opt
do
	case $opt in
	s)sample=$OPTARG;;
	b)batch=$OPTARG;;
	d)sampledir=$OPTARG;;
	l)readLength=$OPTARG;;
	v)vars=$OPTARG;;
	p)panelvars=$OPTARG;;
    h) HELP; exit;;
    \?) echo "ERROR: Invalid option: -$OPTARG" >&2; echo -e "Run $scriptname -h to see usage and help" >&2; exit 1;;
    :) echo "ERROR: Option -$OPTARG requires an argument" >&2; echo -e "Run $scriptname -h to see usage and help" >&2; exit 1;;
    *) echo "ERROR: Unimplemented option: -$OPTARG" >&2; echo -e "Run $scriptname -h to see usage and help" >&2; exit 1;;
  esac
done
shift $(($OPTIND -1))
###
#check to make sure variables files exist
if [ ! -f $vars ] 
then
    echo -e "ERROR: File $vars DOES NOT exist" >&2
	echo -e "Use the -v option to point the script to your file that contains variables for analysis scripts" >&2
    exit 1
fi
source $vars
###
if [ ! -f $panelvars ] 
then
    echo -e "ERROR: File $panelvars DOES NOT exist" >&2
	echo -e "Use the -v option to point the script to your file that contains panel specific variables for analysis scripts" >&2
    exit 1
fi
source $panelvars
#make sure read length is a positive integer
if ! [[ "$readLength" =~ ^[0-9]+$ ]]
    then
        echo -e "ERROR: \"$readLength\" is not a valid argument for [-l] option. Read Length must be an integer."
		exit 1
fi
###
#put sj.out.tab into correct bed format for use
awk 'BEGIN{OFS="\t"}{print $1, $2-20-1, $3+20, "JUNCBJ"NR, $7,
($4 == 1)? "+":"-",$2-20-1, $3+20, "255,0,0", 2, "20,20",
"0,300" }' $sampledir/star_2pass/$sample.$batch.sj.out.tab > $sampledir/psi/$sample.junctions.bed
###
#get IR count
$bedtools/coverageBed -split -s -sorted -g $wholeGenomeFasta.fai -a $exonparts -b $sampledir/star_2pass/$sample.bam | awk 'BEGIN{OFS="\t"}{print $1,$4,$5,$5-$4+1,$9,$10}' | sort -k 5 > $sampledir/psi/$sample.exonic_parts.inclusion
###
#define intronic locations w/in junction file
sed 's/,/\t/g' $sampledir/psi/$sample.junctions.bed | grep -v description | awk '{OFS="\t"}{print $1,$2+$13,$3-$14,$4,$5,$6}' > $sampledir/psi/$sample.intron.bed
rm $sampledir/psi/$sample.junctions.bed
###
#count ER by intersecting introns with exonic parts and counting number of split alignments overlapping each exonic part
$bedtools/intersectBed -wao -f 1.0 -s -sorted -g $wholeGenomeFasta.fai -a $exonparts -b $sampledir/psi/$sample.intron.bed | awk 'BEGIN{OFS="\t"}{$16 == 0? s[$9] += 0:s[$9] += $14}END{for (i in s) {print i,s[i]}}' | sort -k 1 > $sampledir/psi/$sample.exonic_parts.exclusion
rm $sampledir/psi/$sample.intron.bed
###
#psi calculation
paste $sampledir/psi/$sample.exonic_parts.inclusion $sampledir/psi/$sample.exonic_parts.exclusion | awk -v "len=$readLength" 'BEGIN{OFS="\t"; print "exon_ID" , "length" , "inclusion" , "exclusion" , "PSI"}{NIR=$6/($4+len-1) ; NER=$8/(len-1)}{print $5,$4,$6,$8,(NIR+NER<=0)? "NA":NIR / (NIR + NER)}' > $sampledir/psi/$sample.$batch.exonic_parts.psi
rm $sampledir/psi/$sample.exonic_parts.inclusion
rm $sampledir/psi/$sample.exonic_parts.exclusion
