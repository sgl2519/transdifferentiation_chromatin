#!/bin/bash

#********
# USAGE *
#********

display_usage() { 
	echo -e "\nDESCRIPTION: provided a list of bw files and a .bed file of regions of interest, it returns an aggregation plot for each bw\n"
	echo -e "\t--bw <path_to_bw_files> (i.e. a .txt file with the path to the bw files to be used by bwtool)\n"
	echo -e "\t--bedfile <bedFile of regions of interest>\n"
	echo -e "\t--type <gene_body/tss/enhancer>\n"
	echo -e "\t--meta <meta value for bwtool aggregate> (required only if 'type'=='gene_body')\n"
        echo -e "\t--outFile <output file> (default: aggregation.plot.tsv)\n"
	echo -e "\t--keep <TRUE/FALSE> (default: FALSE)\n"
	echo -e "\t--verbose <TRUE/FALSE> (default: FALSE)\n"
} 


if [[  $1 == "--help" ||  $1 == "-h" ]]
then
    	display_usage
        exit 0
fi

if [  $# -le 3  ]
then
	echo -e "ERROR: insufficient number of arguments\n"
    	display_usage
        exit 1
fi

#******************
# READING OPTIONS *
#******************

while [[ $# -gt 1 ]]; do

	key="$1"
	
	case $key in

	--bw)
	bw="$2"
	shift # past argument
	;;
    	
	--bedfile)
	bedFile="$2"
	shift # past argument
	;;

	--type)
	type="$2"
	shift # past argument
	;;
	
	--meta)
	meta="$2"
	shift
	;;

	--outFile)
	outFile="$2"
	shift
	;;

	--keep)
	keep="$2"
	shift
	;;

	--verbose)
	verbose="$2"
	;;
	*)
	
	;;
	esac
	shift
done


: ${outFile:="aggregation.plot.tsv"}
: ${keep:="FALSE"}
: ${verbose:="FALSE"}


if [[ "$type" != "gene_body" && "$type" != "tss" && "$type" != "enhancer" ]]
then
	echo "ERROR: 'type' must be either 'gene_body' or 'tss' or 'enhancer'"
	exit 1
fi 

re='^[0-9]+$'
if [[ "$type" == "gene_body" ]] && ! [[ $meta =~ $re ]]
then
	echo "ERROR: 'meta' option is required when 'type'=='gene_body' and must be a number"
	exit 1
fi

if [[ "$verbose" == "TRUE" ]]
then
	echo -e "\nReading options .."
	echo -e "\tpath to bw files =" "${bw}"
	echo -e "\tbedFile =" "${bedFile}"
	echo -e "\ttype =" "${type}"
	echo -e "\tmeta option =" "${meta}"
	echo -e "\toutput file =" "${outFile}"
	echo -e "\tkeep tmp files: " "${keep}\n"
fi

#********
# BEGIN *
#********

if [[ "$verbose" == "TRUE" ]]
then
	echo -e "Generating aggregation table ..\n"
fi

cat $bw | while read file; do
	sampleId="$(basename $file)"
	sampleId="${sampleId%.*}"
	
	if [[ "$type" == "tss" ]]
	then	
		bwtool aggregate 5000:5000 -starts -keep-bed $bedFile $file "$sampleId".bwtool.aggregate.tsv

	elif [[ "$type" == "enhancer" ]]
	then
		bwtool aggregate 10000:0 -starts -keep-bed $bedFile $file "$sampleId".bwtool.aggregate.tsv

	else
		bwtool aggregate 5000:$meta:5000 -keep-bed $bedFile $file "$sampleId".bwtool.aggregate.tsv

	fi
	
	awk -v sampleId="$sampleId" 'BEGIN{OFS=FS="\t"}{
	print $0, sampleId
	}' "$sampleId".bwtool.aggregate.tsv  		
		
	if [[ "$keep" == "FALSE" ]]
	then
                   rm "$sampleId".bwtool.aggregate.tsv
	fi
        
done > "$outFile"

if [[ "$type" == "gene_body" ]]
then 	
	grep -v "#" "$outFile" > tmp
	mv tmp "$outFile"
fi

# remove duplicated lines

if [[ "$verbose" == "TRUE" ]]
then
	echo -e "Removing duplicated lines ..\n"
fi

sort -u "$outFile" | sort -k3,3 -k1,1g > tmp
mv tmp "$outFile" 



if [[ "$verbose" == "TRUE" ]]
then
	echo "Headerizing .."
fi

sed -i '1iposition\tvalue\tsample' "$outFile"
