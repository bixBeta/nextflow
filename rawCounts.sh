#!/bin/bash

if [ "$1" = "help" ] || [ -z "$1" ]
    then
    echo ""
    echo "--------------------------------------------------------------------------------------"
    echo " To run this script, use the following syntax:"
    echo "   bash" $0 "<1> or <2> or <0> or <salmon>"
    echo "    where 1 is for first strand and 2 is for reverse strands, 0 for unstranded counts"
    echo "    use 'salmon' to extract raw counts (NumReads) from Salmon quant.sf files"
    echo "    Note: for salmon, run this script from within the salmon_counts/ directory"
    echo "--------------------------------------------------------------------------------------"
    echo ""
    exit 1

elif [ "$1" = "salmon" ]
    then
        for i in *.sf
        do
        awk 'NR > 1 {printf "%s\t%d\n", $1, int($5 + 0.5)}' $i > ${i%.sf}.rawCounts
        done
        mkdir -p rawCounts
        mv *.rawCounts rawCounts

elif [ "$1" = "1" ]
    then
        for i in *ReadsPerGene.out.tab
        do
        awk 'NR > 4 {print $1 "\t" $3}' $i > $i.rawCounts
        done
        mkdir rawCounts
        mv *.rawCounts rawCounts

elif [ "$1" = "2" ]
    then
        for i in *ReadsPerGene.out.tab
        do
        awk 'NR > 4 {print $1 "\t" $4}' $i > $i.rawCounts
        done
        mkdir rawCounts
        mv *.rawCounts rawCounts

else
    for i in *ReadsPerGene.out.tab
    do
    awk 'NR > 4 {print $1 "\t" $2}' $i > $i.rawCounts
    done
    mkdir rawCounts
    mv *.rawCounts rawCounts

fi
