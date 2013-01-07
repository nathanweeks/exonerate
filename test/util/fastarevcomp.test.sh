#!/bin/sh
# test for fastadiff utility

FASTAREVCOMP="../../src/util/fastarevcomp"
FASTADIFF="../../src/util/fastadiff"
FASTALENGTH="../../src/util/fastalength"
CALM="../data/cdna/calm.human.dna.fasta"
CALM_RC="fastarevcomp.rc.test.fasta"
CALM_RCRC="fastarevcomp.rc.rc.test.fasta"

clean_exit(){
    rm -f $CALM_RC $CALM_RCRC
    exit $1
    }

run_revcomp(){
    $FASTAREVCOMP $1 > $2
    if [ $? -eq 0 ]
    then
        echo Reverse complement created : $1
    else
        echo Problem generating reverse complement : $1
        clean_exit 1
    fi
    }

run_revcomp $CALM $CALM_RC
run_revcomp $CALM_RC $CALM_RCRC

check_length(){
    $FASTALENGTH $1 | (read LENGTH IDENTIFIER
        if [ $LENGTH -eq 2175 ]
        then
            echo RC length is the same
        else
            echo RC length changed
            clean_exit 1
        fi
        )
    return
    }

check_length $CALM_RC
check_length $CALM_RCRC

$FASTADIFF -c no $CALM $CALM_RCRC
if [ $? -eq 0 ]
then
    echo 2nd reverse complement same as original
else
    echo 2nd reverse complement not same as original
    clean_exit 1
fi

clean_exit 0

