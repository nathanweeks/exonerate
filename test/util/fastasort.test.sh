#!/bin/sh
# test for fastasoft utility

FASTASORT="../../src/util/fastasort"
FASTALENGTH="../../src/util/fastalength"
INPUTDIR="../data/protein"
INPUTFILE="fastasort.test.fasta"
OUTPUTFILE="fastasort.sorted.test.fasta"

cat $INPUTDIR/*.fasta > $INPUTFILE

$FASTASORT $INPUTFILE --key len > $OUTPUTFILE

PREV_ID="NONE"
PREV_LEN=0

clean_exit(){
    rm -f $INPUTFILE $INDEXFILE $OUTPUTFILE
    exit $1
    }

$FASTALENGTH $OUTPUTFILE | while read LENGTH IDENTIFIER
do
    if [ $LENGTH -ge $PREV_LEN ]
    then
       echo Sorted $IDENTIFIER $LENGTH
    else
       echo Sorting failure: $PREV_ID $PREV_LEN : $IDENTIFIER $LENGTH
       clean_exit 1
    fi
    PREV_ID=$IDENTIFIER
    PREV_LEN=$LENGTH
done

clean_exit 0

