#!/bin/sh
# test for fastaoverlap utility

PROGRAM="../../src/util/fastaoverlap"
INPUTFILE="../data/protein/calm.human.protein.fasta"
OUTPUTFILE="fastaoverlap.fasta"

clean_exit(){
    rm -f $OUTPUTFILE
    exit $1
    }

$PROGRAM $INPUTFILE --chunk 10 --jump 10 > $OUTPUTFILE
if [ $? -eq 0 ]
then
    echo $PROGRAM working OK
else
    echo Error with $PROGRAM
    clean_exit 1
fi

TOTAL_SEQS=`grep -c '^>' $OUTPUTFILE`
if [ $TOTAL_SEQS -eq 15 ]
then
    echo Generated $TOTAL_SEQS
else
    echo Unexpected number of overlap seqs: $TOTAL_SEQS
    clean_exit 1
fi

clean_exit 0

