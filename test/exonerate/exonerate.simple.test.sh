#!/bin/sh
# simple test for exonerate

EXONERATE="../../src/program/exonerate"
SEQUENCEFILE="../data/cdna/calm.human.dna.fasta"

OUTPUTFILE="exonerate.simple.test.out"

clean_exit(){
    rm -f $OUTPUTFILE
    exit $?
    }

$EXONERATE --bestn 1 --showvulgar yes \
           $SEQUENCEFILE $SEQUENCEFILE > $OUTPUTFILE
if [ $? -eq 0 ]
then
    echo "Exonerate simple test OK"
else
    echo "Problem running simple test for exonerate"
    clean_exit 1
fi

SCORE=`grep '^vulgar:' $OUTPUTFILE | cut -d' ' -f10`
if [ $SCORE -eq 10875 ]
then
    echo Score as expected: $SCORE
else
    echo Unexpected score: $SCORE
    clean_exit 1
fi

clean_exit 0

