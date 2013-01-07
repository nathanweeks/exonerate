#!/bin/sh
# test for fastasubseq utilility

PROGRAM="../../src/util/fastasubseq"
INPUTFILE="../data/protein/calm.human.protein.fasta"

SUBSEQ=`$PROGRAM $INPUTFILE --start 10 --length 10 | tail -1`
if [ $SUBSEQ = "AEFKEAFSLF" ]
then
    echo Generated correct subseq
else
    echo Wrong subseq generated $SUBSEQ
    exit 1
fi

exit $?

