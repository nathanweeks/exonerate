#!/bin/sh
# test for fastadiff utility

FASTADIFF="../../src/util/fastadiff"
CALM="../data/protein/calm.human.protein.fasta"
P53="../data/protein/p53.human.protein.fasta"

$FASTADIFF $CALM $P53
if [ $? -eq 0 ]
then
    echo Diff test on different files failed
    exit 1
else
    echo Different seqs recognised as different
fi

$FASTADIFF $CALM $CALM
if [ $? -eq 0 ]
then
    echo Identity test OK
else
    echo Identity failed
    exit 1
fi

exit 0

