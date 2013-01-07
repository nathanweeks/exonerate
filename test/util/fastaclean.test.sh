#!/bin/sh
# test for fastaclean utility

PROGRAM="../../src/util/fastaclean"
INPUTFILE="../data/protein/calm.human.protein.fasta"

$PROGRAM -p TRUE $INPUTFILE

exit $?

