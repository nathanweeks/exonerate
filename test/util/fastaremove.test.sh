#!/bin/sh
# test for fastaremove utilility

FASTAREMOVE="../../src/util/fastaremove"
INPUTDIR="../data/protein/"
INPUTFILE="fastaremove.test.in.fasta"
OUTPUTFILE="fastaremove.test.out.fasta"

cat $INPUTDIR/*.fasta > $INPUTFILE

clean_exit(){
    rm -f $INPUTFILE $OUTPUTFILE
    exit $1
    }

HAVE_CALM_HUMAN=`grep -c '^>CALM_HUMAN' $INPUTFILE`
if [ $HAVE_CALM_HUMAN -eq 1 ]
then
    echo Have calmodulin in input
else
    echo Calmodulin missing from input
    clean_exit 1
fi

$FASTAREMOVE $INPUTFILE CALM_HUMAN > $OUTPUTFILE

HAVE_CALM_HUMAN=`grep -c '^>CALM_HUMAN' $OUTPUTFILE`
if [ $HAVE_CALM_HUMAN -eq 0 ]
then
    echo Calmodulin successfully removed
else
    echo Failed to remove calmodulin from input $HAVE_CALM_HUMAN
    clean_exit 1
fi

clean_exit 0

