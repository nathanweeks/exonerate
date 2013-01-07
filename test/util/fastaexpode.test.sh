#!/bin/sh
# test for fastaexpode utility

FASTAEXPLODE="../../src/util/fastaexplode"
INPUTDIR="../data/protein"
INPUTFILE="fastaexplode.test.fasta"
OUTPUTDIR="fastaexplode.out_dir.test.fasta"

clean_exit(){
    rm -f $INPUTFILE
    rm -rf $OUTPUTDIR
    exit $1
    }

cat ${INPUTDIR}/*.fasta > $INPUTFILE
if [ $? -eq 0 ]
then
    echo Made input file for fastaexplode
else
    echo Problem making input file for fastaexplode
    clean_exit 1
fi

mkdir $OUTPUTDIR
if [ $? -eq 0 ]
then
    echo Make output directory for fastaexplode
else
    echo Problem making output directory for fastaexplode
    clean_exit 1
fi

$FASTAEXPLODE $INPUTFILE --directory $OUTPUTDIR
if [ $? -eq 0 ]
then
    echo Successfully ran fastaexplode on: $INPUTFILE
else
    echo Problem running fastaexplode on: $INPUTFILE
    clean_exit 1
fi

TOTALINPUT=`grep -c '^>' $INPUTFILE`
TOTALOUTPUT=`ls -1 ${OUTPUTDIR}/*.fa | wc -l`
if [ $TOTALINPUT -eq $TOTALOUTPUT ]
then
    echo Expected number of seqs in output: $TOTALINPUT
else
    echo Mismathed number of seqs in output $TOTALINPUT,$TOTALOUTPUT
    clean_exit 1
fi

clean_exit 0

