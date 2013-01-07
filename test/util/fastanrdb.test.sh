#!/bin/sh
# test for fastanrdb utility

FASTANRDB="../../src/util/fastanrdb"
INPUTDIR="../data/cdna"
INPUTFILE="fastanrdb.input.test.fasta"
OUTPUTFILE="fastanrdb.output.test.fasta"

clean_exit(){
    rm -f $INPUTFILE $OUTPUTFILE
    exit $1
    }

cat ${INPUTDIR}/*.fasta ${INPUTDIR}/*.fasta > $INPUTFILE
if [ $? -eq 0 ]
then
    echo Made input file for fastanrdb
else
    echo Problem making input file for fastanrdb
    clean_exit 1
fi

$FASTANRDB $INPUTFILE > $OUTPUTFILE
if [ $? -eq 0 ]
then
    echo Successfully ran fastanrdb on: $INPUTFILE
else
    echo Problem running fastanrdb on: $INPUTFILE
    clean_exit 1
fi

TOTALINPUT=`ls -1 ${INPUTDIR}/*.fasta | wc -l`
TOTALOUTPUT=`grep -c '^>' $OUTPUTFILE`
if [ $TOTALINPUT -eq $TOTALOUTPUT ]
then
    echo Expected number of seqs in output: $TOTALINPUT
else
    echo Mismatched number of seqs in output $TOTALINPUT,$TOTALOUTPUT
    clean_exit 1
fi

clean_exit 0

