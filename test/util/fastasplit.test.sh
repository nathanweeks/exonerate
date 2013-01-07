#!/bin/sh
# test for fastasplit utility

FASTASPLIT="../../src/util/fastasplit"
FASTALENGTH="../../src/util/fastalength"

INPUTDIR="../data/protein"
INPUTFILE="fastasplit.test.fasta"
SPLITDIR="fastasplit.test.dir"

OUTPUTFILE_0="$SPLITDIR/fastasplit.test.fasta_chunk_0000000"
OUTPUTFILE_1="$SPLITDIR/fastasplit.test.fasta_chunk_0000001"

clean_exit(){
    rm -f $INPUTFILE
    rm -rf $SPLITDIR
    exit $1
    }

cat $INPUTDIR/* > $INPUTFILE

rm -rf $SPLITDIR
mkdir $SPLITDIR

$FASTASPLIT $INPUTFILE $SPLITDIR --chunk 2
if [ $? -eq 0 ]
then
    echo Split $INPUTFILE OK
else
    echo Problem splitting $INPUTFILE
    clean_exit 1
fi

TOTALFILES=`ls -1 $SPLITDIR | wc -l`
if [ $TOTALFILES -eq 2 ]
then
    echo Split into two files as expected
else
    echo Did not expect $TOTALFILES files
    clean_exit 1
fi

total_seq_len(){
    TOTAL_LEN=0
    for LEN in `$FASTALENGTH $1 | cut -d' ' -f1`
    do
        TOTAL_LEN=`expr $TOTAL_LEN + $LEN`
    done
    echo $TOTAL_LEN
    return
    }

# check total lengths
INPUT_LEN=`total_seq_len $INPUTFILE`
echo Input len $INPUT_LEN
OUTPUT_LEN_0=`total_seq_len $OUTPUTFILE_0`
OUTPUT_LEN_1=`total_seq_len $OUTPUTFILE_1`
OUTPUT_LEN=`expr $OUTPUT_LEN_0 + $OUTPUT_LEN_1`
echo Output len $OUTPUT_LEN
if [ $INPUT_LEN -eq $OUTPUT_LEN ]
then
    echo Input and output lengths match
else
    echo Output length $OUTPUT_LEN length, but expected $INPUT_LEN
    clean_exit 1
fi

clean_exit 0

