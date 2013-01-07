#!/bin/sh
# test for fastaindex and fastafetch utilities

FASTAINDEX="../../src/util/fastaindex"
FASTAFETCH="../../src/util/fastafetch"
INPUTDIR="../data/protein"
INPUTFILE="fastafetch.test.fasta"
INDEXFILE="fastafetch.test.idx"

rm -f $INDEXFILE

cat $INPUTDIR/*.fasta > $INPUTFILE

# test fastaindex
$FASTAINDEX $INPUTFILE $INDEXFILE
if [ $? -ne 0 ]
then
    echo Error with $FASTAINDEX
    exit $?
else
    echo Made index $INDEXFILE
fi

clean_exit(){
    rm -f $INPUTFILE $INDEXFILE
    exit $1
    }

run_fastafetch(){
   IDENTIFIER=$1
   ANSWER=$2
   echo $FASTAFETCH $INPUTFILE $INDEXFILE $IDENTIFIER
   echo expect $ANSWER
   $FASTAFETCH $INPUTFILE $INDEXFILE $IDENTIFIER
   RESULT=$?
   if [ $RESULT -ne $ANSWER ]
    then
        echo Error running $FASTAINDEX $INPUTFILE $INDEXFILE
        clean_exit 1
    fi
    return $RESULT
    }

# Prevent fatal debugging errors
unset EXONERATE_DEBUG

# look for IDs which are present
run_fastafetch CALM_HUMAN 0
run_fastafetch P53_HUMAN 0

# look for IDs which are missing
run_fastafetch A_MISSING_FROM_START 1
run_fastafetch M_MISSING_FROM_MIDDLE 1
run_fastafetch z_MISSING_FROM_END 1

echo fastaindex fastafetch test OK
clean_exit 0

