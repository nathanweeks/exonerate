#!/bin/sh
# test for fastasoft utility

FASTACLIP="../../src/util/fastaclip"
FASTALENGTH="../../src/util/fastalength"
UNCLIPPEDFILE="fastaclip.test.fasta"
CLIPPEDFILE="fastaclip.clipped.test.fasta"

clean_exit(){
    rm -f $UNCLIPPEDFILE $CLIPPEDFILE
    exit $1
    }

cat << EOF > $UNCLIPPEDFILE
>test
ACGTNNNN
EOF

$FASTACLIP $UNCLIPPEDFILE > $CLIPPEDFILE
if [ $? -eq 0 ]
then 
    echo Clipped OK
else
    echo Problem clipping $CLIPPEDFILE
    clean_exit 1
fi

cat $CLIPPEDFILE

$FASTALENGTH $CLIPPEDFILE | (read LENGTH IDENTIFIER
    if [ $LENGTH -eq 4 ]
    then
        echo Clipped input correctly
    else
        echo "Clipping error (length $LENGTH)"
        clean_exit 1
    fi
    )

clean_exit 0

