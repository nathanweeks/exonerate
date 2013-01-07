#!/bin/sh
# test for fastareformat utility

FASTAREFORMAT="../../src/util/fastareformat"
FASTALENGTH="../../src/util/fastalength"

UNFORMATTEDFILE="fastareformat.unformatted.test.fasta"
FORMATTEDFILE="fastareformat.formatted.test.fasta"

clean_exit(){
    rm -f $UNFORMATTEDFILE $FORMATTEDFILE
    exit $1
    }

cat << EOF > $UNFORMATTEDFILE
>test
ACGT
AACCGGTT
AAACCCGGGTTT
AAAACCCCGGGGTTTT
EOF

$FASTAREFORMAT $UNFORMATTEDFILE > $FORMATTEDFILE
if [ $? -eq 0 ]
then
    echo Reformatted OK
else
    echo Problem reformatting $UNFORMATTEDFILE
    clean_exit 1
fi

cat $FORMATTEDFILE

$FASTALENGTH $UNFORMATTEDFILE | ( read UNFORMATTEDLEN ID
    $FASTALENGTH $FORMATTEDFILE | (read FORMATTEDLEN ID
        if [ $UNFORMATTEDLEN -eq $FORMATTEDLEN ]
        then
            echo Reformatted length consistent
        else
            echo Reformatted length different
            clean_exit 1
        fi
        )
    )

clean_exit 0

