#!/bin/sh
# test for fastareformat utility

FASTAVALIDCDS="../../src/util/fastavalidcds"

INPUTFILE="fastavalidcds.input.test.fasta"
OUTPUTFILE="fastavalidcds.output.test.fasta"

clean_exit(){
    rm -f $INPUTFILE $OUTPUTFILE
    exit $1
    }

cat << EOF > $INPUTFILE
>odd_length
ATGATAA
>too_short
ATAA
>no_start
AAACCCGGGTGA
>no_end
ATGCCCGGGTTT
>in_frame_stop
ATGAAACCCTAGGGGTTTTAG
>non_acgt_base
ATGAANTGA
>valid_seq
ATGAAACCCGGGTTTTAA
EOF

$FASTAVALIDCDS -e yes $INPUTFILE > $OUTPUTFILE
if [ $? -eq 0 ]
then
    echo Validiated CDS sequences OK
else
    echo Problem validating CDS in $INPUTFILE
    clean_exit 1
fi

cat $OUTPUTFILE

TOTAL=`grep -c '^>' $OUTPUTFILE`
if [ $TOTAL -eq 1 ]
then
    echo Validated correctly a single seq from input
else
    echo "Missed bad CDS $INPUTFILE (found: $TOTAL)"
    clean_exit 1
fi

clean_exit 0

