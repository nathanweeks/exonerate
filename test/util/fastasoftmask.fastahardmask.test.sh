#!/bin/sh
# test for fastasoftmask and fastahardmask utilities

FASTASOFTMASK="../../src/util/fastasoftmask"
FASTAHARDMASK="../../src/util/fastahardmask"
FASTADIFF="../../src/util/fastadiff"

UNMASKEDFILE="fastasoftmask.fastahardmask.unmasked.test.fasta"
MASKEDFILE="fastasoftmask.fastahardmask.masked.test.fasta"
SOFTMASKEDFILE="fastasoftmask.fastahardmask.softmasked.test.fasta"
HARDMASKEDFILE="fastasoftmask.fastahardmask.hardmasked.test.fasta"

clean_exit(){
    rm -f $UNMASKEDFILE $MASKEDFILE $SOFTMASKEDFILE $HARDMASKEDFILE
    exit $1
    }

cat << EOF > $UNMASKEDFILE
>test
ACGTAACCGGTTAAACCCGGGTTTAAAACCCCGGGGTTTT
EOF

cat << EOF > $MASKEDFILE
>test
ACGNAACCGGNNAAACCCGGGNNNAAAACCCCGGGGNNNN
EOF

$FASTASOFTMASK $UNMASKEDFILE $MASKEDFILE > $SOFTMASKEDFILE
if [ $? -eq 0 ]
then
    echo Softmasked OK
else
    echo Problem softmasking $UNMASKEDFILE
    clean_exit 1
fi

cat $SOFTMASKEDFILE

$FASTAHARDMASK $SOFTMASKEDFILE > $HARDMASKEDFILE
if [ $? -eq 0 ]
then
    echo Hardmasked OK
else
    echo Problem hardmasking $SOFTMASKEDFILE
    clean_exit 1
fi

$FASTADIFF -c no $MASKEDFILE $HARDMASKEDFILE
if [ $? -eq 0 ]
then
    echo Hardmasked file is same as original masked file
else
    echo Hardmasked file is not the same as masked
    clean_exit 1
fi

clean_exit 0

