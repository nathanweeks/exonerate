#!/bin/sh
# simple test for ipcress

IPCRESS="../../src/program/ipcress"
SEQUENCEFILE="../data/cdna/calm.human.dna.fasta"

IPCRESSFILE="test.ipcress"
OUTPUTFILE="test.ipcress.out"

clean_exit(){
    rm -f $IPCRESSFILE $OUTPUTFILE
    exit $?
    }

cat > $IPCRESSFILE << EOF
test_primer CGCGGACGCGCG GTATTTTATTGG 2000 2500
EOF

$IPCRESS $IPCRESSFILE $SEQUENCEFILE > $OUTPUTFILE
if [ $? -eq 0 ]
then
    echo Ipcress run successfully
else
    echo Problem running ipcress
    clean_exit 1
fi

TOTAL_PRODUCTS=`grep -c '^ipcress:' $OUTPUTFILE`
if [ $TOTAL_PRODUCTS -eq 1 ]
then
    echo Found 1 product, as expected
else
    echo Error: unexpectely found $TOTAL_PRODUCTS products
    clean_exit 1
fi

clean_exit 0

