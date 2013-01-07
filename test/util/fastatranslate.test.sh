#!/bin/sh
# test for fastatranslate utility

FASTATRANSLATE="../../src/util/fastatranslate"
FASTASUBSEQ="../../src/util/fastasubseq"
FASTADIFF="../../src/util/fastadiff"

CDNA="../data/cdna/calm.human.dna.fasta"
PROTEIN="../data/protein/calm.human.protein.fasta"

CDS="fastatranslate.test.cds.fasta"
TRANSLATED="fastatranslate.test.translated.fasta"

clean_exit(){
    rm -rf $CDS $TRANSLATED
    exit $1
    }

$FASTASUBSEQ $CDNA --start 103 --length 447 > $CDS
if [ $? -eq 0 ]
then
    echo Extraced CDS for translation
else
    echo Problem extracting CDS for tranlation
    clean_exit 1
fi

$FASTATRANSLATE $CDS --frame 1 > $TRANSLATED
if [ $? -eq 0 ]
then
    echo Tranlated sequence
else
    echo Failed to tranlate sequence
    clean_exit 1
fi

$FASTADIFF --checkids no $TRANSLATED $PROTEIN
if [ $? -eq 0 ]
then
    echo Tranlsated sequence correct
else
    echo Translated sequence is wrong
    clean_exit 1
fi

clean_exit 0

