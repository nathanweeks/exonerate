#!/bin/sh
# Test that the fasta utilities skip ';' comment lines (the original
# Pearson FASTA comment convention, also used by FAIR-bioHeaders),
# both before the first record and between records.

FASTALENGTH="../../src/util/fastalength"
FASTAINDEX="../../src/util/fastaindex"
FASTAFETCH="../../src/util/fastafetch"
FASTA2ESD="../../src/util/fasta2esd"
ESD2ESI="../../src/util/esd2esi"
INPUTFILE="fasta.comment.test.fasta"
INDEXFILE="fasta.comment.test.idx"
ESDFILE="fasta.comment.test.esd"
ESIFILE="fasta.comment.test.esi"
OUTPUTFILE="fasta.comment.test.out"

clean_exit(){
    rm -f $INPUTFILE $INDEXFILE $ESDFILE $ESIFILE $OUTPUTFILE
    exit $1
    }

# Prevent fatal debugging errors
unset EXONERATE_DEBUG

( echo ';comment'
  echo ';~schema: https://example.org/fhr.json'
  cat ../data/protein/calm.human.protein.fasta
  echo ';comment between records with LETTERS'
  cat ../data/protein/tube.drome.protein.fasta ) > $INPUTFILE

# fastalength: comment letters must not be counted
$FASTALENGTH $INPUTFILE > $OUTPUTFILE
if [ $? -ne 0 ]
then
    echo "Error running $FASTALENGTH"
    clean_exit 1
fi
if [ "`cat $OUTPUTFILE | tr '\n' ' '`" = "149 CALM_HUMAN 462 TUBE_DROME " ]
then
    echo "fastalength comment test OK"
else
    echo "Unexpected fastalength output:"
    cat $OUTPUTFILE
    clean_exit 1
fi

# fastaindex/fastafetch: offsets must still point at the records
$FASTAINDEX $INPUTFILE $INDEXFILE
if [ $? -ne 0 ]
then
    echo "Error running $FASTAINDEX"
    clean_exit 1
fi
for ID in CALM_HUMAN TUBE_DROME
do
    $FASTAFETCH $INPUTFILE $INDEXFILE $ID > $OUTPUTFILE
    if [ $? -ne 0 ] || [ "`head -1 $OUTPUTFILE | cut -d' ' -f1`" != ">$ID" ] \
                    || grep -q ';' $OUTPUTFILE
    then
        echo "Unexpected fastafetch output for $ID:"
        cat $OUTPUTFILE
        clean_exit 1
    fi
done
echo "fastaindex fastafetch comment test OK"

# fasta2esd/esd2esi: the comments must not make the file irregular
$FASTA2ESD $INPUTFILE $ESDFILE > /dev/null 2>&1 \
    && $ESD2ESI $ESDFILE $ESIFILE > /dev/null 2>&1
if [ $? -ne 0 ]
then
    echo "Error running $FASTA2ESD / $ESD2ESI"
    clean_exit 1
fi
echo "fasta2esd esd2esi comment test OK"

# A comment inside a record's sequence makes the file positions of the
# sequence irregular, so fasta2esd must refuse it.
rm -f $ESDFILE
printf '>a\nACGTACGTAC\n;comment\nACGTACGTAC\nACG\n' > $INPUTFILE
$FASTA2ESD $INPUTFILE $ESDFILE > /dev/null 2>&1
if [ $? -eq 0 ]
then
    echo "fasta2esd accepted a comment inside a sequence"
    clean_exit 1
fi
echo "fasta2esd irregular comment test OK"

clean_exit 0
