#!/bin/sh
# test that --showalignment coordinates stay correct when a codon is
# split across the row-wrap boundary in a translated (dna2protein)
# alignment (https://github.com/nathanweeks/exonerate/issues/26)

EXONERATE="../../src/program/exonerate"

QUERYFILE="exonerate.dna2protein.alignment.test.query.fasta"
TARGETFILE="exonerate.dna2protein.alignment.test.target.fasta"
OUTPUTFILE="exonerate.dna2protein.alignment.test.out"

clean_exit(){
    rm -f $QUERYFILE $TARGETFILE $OUTPUTFILE
    exit $1
    }

cat > $QUERYFILE << SEQEOF
>q
ATGGGCCGGGCCCGGCCGGGCCAACGCGGGCCGCCCAGCCCCGGCCCCGCCGCGCAGCCT
SEQEOF

cat > $TARGETFILE << SEQEOF
>t
MGRARPGQRGPPSPGPAAQP
SEQEOF

# --alignmentwidth 25 gives a wrapped row width of 11 columns, which is
# not a multiple of 3, so the 4th codon (GCC) is split across the first
# two displayed rows: only its first 2 bases fit on row 1.
$EXONERATE --model ungapped --showalignment yes --showvulgar no \
           --alignmentwidth 25 $QUERYFILE $TARGETFILE > $OUTPUTFILE
if [ $? -eq 0 ]
then
    echo "Exonerate dna2protein alignment test OK"
else
    echo "Problem running dna2protein alignment test for exonerate"
    clean_exit 1
fi

QUERY_ROW1_END=`grep -m1 -E '^ *1 : [ACGT]+ : ' $OUTPUTFILE | awk '{print $NF}'`
QUERY_ROW2_START=`grep -m1 -E '^ *12 : [ACGT]+ : ' $OUTPUTFILE | awk '{print $1}'`
TARGET_ROW1_END=`grep -m1 -E '^ *1 : MetGlyArgAl : ' $OUTPUTFILE | awk '{print $NF}'`
TARGET_ROW2_START=`grep -m1 -E '^ *[0-9]+ : aArgProGlyG : ' $OUTPUTFILE | awk '{print $1}'`

# The row should end/resume showing the actual base at that column (a
# split codon's bases still belong to consecutive, non-overlapping
# nucleotide positions): row 1 ends at base 11, row 2 resumes at base 12.
# The target is shown as three-letter amino acid names.  The same residue
# can span the row boundary, so row 2 starts with the remainder of target
# residue 4, not residue 5.
if [ "$QUERY_ROW1_END" = "11" ] && [ "$QUERY_ROW2_START" = "12" ] \
   && [ "$TARGET_ROW1_END" = "4" ] && [ "$TARGET_ROW2_START" = "4" ]
then
    echo "Coordinates as expected: query $QUERY_ROW1_END / $QUERY_ROW2_START target $TARGET_ROW1_END / $TARGET_ROW2_START"
else
    echo "Unexpected coordinates: query row1 end=$QUERY_ROW1_END row2 start=$QUERY_ROW2_START target row1 end=$TARGET_ROW1_END row2 start=$TARGET_ROW2_START"
    clean_exit 1
fi

clean_exit 0
