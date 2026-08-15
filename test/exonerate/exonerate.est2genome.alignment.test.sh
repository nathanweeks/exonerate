#!/bin/sh
# Test that --showalignment coordinates are correct on the sequence that
# does NOT advance across a target intron (the cDNA/query in est2genome).
# The label spanning the intron is a zero-advance unit for the query, and
# its row-boundary coordinate must report the last emitted query base, not
# the next one (https://github.com/nathanweeks/exonerate/issues/26).

EXONERATE="../../src/program/exonerate"

QUERYFILE="exonerate.est2genome.alignment.test.query.fasta"
TARGETFILE="exonerate.est2genome.alignment.test.target.fasta"
OUTPUTFILE="exonerate.est2genome.alignment.test.out"

clean_exit(){
    rm -f $QUERYFILE $TARGETFILE $OUTPUTFILE
    exit $1
    }

# Query: two 40 bp exons joined (cDNA).  Target: the same two exons in a
# genome, separated by a 50 bp GT..AG intron.  exon1 = query bases 1-40,
# exon2 = query bases 41-80.
cat > $QUERYFILE << SEQEOF
>cdna
ATGCAAGGTCTAGCTAGGCATTACGGATCCATGAACCTAGGTTCCAAGACTGCATGGACTAGCATCGATCAGCTAGCATG
SEQEOF

cat > $TARGETFILE << SEQEOF
>genome
ATGCAAGGTCTAGCTAGGCATTACGGATCCATGAACCTAGGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACAGGTTCCAAGACTGCATGGACTAGCATCGATCAGCTAGCATG
SEQEOF

# --alignmentwidth 46 -> a display width of 30 columns, so a row boundary
# falls inside the intron label: the query row showing the last exon-1
# bases plus the start of the intron label ends at query base 40, and the
# query resumes exon 2 at base 41 on the next row.  (Before the fix these
# were reported as 41 and 42.)
$EXONERATE --model est2genome --showalignment yes --showvulgar no \
           --alignmentwidth 46 $QUERYFILE $TARGETFILE > $OUTPUTFILE
if [ $? -eq 0 ]
then
    echo "Exonerate est2genome alignment test OK"
else
    echo "Problem running est2genome alignment test for exonerate"
    clean_exit 1
fi

# Query row that ends exon 1 (contains the intron label): last field is the
# query end coordinate, which must be 40 (10 bases 31-40 shown).
EXON1_END=`grep 'Target Intron' $OUTPUTFILE | awk '{print $NF}'`
# Query row that resumes exon 2 (has both the label arrows and exon-2
# sequence): first field is the query start coordinate, which must be 41.
EXON2_START=`grep '>>>>' $OUTPUTFILE | grep 'GTTCCAAGACTGCATGGACTA' | awk '{print $1}'`

if [ "$EXON1_END" = "40" ] && [ "$EXON2_START" = "41" ]
then
    echo "Query coordinates as expected: $EXON1_END / $EXON2_START"
else
    echo "Unexpected query coordinates: exon1 end=$EXON1_END exon2 start=$EXON2_START"
    clean_exit 1
fi

# The target (genome) is contiguous across the intron, so its row
# coordinates must stay continuous even though the intron is drawn
# compressed as dots: the row ending inside the intron ends at target
# base 74, and the next row must resume at 75 (not skip ahead).
TARGET_INTRON_END=`grep 'ATGAACCTAGgt' $OUTPUTFILE | awk '{print $NF}'`
TARGET_RESUME_START=`grep 'agGTTCCAAGACTGCATGGACTA' $OUTPUTFILE | awk '{print $1}'`

if [ "$TARGET_INTRON_END" = "74" ] && [ "$TARGET_RESUME_START" = "75" ]
then
    echo "Target coordinates as expected: $TARGET_INTRON_END / $TARGET_RESUME_START"
else
    echo "Unexpected target coordinates: intron end=$TARGET_INTRON_END resume start=$TARGET_RESUME_START"
    clean_exit 1
fi

clean_exit 0
