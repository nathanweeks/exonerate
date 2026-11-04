#!/bin/sh
# Test that FASTA input starting with ';' comment lines (the original
# Pearson FASTA comment convention, also used by FAIR-bioHeaders) is read
# as FASTA, rather than as a list of file or server names.

EXONERATE="../../src/program/exonerate"

QUERYFILE="exonerate.comment.test.query.fasta"
TARGETFILE="exonerate.comment.test.target.fasta"
CQUERYFILE="exonerate.comment.test.cquery.fasta"
CTARGETFILE="exonerate.comment.test.ctarget.fasta"
CRLFFILE="exonerate.comment.test.crlf.fasta"
NONLFILE="exonerate.comment.test.nonl.fasta"
ONLYFILE="exonerate.comment.test.only.fasta"
LISTFILE="exonerate.comment.test.list"
OUTPUTFILE="exonerate.comment.test.out"
EXPECTFILE="exonerate.comment.test.expect"

clean_exit(){
    rm -f $QUERYFILE $TARGETFILE $CQUERYFILE $CTARGETFILE $CRLFFILE \
          $NONLFILE $ONLYFILE $LISTFILE $OUTPUTFILE $EXPECTFILE
    exit $1
    }

cat > $QUERYFILE << SEQEOF
>cdna
ATGCAAGGTCTAGCTAGGCATTACGGATCCATGAACCTAGGTTCCAAGACTGCATGGACTAGCATCGATCAGCTAGCATG
SEQEOF

cat > $TARGETFILE << SEQEOF
>genome
ATGCAAGGTCTAGCTAGGCATTACGGATCCATGAACCTAGGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACAGGTTCCAAGACTGCATGGACTAGCATCGATCAGCTAGCATG
SEQEOF

# Leading comment block, including a key:value line (which must not be
# taken as a server name) and a '>' that is not at the start of a line.
( echo ';comment'
  echo ';~schema: https://example.org/fhr.json'
  echo ';> not a header'
  echo
  cat $QUERYFILE ) > $CQUERYFILE
# The extra short record is too short to align, so the expected
# output is the same as for $TARGETFILE.
( echo ';comment'
  cat $TARGETFILE
  echo ';comment between records'
  echo '>short'
  echo 'ACGTACGTAC' ) > $CTARGETFILE
# CRLF line endings
printf ';comment\r\n>cdna\r\nATGCAAGGTCTAGCTAGGCATTACGGATCCATGAACCTAGGTTCCAAGACTGCATGGACTAGCATCGATCAGCTAGCATG\r\n' > $CRLFFILE

run_exonerate(){
    $EXONERATE --model est2genome --showalignment no --showvulgar yes \
               $1 $2 > $OUTPUTFILE
    if [ $? -ne 0 ]
    then
        echo "Problem running exonerate on $1 $2"
        clean_exit 1
    fi
    grep '^vulgar:' $OUTPUTFILE > $OUTPUTFILE.vulgar
    mv $OUTPUTFILE.vulgar $OUTPUTFILE
    }

run_exonerate $QUERYFILE $TARGETFILE
mv $OUTPUTFILE $EXPECTFILE
if [ ! -s $EXPECTFILE ]
then
    echo "No alignment for uncommented input"
    clean_exit 1
fi

check(){
    run_exonerate $1 $2
    if cmp -s $OUTPUTFILE $EXPECTFILE
    then
        echo "Comment test OK: $3"
    else
        echo "Unexpected result: $3"
        cat $OUTPUTFILE
        clean_exit 1
    fi
    }

check $CQUERYFILE $TARGETFILE "commented query"
check $QUERYFILE $CTARGETFILE "commented target"
check $CQUERYFILE $CTARGETFILE "commented query and target"
check $CRLFFILE $TARGETFILE "CRLF commented query"

# A file of sequence names listing commented FASTA files: the comments at
# the start of the second file must be skipped by the FASTA reader.
echo $TARGETFILE > $LISTFILE
echo $CTARGETFILE >> $LISTFILE
run_exonerate $QUERYFILE $LISTFILE
if [ `wc -l < $OUTPUTFILE` -eq 2 ]
then
    echo "Comment test OK: list of commented files"
else
    echo "Unexpected result: list of commented files"
    cat $OUTPUTFILE
    clean_exit 1
fi

# The same, where the first file listed has no trailing newline.
printf '>cdna\nATGCAAGG' > $NONLFILE
echo $NONLFILE > $LISTFILE
echo $CTARGETFILE >> $LISTFILE
run_exonerate $QUERYFILE $LISTFILE
if [ `wc -l < $OUTPUTFILE` -eq 1 ]
then
    echo "Comment test OK: list with file lacking trailing newline"
else
    echo "Unexpected result: list with file lacking trailing newline"
    cat $OUTPUTFILE
    clean_exit 1
fi

# A comments-only file is a fasta file without sequences
# (not a list of file names).
printf ';only\n;comments\n' > $ONLYFILE
$EXONERATE $QUERYFILE $ONLYFILE > $OUTPUTFILE 2>&1
if [ $? -ne 0 ] && grep -q 'No sequences found' $OUTPUTFILE
then
    echo "Comment test OK: comments-only file"
else
    echo "Unexpected result: comments-only file"
    cat $OUTPUTFILE
    clean_exit 1
fi

clean_exit 0
