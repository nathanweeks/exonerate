#!/bin/sh
# test for fastacomposition utilility

PROGRAM="../../src/util/fastacomposition"
INPUTFILE="../data/cdna/calm.human.dna.fasta"

$PROGRAM $INPUTFILE | (read FILENAME A ACOMP C CCOMP G GCOMP T TCOMP
   if [ $ACOMP -eq 430 ] \
   && [ $CCOMP -eq 626 ] \
   && [ $GCOMP -eq 592 ] \
   && [ $TCOMP -eq 527 ]
   then
      echo Calmodulin cDNA composition correct
   else
      echo Calmodulin cDNA composition wrong $ACOMP $CCOMP $GCOMP $TCOMP
      exit 1
   fi
   )

exit $?

