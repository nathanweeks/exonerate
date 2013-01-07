#!/bin/sh
# test for fastalength utility

PROGRAM="../../src/util/fastalength"
INPUTFILE="../data/protein/calm.human.protein.fasta"

$PROGRAM $INPUTFILE | (read LENGTH IDENTIFIER
   echo $LENGTH $IDENTIFIER
   if [ $LENGTH -eq 149 ]
   then
       echo Calmodulin length OK
   else
       echo Bad calmodulin length: $LENGTH
       exit 1
   fi
   if [ $IDENTIFIER = "CALM_HUMAN" ]
   then
       echo Calmodulin identifier OK
   else
       echo bad calmodulin identifier: $IDENTIFIER
       exit 1
   fi
   )

exit $?

