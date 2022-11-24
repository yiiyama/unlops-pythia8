#!/bin/bash

# Usage: execute.sh JOBID NEVENTS

# This script is set up to be run in a batch environment, with an integer uniquely defining the job (JOBID)
# and the number of events to run over in each job (NEVENTS) as two arguments. LHE files must be generated
# and stored somewhere beforehand, with a fixed number of events (must equal EVENTSPERFILE below) per file.
# The list of LHE files should be compiled per process (Z_[N]LO_Inclusive_(0|1|2)j.list) and placed in LISTDIR.
# The script then calculates the line number within the lists corresponding to the files that the job should
# process, copies those LHE files to $TMPDIR, extracts the specific events to process, and run the UNLOPS program.

LISTDIR=lists
EVENTSPERFILE=40000
OUTDIR=
#export X509_USER_PROXY=...
#alias copycmd=xrdcp
alias copycmd=cp

JOBID=$1
NEVENTS=$2

JOBSPERFILE=$(($EVENTSPERFILE/$NEVENTS))

echo "JOBID $JOBID"
echo "NEVENTS $NEVENTS"

ISOURCE=$(($JOBID/$JOBSPERFILE+1))

echo "SOURCE $ISOURCE"

REPO=$(dirname $(readlink -f $0))

source $REPO/env.sh

NLO_0J=$(sed -n ${ISOURCE}p $LISTDIR/Z_NLO_Inclusive_0j.list)
NLO_1J=$(sed -n ${ISOURCE}p $LISTDIR/Z_NLO_Inclusive_1j.list)
NLO_2J=$(sed -n ${ISOURCE}p $LISTDIR/Z_NLO_Inclusive_2j.list)
LO_1J=$(sed -n ${ISOURCE}p $LISTDIR/Z_LO_Inclusive_1j.list)
LO_2J=$(sed -n ${ISOURCE}p $LISTDIR/Z_LO_Inclusive_2j.list)
LO_3J=$(sed -n ${ISOURCE}p $LISTDIR/Z_LO_Inclusive_3j.list)

for FILE in $NLO_0J $NLO_1J $NLO_2J $LO_0J $LO_1J $LO_2J
do
  copycmd $FILE $TMPDIR/
done

NLO_0J=$TMPDIR/$(basename $NLO_0J)
NLO_1J=$TMPDIR/$(basename $NLO_1J)
NLO_2J=$TMPDIR/$(basename $NLO_2J)
LO_1J=$TMPDIR/$(basename $LO_1J)
LO_2J=$TMPDIR/$(basename $LO_2J)
LO_3J=$TMPDIR/$(basename $LO_3J)

START=$((($JOBID%$JOBSPERFILE)*$NEVENTS))

for FILE in $NLO_0J $NLO_1J $NLO_2J $LO_0J $LO_1J $LO_2J
do
  mv $FILE $TMPDIR/tmp.lhe
  $REPO/tools/extract_lhe.py $TMPDIR/tmp.lhe $START $NEVENTS $FILE
  rm $TMPDIR/tmp.lhe
done

sed -e "s|NLO_0J_LHE|$NLO_0J|" -e "s|NLO_1J_LHE|$NLO_1J|" -e "s|NLO_2J_LHE|$NLO_2J|" -e "s|LO_1J_LHE|$LO_1J|" -e "s|LO_2J_LHE|$LO_2J|" -e "s|LO_3J_LHE|$LO_3J|" $REPO/cfg/run_unlops_batch.cfg > $TMPDIR/run_unlops.cfg

cat $TMPDIR/run_unlops.cfg

$REPO/bin/run_unlops $TMPDIR/run_unlops.cfg $TMPDIR/unlops_$JOBID.root

copycmd unlops_$JOBID.root $OUTDIR/unlops_${JOBID}.root
