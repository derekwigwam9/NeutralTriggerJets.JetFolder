#!/bin/bash
# 'GenerateXML.sh'
#
# This generates the job description file for star-submit.  Is called by the
# script 'SubmitFolding.sh'.

sim=$1
cwd=$2
ver=$3
mac=$4
arg=$5
job=$6


name=$job".job.xml"
if [ -f $name ]; then
  printf "  GenerateXML.sh: WARNING! $name already exists.\n"
  printf "                  Deleting and recreating...\n"
  rm $name
fi

touch $name
printf "<?xml version=\"1.0\" encoding=\"utf-8\" ?>\n" >> $name
printf "<job simulateSubmission=\"$sim\" fileListSyntax=\"paths\">\n" >> $name
printf "\n" >> $name
printf "  <command>\n" >> $name
printf "    cd $cwd\n" >> $name
printf "    starver $ver\n" >> $name
printf "    root -b -q $mac$arg\n" >> $name
printf "  </command>\n" >> $name
printf "\n" >> $name
printf "  <stdout URL=\"file:$cwd/$job.out\" />\n" >> $name
printf "  <stderr URL=\"file:$cwd/$job.err\" />\n" >> $name
printf "\n" >> $name
printf "</job>" >> $name
