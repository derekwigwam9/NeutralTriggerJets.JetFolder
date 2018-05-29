#!/bin/bash
# 'GenerateDir.sh'
#
# This generates directories for the output of 'SubmitFolding.sh'.

topDir=$1
data=$2
meth=$3
cwd=$PWD


# create top output dir. if necessary
outDir=$cwd"/"$topDir
if [ ! -d $outDir ]; then
  printf "  GenerateDir.sh: Creating directory './$topDir'\n"
  mkdir $outDir
fi

outData=$outDir"/"$data
if [ ! -d $outData ]; then
  printf "  GenerateDir.sh: Creating directory './$topDir/$data'\n"
  mkdir $outData
fi

outMeth=$outDir"/"$data"/"$meth
if [ ! -d $outMeth ]; then
  printf "  GenerateDir.sh: Creating directory './$topDir/$data/$meth'\n"
  mkdir $outMeth
fi
