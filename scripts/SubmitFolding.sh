#!/bin/bash
# 'SubmitFolding.sh'
# Dere Anderson
#
# Use this to create xml scripts and submit them via 'star-submit'.
#
# NOTE: This script requires 'GenerateXML.sh' and 'GenerateDir.sh'


# submission parameters
sim="false"
ver="pro"
mac="VaryFolding.C"

# output parameters
top="VaryFolding"
dat="ChargedParticlePt"

# folding parameters
P=(0)
M=(0)
K=(0)
N=(5.8)
T=(0.4)



cwd=$PWD
printf "Running submission script...\n"
for p in ${P[@]}; do
  for m in ${M[@]}; do
    for k in ${K[@]}; do
      for n in ${N[@]}; do
        for t in ${T[@]}; do

          # create sub directory (if necessary)
          if [ $m == "0" ]; then
            met="NoUnfolding"
          elif [ $m == "1" ]; then
            met="Bayes"
          elif [ $m == "2" ]; then
            met="Svd"
          fi
          ./GenerateDir.sh $top $dat $met
          sub=$cwd"/"$top"/"$dat"/"$met

          nTxt=$(echo "scale=0; $n*10/1" | bc)
          tTxt=$(echo "scale=0; $t*10/1" | bc)
          job="p"$p"m"$m"k"$k"n"$nTxt"t"$tTxt
          arg="'($p,$m,$k,$n,$t,\"$sub\")'"
          ./GenerateXML.sh $sim $sub $ver $mac $arg $job

          cp $mac $sub
          mv $job".job.xml" $sub
          cd $sub
          printf "  Submitting job '$job'...\n"
          printf "\n"
          star-submit $job".job.xml"
          printf "\n"
          cd $cwd

        done  # end T loop
      done  # end N loop
    done  # end K loop
  done  # end M loop
done  # end P loop

printf "Finished submitting!\n"
