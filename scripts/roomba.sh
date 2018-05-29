#!/bin/bash
# 'roomba.sh'
#
# Clean up your sub-directories!


# dirty directories
Top="VaryFolding"
Data=("ChargedParticlePt")
Meth=("NoUnfolding")

echo "Activating roomba..."

if [ $# -gt 1 ]
then
  echo "Hey! Only pass ONE argument!"
else

  cwd=$PWD
  for data in ${Data[@]}; do
    for meth in ${Meth[@]}; do

      dirtyDir=$Top"/"$data"/"$meth
      if [ ! -d $dirtyDir ]; then
        continue
      fi

      cd $dirtyDir
      cp $cwd"/"clean.sh ./cleaning.sh
      ./cleaning.sh $1
      rm cleaning.sh
      cd $cwd

    done
  done

fi

echo "roomba finished."
