#!/bin/bash
# 'clean.sh'
#
# Clean up your directory!

if [ $# -gt 1 ]
then
  echo "Hey! Only pass ONE argument!"
else

  echo "Cleaning '$PWD'..."
  if [ $1 == "deep" ]
  then
    echo "  Time for a deep scrub..."
    rm *.log
    rm *.out
    rm *.err
    rm *.report
    rm *.list
    rm sched*.csh
    rm array*.csh*
    rm *.session.xml
    echo "  All clean!"
  else
    echo "  Just a quick scrub..."
    rm *.log
    rm *.out
    rm *.err
    rm *.report
    echo "  All clean!"
  fi

fi
