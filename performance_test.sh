#!/bin/bash

# Author: Raphael Ruban
# Runs either the sequential or parallel version of the MPAS Regent code.
# Prints out the time of each run into a file called sequential_times.txt or parallel_times.txt.
# Defaults to parallel. Use -s for sequential.

############################################################
# Help                                                     #
############################################################
Help()
{
  # Display Help
  echo "This script runs either the sequential or parallel version of the MPAS Regent code."
  echo "Defaults to running the parallel version and saves the results to parallel_times.txt."
  echo
  echo "Syntax: scriptTemplate [-h|n|s]"
  echo "options:"
  echo "h     Print this Help."
  echo "n     Set the number of times MPAS should be run. Defaults to 30."
  echo "s     Run MPAS without CUDA support and save times in sequential_times.txt."
  echo
}

############################################################
# Main program                                             #
############################################################
# Set default variables.
num=30
cuda=true
file_name=parallel_times.txt

# Get the options
while getopts ":hn:s" option; do
  case $option in
    h) # display Help
      Help
      exit;;
    n) # Enter number of runs.
      num=$OPTARG;;
    s) # Set to sequential.
      file_name=sequential_times.txt
      cuda=false;;
    \?) # Invalid option
      echo "Error: Invalid option. Use -h for more info."
      exit;;
  esac
done

# Performance test
if [ "$cuda" = true ]
then
  echo "Starting parallel performance test..."
else
  echo "Starting sequential performance test..."
fi

echo "Start of performance script." >> $file_name
for i in {1..$num}
do
    echo "$i"
    echo "Run $i:" >> $file_name
    if [ "$cuda" = true ]
    then
      (time -p ( ../legion/language/regent.py ~/mpas-regent/main.rg -fcuda 1 -ll:gpu 4 &> /dev/null)) 2>> parallel_times.txt
    else
      (time -p (/home/zengcs/regent ~/mpas-regent/main.rg &> /dev/null)) 2>> sequential_times.txt
    fi
    echo "" >> $file_name
done

echo "All done!"
