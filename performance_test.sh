1   #!/bin/bash
  1
  2 # Author: Raphael Ruban
  3 # Runs either the sequential or parallel version of the MPAS Regent code.
  4 # Prints out the time of each run into a file called sequential_times.txt or parallel_times.txt.
  5 # Defaults to parallel. Use -s for sequential.
  6
  7 ############################################################
  8 # Help                                                     #
  9 ############################################################
 10 Help()
 11 {
 12   # Display Help
 13   echo "This script runs either the sequential or parallel version of the MPAS Regent code."
 14   echo "Defaults to running the parallel version and saves the results to parallel_times.txt."
 15   echo
 16   echo "Syntax: scriptTemplate [-h|n|s]"
 17   echo "options:"
 18   echo "h     Print this Help."
 19   echo "n     Set the number of times MPAS should be run. Defaults to 30."
 20   echo "s     Run MPAS without CUDA support and save times in sequential_times.txt."
 21   echo
 22 }
 23
 24 ############################################################
 25 # Main program                                             #
 26 ############################################################
 27 # Set default variables.
 28 num=30
 29 cuda=true
 30 file_name=parallel_times.txt
 31
 32 # Get the options
 33 while getopts ":hn:s" option; do
 34   case $option in
 35     h) # display Help
 36       Help
 37       exit;;
 38     n) # Enter number of runs.
 39       num=$OPTARG;;
 40     s) # Set to sequential.
 41       file_name=sequential_times.txt
 42       cuda=false;;
