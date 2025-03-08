#!/bin/bash

# executable for running WISC_MVPA with Apptainer on CHTC. 

# DEFINE FUNCTIONS

# this function removes files from the compute node after the script finishes (whether it finishes successfully, is manually stopped, or crashes). As input it takes a list of all arguments "$@"
cleanup(){
	# for all data files (i.e. every argument in the list of arguments), remove it
	for x in "$@"; do
		rm ./$(basename "$x")
	done
	# declare finished
	echo "all clean"	
}

# this function handles thing when the script stops because of an error
abort(){
	# print this message to the file indicated by 2, which is the .err file
	echo >&2 '
	*************
	** ABORTED **
	*************
	'
	# also print which files are present on the compute node after finishing. List detailed information about each (-l)
	echo >&2 "Files at time of error"
	echo >&2 "----------------------"
	ls >&2 -l

	# cleanup
	cleanup "$@"

	# exit with code 1, to indicate unsuccessful
	echo "An error occurred. Exiting..."
	exit 1	
}

# this function handles things when you stop the script manually
terminated(){
	# print this message to the file indicated by 2, which is the .err file
	echo >&2 '
	****************
	** TERMINATED **
	****************
	'
	# also print which files are present on the compute node after finishing. List detailed information about each (-l)
	echo >&2 "Files at time of interrupt"
	echo >&2 "--------------------------"
	ls >&2 -l

	# cleanup
	cleanup "$@"

	# exit with code 1, to indicate unsuccessful
	echo "Job terminated. Exiting..."
	exit 1	
}

# this function handles things when the script runs successfully
success(){
	# print this message to the file indicated by 2, which is the .err file
	echo >&2 '
	*************
	** SUCCESS **
	*************
	'
	
	# cleanup
	cleanup "$@"

	# exit with code 0, to indicate successful
	exit 0
	
}

# MAIN SCRIPT

# if any command fails, exit immediately (with exit code 1 for unsuccessful)
set -e
# on exit, consider the exit code. If it is 1, call the abort function. If it is 1, call the success function
trap '[[ $? -ne 0]] && abort "$@" || success "$@" ' EXIT
# if the script is instead terminated manually, call the terminated function
trap 'terminated "$@" ' SIGTERM
# print every command run. This helps with debugging
set -x

# For all data files (i.e. every argument in the list of arguments), copy them to the submit node
for x in "$@"; do
	cp "$x" .
done

# run the analysis
/WISC_MVPA

# clean up
cleanup "$@"






