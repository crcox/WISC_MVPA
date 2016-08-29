#!/bin/bash

reverse=''

function usage()
{
  printf "identify best soslasso paramaters from csv\n"
  printf "\n"
  printf "./pick_alpha_lambda.sh [options] <results.csv>\n"
  printf "\t-h --help\n"
  printf "\t-r,--reverse  Default is to search for best alpha and lambda from\n"
  printf "\t              smallest to largest. Reverse means to search from\n"
  printf "\t              largest to smallest.  Matters in case of ties.\n"
  printf "\n"
}

SHORT=hr
LONG=help,reverse
PARSED=`getopt --options $SHORT --longoptions $LONG --name "$0" -- "$@"`
if [[ $? != 0 ]]; then
  # e.g. $? == 1
  #  then getopt has complained about wrong arguments to stdout
  exit 2
fi

eval set -- "$PARSED"

while true; do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    -r|--reverse)
      reverse='r'
      shift
      ;;
    --)
      shift
      break
      ;;
    *)
      echo "Programming error"
      exit 3
      ;;
  esac
done

# handle required positional arguments
if [[ $# != 1 ]]; then
    echo "$0: A single input file is required."
    echo ""
    usage
    exit 4
fi

results=$1
mean_by="mean_diff_by_alpha_lambda_cv.csv"
sort_cmd="sort -n${reverse}"

awk -f mean_diff_by_alpha_lambda_cv.awk "$results"|${sort_cmd}> "$mean_by"
awk -f best_alpha_lambda_by_cv.awk "$mean_by"
