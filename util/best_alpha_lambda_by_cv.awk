#!/usr/bin/awk

BEGIN {
  FS=",";
  OFS=","
}

{
  if (arr[$1] < $4) {
    arrOUT[$1]=$2","$3","$4;
    arr[$1]=$4;
  }
}

END {
  for (a in arr) print a, arrOUT[a]
}
