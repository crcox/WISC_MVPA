#!/usr/bin/awk

BEGIN{
  FS=",";
  OFS=",";
}

NR>1 {
  arr[$4","$6","$7]+=$11;
  arrC[$4","$6","$7]++
}

END {
  for (a in arr) {
    print a, arr[a]/arrC[a]
  }
}
