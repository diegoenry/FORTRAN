#!/bin/bash  
awk '!/AA30/ {print $2}' apagar.ifp > clean.ifp

rm -rf ifp_split.ifp
while read line ; do
  echo $line | awk 'BEGIN {ORS=","} {for (i=1;i<=length($0);i+=5) print substr( $0, i, 5 )}' >> ifp_split.ifp
done <clean.ifp
