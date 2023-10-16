#!/bin/bash


for file in *.txt; do
  awk '/StartTrial/ { show=1 } show; /Trial End/ { show=0 }' "$file" >> "${file%.txt}_Z.txt"
done

for file in *_Z.txt; do
   awk '{print $1, $2, $3, $4, $5, $6}'  "$file" >> "${file%.txt}_A.txt"
done

for file in *_A.txt; do
  sed '/[a-zA-Z]/d' "$file" >> "${file%.txt}_B.txt"
  sed -n -e '/StartTrial/p' "$file" >  "${file%.txt}1T.txt"
  sed -n -e '/Fix ON/p' "$file" >  "${file%.txt}2T.txt"
  sed -n -e '/Adaptor ON/p' "$file" >  "${file%.txt}3T.txt"
  sed -n -e '/Stim ON/p' "$file" >  "${file%.txt}4T.txt"
  sed -n -e '/Stim OFF/p' "$file" >  "${file%.txt}5T.txt"
  sed -n -e '/Trial End/p' "$file" >  "${file%.txt}6T.txt"
done

for file in *_B.txt; do
  awk '{print $1, $2, $3, $4}' "$file" >> "${file%.txt}_C.txt"
done

for file in *_C.txt; do
  sed -n '/. . 0.0/!p' "$file" >> "${file%.txt}_D.txt"
done

for file in *T.txt; do
  awk '{print $2, $5, $6}'  "$file" >> "${file%.txt}T.txt"
done


