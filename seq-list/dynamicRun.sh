#!/bin/bash
i=16
echo $i;
until [ $i -gt 19 ]
do
  echo $i
  sed -i.bak "s,.*#define default_SCALE.*,#define default_SCALE ((int64_t)$i),g" ../options.h &
  sed -i.bak "s,.*#define default_edgefactor.*,#define default_edgefactor ((int64_t)$i),g" ../options.h &
  wait
  cd .. && make clean && make && cd seq-list && ./seq-list | egrep "edgefactor|median_time|median_TEPS" >> out.txt &
  wait
  i=$(( i+1 ))
done
