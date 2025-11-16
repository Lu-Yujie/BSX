#!/bin/bash
# bash script for real graph query

exefile="/home/yujielu/subgraph/bsx/code/BSX/build/bin/BS"
file_suffix=".txt"
time_limit=300000
max_num="MAX"
num=100
parallel_num=25
# before change the engine, move the experimental results to other directory. Otherwise, the results will append to the exist files.
engine="BSX"  # "BS1"

for graph in "maayan-figeys" "YeastS" "citeseer"; do
  dataspace="/var/lib/docker/subgraph/dataset/slabel/nlabel/${graph}"
  for lsize in 15; do
    data_graph="${dataspace}/L${lsize}/${graph}-${lsize}${file_suffix}"
    for size in 10; do
      for ((j1=0; j1<=num/parallel_num; j1+=1)); do
        for ((j2=1; j2<=parallel_num&&j1*parallel_num+j2<=num; j2+=1)); do
          i=$((j1*parallel_num+j2));
          query_graph="L${lsize}/Q${size}/Q${size}-${i}"
          timeout 20 ${exefile} -d ${data_graph} -q "${dataspace}/${query_graph}${file_suffix}" \
                         -num $max_num -time_limit $time_limit -filter LDF -order GQL -engine ${engine}
        done
        wait
      done
    done
  done
done

# how to stop
# ps -ef | grep -E "SubgraphMatching|duplicate.sh" | grep -v grep | cut -c 9-16 | xargs kill -9

# clear old outputs
# find . -type f -name "*.csv" -delete
