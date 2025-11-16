#!/bin/bash
# bash script for real graph query

exefile="/home/yujielu/subgraph/bsx/code/BSX/build/bin/BS"     # build/bin/BS
file_suffix=".txt"                     # suffix of data file
time_limit=1000                        # time limit for each query (ms) (for each method: filter-order-engine)
max_num="MAX"                          # output limit, "'MAX' or integer"
num=1000                               # the number of query graphs
parallel_num=40                        # the number of threads
con_path="/home/yujielu/subgraph/bsx/code/BSX/test/conf"           # path to your selected methods

for graph in "maayan-figeys" "YeastS" "citeseer" "HPRD" "human" "wordnet-words" "dblp" "web-Stanford" "youtube" "twitch"; do

  dataspace="/var/lib/docker/subgraph/dataset/slabel/nlabel/${graph}"    # path to your data graph
  outputspace="/var/lib/docker/subgraph/output/bsx/empty/${graph}"       # path to your output directory

  for lsize in 15 30 45 60; do                      # graph label sizes: 15 30 45 60
    data_graph="${dataspace}/L${lsize}/${graph}-${lsize}${file_suffix}"
    for size in 10 20 30 40 50; do                  # query sizes: 10 20 30 40 50
      for ((j1=4; j1<=num/parallel_num; j1+=1)); do # j1 & j2 used for parallel
        for ((j2=1; j2<=parallel_num&&j1*parallel_num+j2<=num; j2+=1)); do
          i=$((j1*parallel_num+j2));
          query_graph="L${lsize}/Q${size}/Q${size}-${i}"
          output_file="L${lsize}/Q${size}/${j2}.csv"
          timeout 20 ${exefile} -d ${data_graph} -q "${dataspace}/${query_graph}${file_suffix}" \
                         -num $max_num -time_limit $time_limit -o "${outputspace}/${output_file}" -conf "${con_path}" &
        done
        wait
      done
    done
  done
done

# how to stop
# ps -ef | grep -E "BS|overall.sh" | grep -v grep | cut -c 9-16 | xargs kill -9

# clear old outputs
# find . -type f -name "*.csv" -delete
