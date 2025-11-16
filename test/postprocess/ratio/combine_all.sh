#!/bin/bash
# combine parellel exp result

workspace="/home/yujielu/subgraph/bsx/code/BSX/test/postprocess/overall"
dataspace="/var/lib/docker/subgraph/output/bsx/ratio_label"

for graph in "citeseer", "dblp", "youtube", "HPRD", "human", "maayan-figeys", "twitch", "web-Stanford", "wordnet-words", "YeastS"; do
  for label in 10 50 100 200 1000 10000 100000; do
    python ${workspace}/combine.py ${dataspace}/${graph}/B${label}
  done
done
