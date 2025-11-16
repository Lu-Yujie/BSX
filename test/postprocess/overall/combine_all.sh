#!/bin/bash
# combine parellel exp result

workspace="/home/yujielu/subgraph/bsx/code/BSX/test/postprocess/overall"
dataspace="/var/lib/docker/subgraph/output/bsx/nlabel"
for graph in "citeseer", "dblp", "youtube", "HPRD", "human", "maayan-figeys", "twitch", "web-Stanford", "wordnet-words", "YeastS"; do
  for label in 15 30 45 60; do
    python ${workspace}/combine.py ${dataspace}/${graph}/L${label}
  done
done
