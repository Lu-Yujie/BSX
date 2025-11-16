#!/bin/bash
# combine parellel exp result

workspace="/home/yujielu/subgraph/bsx/code/BSX/test/postprocess/overall"
dataspace="/var/lib/docker/subgraph/output/bsx/similar"
for graph in "citeseer" "dblp" "HPRD" "human" "maayan-figeys" "twitch" "web-Stanford" "wordnet-words" "YeastS" "youtube"; do
  for sim in 0 25 50 75 100; do
    for label in 15 30 45 60; do
      python ${workspace}/combine.py ${dataspace}/s${sim}/${graph}/L${label}
    done
  done
done
