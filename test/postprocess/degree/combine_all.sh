#!/bin/bash
# combine parellel exp result

dataspace="/var/lib/docker/subgraph/output/bsx/degree"
for graph in "citeseer" "dblp" "HPRD" "human" "maayan-figeys" "twitch" "web-Stanford" "wordnet-words" "YeastS" "youtube"; do
  for label in 15 30 45 60; do
    python combine.py ${dataspace}/${graph}/L${label}
  done
done
