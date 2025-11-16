#!/bin/bash
# combine parellel exp result

workspace="/home/yujielu/subgraph/bsx/code/BSX/test/postprocess/overall"
dataspace="/var/lib/docker/subgraph/output/bsx/synthetic/EvoGraph"
for graph in "citeseer"; do
  for feature in 5 ; do
    python ${workspace}/combine.py ${dataspace}/${graph}/${graph}_${feature}
  done
done
