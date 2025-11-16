'''
avg query_time & result_num directly
avg on three level: query_size, label_size, graph
'''
import os
import csv
import sys

import pandas as pd

graphs = ["citeseer", "dblp", "HPRD", "human", "maayan-figeys", "twitch", "web-Stanford", "wordnet-words", "YeastS", "youtube"]
datapath = "/var/lib/docker/subgraph/output/bsx/global/o_10_5/"
labelsizes = ["L15","L30","L45","L60"]

glsout = pd.DataFrame(columns=['graph', 'label', 'graph_id', 'method', 'time', 'num', 'eps'])
gls_idx = 0

for graph in graphs:
  for labelsize in labelsizes:
    for graph_id in range(1,9):
      input_file = datapath+graph+'/'+labelsize+'/q'+str(graph_id)+'.csv'
      with open(input_file, 'r') as file:
        reader = csv.reader(file)
        for line in reader:
          if len(line) <= 0 or not line[0].startswith('q'):
            continue
          row = [graph, labelsize, graph_id, line[4], line[5], line[6]]
          # Converting each column to numeric type
          row[-1] = pd.to_numeric(row[-1], errors='coerce')
          row[-2] = pd.to_numeric(row[-2], errors='coerce')
          # mean total_num and total_time, then compute eps
          row.append(row[-1]/row[-2])
          glsout.loc[gls_idx] = row
          gls_idx += 1
glsout.to_csv('o_10_5.csv', index=False, float_format='%.6f')
