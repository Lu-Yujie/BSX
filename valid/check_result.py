# a scripts used for comparing the results of different methods
# assume your results file organized as the following format:
# size_,id,unsolved,filter_,order,engine,total_time,result_num
# 10,1,1,RM,RM,RM,0.999981,44235106
# 10,1,1,DPiso,DPiso,KSS,0.999988,4445064
# 10,1,0,CFL,null,BSX,0.0904933,242282283
# 10,1,0,VEQ,VEQ,VEQ,0.227992,242282283
# size_: query graph size
# id: query graph id
# unsolved: whether find all results in given time
# filter_,order,engine: your method
# total_time,result_num: time and #embeddings

import pandas as pd

graphs = ["citeseer", "dblp", "HPRD", "human", "maayan-figeys", "web-Stanford", "wordnet-words", "YeastS", "youtube", "twitch"]
graphpath = "/path/to/your/outputs"
labelsizes = ["L15","L30","L45","L60"]
querysizes = [10, 20, 30, 40 ,50]
results = []

for graph in graphs:
  results.clear()
  for labelsize in labelsizes:
    data = pd.read_csv(graphpath+graph+'/'+labelsize+'/result.csv')
    data_filtered = data[data['unsolved'] == 0]
    grouped = data_filtered.groupby(['size_', 'id'])
    for name, group in grouped:
      group['label'] = labelsize
      if group['result_num'].nunique() > 1:
          results.append(group)
  if results:
    result_df = pd.concat(results)
    result_df.to_csv(graph+'.csv', index=False)
    print(graph+" has error results")
  else:
    print(graph+" check ok")
