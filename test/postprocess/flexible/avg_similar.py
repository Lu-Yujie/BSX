'''
avg query_time & result_num directly
avg on three level: query_size, label_size, graph
'''

import pandas as pd

graphs = ["citeseer", "dblp", "HPRD", "human", "maayan-figeys", "twitch", "web-Stanford", "wordnet-words", "YeastS", "youtube"]
datapath = "/var/lib/docker/subgraph/output/bsx/similar/"
labelsizes = ["L15", "L30", "L45", "L60"]
querysizes = [10, 20, 30, 40 ,50]
sims = ["s0", "s25", "s50", "s75", "s100"]
glsout = pd.DataFrame(columns=['sim', 'graph', 'label', 'size_', 'time', 'num', 'eps'])
gls_idx = 0

for sim in sims:
  for graph in graphs:
    # 分label, size, 先构建整体数据
    for labelsize in labelsizes:
      data = pd.read_csv(datapath+sim+'/'+graph+'/'+labelsize+'/result.csv', na_values=['NULL'], keep_default_na=False)
      ################# construct avg data ##################
      for size_ in querysizes:
        tmp = data[(data.filter_ == 'CFL')&(data.size_ == size_)]
        if tmp.shape[0] == 0:
          continue
        print(sim+", "+graph+", "+labelsize+", "+str(size_))
        allrow = [sim, graph, labelsize, size_]
        # tmp = tmp[['total_time', 'call_cnt', 'result_num', 'unsolved', 'avg_candidate_size']]
        tmp = tmp[["total_time", 'result_num']]
        # Converting each column to numeric type
        results = []
        for column in tmp.columns:
          tmp[column] = pd.to_numeric(tmp[column], errors='coerce')
          results.append(sum(tmp[column]))
        # mean total_num and total_time, then compute eps
        results.append(results[1]/(results[0]))
        allrow += results
        glsout.loc[gls_idx] = allrow
        gls_idx += 1
        # print(out.shape)
glsout.to_csv('similar.csv', index=False, float_format='%.6f')
