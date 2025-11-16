'''
avg query_time & result_num directly
avg on query id, features: graph, degree, query_size, label_size
'''

import pandas as pd
import os

graphs = ["citeseer", "dblp", "HPRD", "human", "maayan-figeys", "twitch", "web-Stanford", "wordnet-words", "YeastS", "youtube"]
graphpath = "/var/lib/docker/subgraph/output/bsx/degree/"
labelsizes = ["L15","L30","L45","L60"]
querysizes = [20, 30, 40 ,50]
degrees = [1, 2, 3, 4, 5]
num_method = 6

glsout = pd.DataFrame(columns=['graph', 'label', 'size_', 'degree', 'filter_', 'order', 'engine', 'time', 'num', 'eps'])
gls_idx = 0

for graph in graphs:
  for labelsize in labelsizes:
    if not os.path.exists(graphpath+graph+'/'+labelsize+'/result.csv'):
      continue
    data = pd.read_csv(graphpath+graph+'/'+labelsize+'/result.csv', na_values=['NULL'], keep_default_na=False)
    filters = data['filter_'].drop_duplicates().tolist()
    orders = data['order'].drop_duplicates().tolist()
    engines = data['engine'].drop_duplicates().tolist()
    ################ filter the same results ##############
    same_data = pd.DataFrame(columns=data.columns)
    grouped = data.groupby(['size_', 'degree', 'id'])
    for name, group in grouped:
      if (len(group) != num_method):
        continue
      data_filtered = group[group['unsolved'] == 0]
      if len(data_filtered) >1 and data_filtered['result_num'].nunique() != 1:
        continue
      same_data = pd.concat([same_data, group], ignore_index=True)

    ################# construct avg data ##################
    for size_ in querysizes:
      for degree in degrees:
        for filter_ in filters:
          for order in orders:
            for engine in engines:
              tmp = same_data[(same_data.filter_ == filter_)&(same_data.order == order)&(same_data.engine == engine)&\
                              (same_data.size_ == size_)&(same_data.degree == degree)]
              if tmp.shape[0] == 0:
                continue
              print(graph+", "+labelsize+", "+str(size_)+", "+str(degree)+", "+filter_+", "+order+", "+engine+"\n")
              allrow = [graph, labelsize, size_, degree, filter_, order, engine]
              # tmp = tmp[['total_time', 'call_cnt', 'result_num', 'unsolved', 'avg_candidate_size']]
              tmp = tmp[["total_time", 'result_num']]
              # Converting each column to numeric type
              results = []
              for column in tmp.columns:
                tmp[column] = pd.to_numeric(tmp[column], errors='coerce')
                results.append(sum(tmp[column]))
              # mean total_num and total_time, then compute eps
              # appendvalues.append(appendvalues[2]/appendvalues[0])
              results.append(results[1]/results[0])
              allrow += results
              glsout.loc[gls_idx] = allrow
              gls_idx += 1
              # print(out.shape)
glsout.to_csv('degree.csv', index=False, float_format='%.6f')
