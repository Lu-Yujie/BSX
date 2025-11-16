'''
avg query_time & result_num directly
avg on three level: query_size, label_size, graph
'''

import pandas as pd

graphs = ["citeseer", "dblp", "HPRD", "human", "maayan-figeys", "twitch", "web-Stanford", "wordnet-words", "YeastS", "youtube"]
datapath = "/var/lib/docker/subgraph/output/bsx/ratio_label/"
labelsizes = ["B10", "B50", "B100", "B200", "B1000", "B10000", "B100000"]
querysizes = [10, 20, 30, 40 ,50]

glsout = pd.DataFrame(columns=['graph', 'label', 'size_', 'filter_', 'order', 'engine', 'time', 'num', 'eps'])
gls_idx = 0

for graph in graphs:
  # 分label, size, 先构建整体数据
  for labelsize in labelsizes:
    data = pd.read_csv(datapath+graph+'/'+labelsize+'/result.csv', na_values=['NULL'], keep_default_na=False)
    filters = data['filter_'].drop_duplicates().tolist()
    orders = data['order'].drop_duplicates().tolist()
    engines = data['engine'].drop_duplicates().tolist()
    ################ filter the same results ##############
    same_data = pd.DataFrame(columns=data.columns)
    grouped = data.groupby(['size_', 'id'])
    for name, group in grouped:
      data_filtered = group[group['unsolved'] == 0]
      if len(data_filtered) > 1 and data_filtered['result_num'].nunique() != 1:
        continue
      same_data = pd.concat([same_data, group], ignore_index=True)

    ################# construct avg data ##################
    for size_ in querysizes:
      for filter_ in filters:
        for order in orders:
          for engine in engines:
            tmp = same_data[(same_data.filter_ == filter_)&(same_data.order == order)&(same_data.engine == engine)&(same_data.size_ == size_)]
            if tmp.shape[0] == 0:
                continue
            print(graph+", "+labelsize+", "+str(size_)+", "+filter_+", "+order+", "+engine+"\n")
            allrow = [graph, labelsize, size_, filter_, order, engine]
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
glsout.to_csv('ratio.csv', index=False, float_format='%.6f')
