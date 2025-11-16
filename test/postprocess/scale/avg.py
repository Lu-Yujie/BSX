# combine output limit exp data

import pandas as pd

data_dir = '/var/lib/docker/subgraph/output/bsx/synthetic/EvoGraph'

graphs = ["human", "HPRD", "citeseer", "YeastS", "Figeys"]
sacles = ["5", "10", "20", "50", "100"]
qsizes = [10, 20, 30, 40, 50]
avgout = pd.DataFrame(columns=['graph', 'scale', 'qsize', 'filter_', 'order', 'engine', 'time', 'num', 'eps'])
avgout_index = 0
for graph in graphs:
  for sacle in sacles:
    avgtmp = pd.read_csv(data_dir + '/' + graph + '/' + graph + '_' + sacle + '/result.csv', na_values=['NULL'], keep_default_na=False)
    filters = avgtmp['filter_'].drop_duplicates().tolist()
    orders = avgtmp['order'].drop_duplicates().tolist()
    engines = avgtmp['engine'].drop_duplicates().tolist()
    for qsize in qsizes:
      for filter_ in filters:
        for order in orders:
          for engine in engines:
            avgdata = avgtmp[(avgtmp.filter_ == filter_)&(avgtmp.order == order)&(avgtmp.engine == engine)\
                            &(avgtmp.size_ == qsize)]
            if avgdata.shape[0] == 0:
                continue
            results = []
            row = [graph, sacle, qsize, filter_, order, engine]
            avgdata = avgdata[['total_time', 'result_num']]
            for column in avgdata.columns:
              avgdata[column] = pd.to_numeric(avgdata[column], errors='coerce')
              results.append(sum(avgdata[column]))
            results.append(results[1]/(results[0]))
            row += results
            avgout.loc[avgout_index] = row
            avgout_index += 1
avgout.to_csv('syn.csv',index=False)
