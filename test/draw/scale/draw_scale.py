import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# 创建数据
df = pd.read_csv('scale.csv')
graphs = ["Figeys", "Citeseer"]
scales = [5, 10, 20, 50, 100]
methods = ["KSS", "RM", "GUP", "BICE", "VEQ", r"B$\circledS X$"]

# 使用.replace()方法进行替换
df['method'] = df['method'].replace({
    'DPiso-DPiso-KSS': 'KSS',
    'RM-RM-RM': 'RM',
    'CFL-nan-BSX': r'B$\circledS X$',
    'VEQ-VEQ-VEQ': 'VEQ',
    "bice": "BICE",
    "gup": "GUP",
})
df['graph'] = df['graph'].replace({
    'web-Stanford': 'Stanford',
    'maayan-figeys': 'Figeys',
    'wordnet-words': 'WordNet',
    'human': 'Human',
    'dblp': 'DBLP',
    'twitch': 'Twitch',
    'youtube': 'Youtube',
    'citeseer': 'Citeseer'
})
# 调色, 填充样式
sns.set_style("whitegrid")
palette = sns.color_palette("husl", len(methods), desat=0.6)
hatches = ['//', '\\\\', '||', 'oo', '--', 'xx']
full_df = df
plt.figure(figsize=(12, 4))
for idx, graph in enumerate(graphs):
  plt.subplot(1, 2, idx + 1)
  # filter all data
  df = full_df[(full_df["method"].isin(methods)) &(full_df["scale"].isin(scales)) &(full_df["graph"]==graph)]
  # 计算每个graph中每个method的平均eps
  df = df.groupby(['scale', 'method'])['eps'].mean().unstack()
  print(df)

  # 设置绘图参数
  x = np.arange(len(scales))  # x 轴标签的位置
  width = 0.14  # 柱状图的宽度

  for i, method in enumerate(methods):
      values = df[method]
      plt.bar(x + i * width, values, width, color='none', hatch=hatches[i % len(hatches)], edgecolor=palette[i], label=method)

  plt.xlabel(graph, fontsize=24)
  plt.xlim(-0.15, len(methods) - 1.15)
  if (idx == 0):
    plt.ylabel('EPS', fontsize=24)
  plt.yscale('log')  # 根据数据范围调整刻度类型
  plt.grid(True, which='major', linestyle='--', color='gray', alpha=0.3)  # 设置网格为灰色虚线
  # 调整 x 轴和 y 轴 ticks 的字体大小
  plt.xticks(x+0.35, scales, fontsize=18)
  plt.yticks(fontsize=18)
  plt.rcParams['axes.linewidth'] = 1
  for spine in plt.gca().spines.values():
    spine.set_color('black')
  plt.tight_layout()
plt.legend(loc='upper center',frameon=False, fontsize=20, ncol=6, bbox_to_anchor=(-0.1,1.2))
plt.subplots_adjust(top=0.9)
plt.subplots_adjust(wspace=0.15)
plt.savefig('scale.pdf')
