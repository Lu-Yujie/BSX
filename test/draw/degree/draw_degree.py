import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.ticker import MaxNLocator

# 创建数据
df = pd.read_csv('degree.csv')
graphs = ["Figeys", "Citeseer", "Twitch"]
methods = ["KSS", "RM", "GUP", "BICE", "VEQ", r"B$\circledS X$"]
degrees = [1, 2, 3, 4, 5]
degreesLabel = ["[1,2)", "[2,3)", "[3,4)", "[4,5)", "[5,6)"]

# 使用.replace()方法进行替换
df['method'] = df['method'].replace({
    'DPiso-DPiso-KSS': 'KSS',
    'RM-RM-RM': 'RM',
    'CFL-nan-BSX': r"B$\circledS X$",
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

# filter all data
df = df[(df["method"].isin(methods)) &(df["graph"].isin(graphs)) & (df['degree'].isin(degrees))]

# 计算每个graph中每个method的平均eps
df = df.groupby(['graph', 'degree', 'method'])['eps'].mean().unstack()

# 设置绘图参数
x = np.arange(len(degrees))  # x 轴标签的位置
width = 0.14  # 柱状图的宽度

# 调色, 填充样式
sns.set_style("whitegrid")
palette = sns.color_palette("husl", len(methods), desat=0.6)
hatches = ['//', '\\\\', '||', 'oo', '--', 'xx']

plt.figure(figsize=(12, 4))
# 绘制每个查询图下的结果
for idx, graph in enumerate(graphs):
  plt.subplot(1, len(graphs), idx + 1)
  # 绘制每个方法的柱状图
  g_data = df.loc[graph]
  print(g_data)
  for i, method in enumerate(methods):
    values = g_data[method]
    plt.bar(x + i * width, values, width, color='none', hatch=hatches[i % len(hatches)], edgecolor=palette[i], label=method)
  plt.xlabel(graph, fontsize=24)
  plt.yticks(fontsize=18)
  if idx == 0:
    plt.ylabel('EPS', fontsize=24)
  plt.yscale('log')
  plt.grid(False)
  plt.xlim(-0.2, len(methods)-1.1)
  plt.xticks(x+0.35, degreesLabel, fontsize=18)
  plt.rcParams['axes.linewidth'] = 1
  for spine in plt.gca().spines.values():
    spine.set_color('black')
  # plt.gca().yaxis.set_major_locator(MaxNLocator(nbins=2))
plt.tight_layout()
# 共用一个 legend
plt.legend(frameon=False, fontsize=20, loc='upper center', ncol=6, bbox_to_anchor=(-0.8, 1.2))
plt.subplots_adjust(wspace=0.25)
plt.subplots_adjust(top=0.9)
plt.savefig('degree.pdf')  # 保存图像
