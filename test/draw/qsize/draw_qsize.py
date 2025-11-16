import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# 创建数据
df = pd.read_csv('../eps.csv')
graphs = ["Figeys", "WordNet", "Twitch"]
methods = ["KSS", "RM", "GUP", "BICE", "VEQ", r"B$\circledS X$"]
markers = ['^', '*', 's', 'o', 'v', 'D']
qsizes = [10, 20, 30, 40, 50]
qsizesLabel = ["Q10", "Q20", "Q30", "Q40", "Q50"]

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
df['label'] = df['label'].replace({
    '15': 'L15',
    '30': 'L30',
    '45': 'L45',
    '60': 'L60',
})

# filter all data
df = df[(df["method"].isin(methods)) &(df["graph"].isin(graphs)) & (df['size'].isin(qsizes))]

palette = sns.color_palette("husl", len(methods), desat=0.6)

# 计算每个graph中每个method的平均eps
df = df.groupby(['graph', 'size', 'method'])['eps'].mean().unstack()
print(df)
plt.figure(figsize=(12, 4))
for idx, graph in enumerate(graphs):
    plt.subplot(1, len(graphs), idx + 1)
    df_graph = df.loc[graph]
    for method, marker, color in zip(methods, markers, palette):
        df_method = df_graph[method]
        print(df_method)
        plt.plot(df_method.index, df_method.values+0.1, marker=marker, markersize=8, color=color, linewidth=2.5, label=method)
    plt.xlabel(graph, fontsize=24)
    plt.xticks(qsizes, qsizesLabel, fontsize=26)
    plt.yticks(fontsize=15)
    if idx == 0:
        plt.ylabel('EPS', fontsize=24)
    plt.yscale('log')  # 根据数据范围调整刻度类型
    plt.grid(True, which='major', linestyle='--', color='gray', alpha=0.6)  # 设置网格为灰色虚线
    # 调整 x 轴和 y 轴 ticks 的字体大小
    plt.xticks(fontsize=16)
    plt.tight_layout()
# 共用一个 legend
plt.legend(frameon=False, fontsize=20, loc='upper center', ncol=6, bbox_to_anchor=(-0.8, 1.2))
plt.subplots_adjust(top=0.9)
plt.subplots_adjust(wspace=0.3)
plt.savefig('qsize.pdf')  # 保存图像

