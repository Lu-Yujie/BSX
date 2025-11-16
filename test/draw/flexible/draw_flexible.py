import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 读取数据
df_similar = pd.read_csv('similar.csv')
df_refine = pd.read_csv('refine.csv')

graphs = ["Figeys", "YeastS", "Citeseer", "HPRD", "Human", "WordNet", "DBLP", "Stanford", "Youtube", "Twitch"]
markers = ['^', '*', 's', 'o', 'v', 'D']
labels = ["L15", "L30", "L45", "L60"]
sims = ["s0", "s25", "s50", "s75", "s100"]
refines = ["r0", "r25", "r50", "r75", "r100"]
x_labels = [0, 0.25, 0.5, 0.75, 1]

# 替换 graph 名称
df_similar['graph'] = df_similar['graph'].replace({
    'web-Stanford': 'Stanford',
    'maayan-figeys': 'Figeys',
    'wordnet-words': 'WordNet',
    'human': 'Human',
    'dblp': 'DBLP',
    'twitch': 'Twitch',
    'youtube': 'Youtube',
    'citeseer': 'Citeseer'
})
df_refine['graph'] = df_refine['graph'].replace({
    'web-Stanford': 'Stanford',
    'maayan-figeys': 'Figeys',
    'wordnet-words': 'WordNet',
    'human': 'Human',
    'dblp': 'DBLP',
    'twitch': 'Twitch',
    'youtube': 'Youtube',
    'citeseer': 'Citeseer'
})

# 筛选所需数据
df_similar = df_similar[(df_similar["graph"].isin(graphs)) & (df_similar['label'].isin(labels))]
df_refine = df_refine[(df_refine["graph"].isin(graphs)) & (df_refine['label'].isin(labels))]

# 设置调色板
palette = sns.color_palette("husl", len(graphs), desat=0.6)

# 初始化图形
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# 绘制每个 graph 的折线 - similar.csv
for idx, graph in enumerate(graphs):
    # similar.csv 处理
    df_graph_similar = df_similar[df_similar['graph'] == graph]
    df_s100 = df_graph_similar[df_graph_similar['sim'] == 's100']['eps'].mean()
    df_graph_similar['relative_eps'] = (df_graph_similar['eps'] - df_s100) / df_s100
    df_label_similar = df_graph_similar.groupby('sim')['relative_eps'].mean()
    df_label_similar = df_label_similar.reindex(sims)
    
    # 绘制 similar 图表
    ax1.plot(sims, df_label_similar.values, marker=markers[idx % len(markers)], markersize=8, 
             linewidth=2, label=graph, color=palette[idx])

# 设置 similar 图的参数
ax1.set_xlabel(r'Similarity Threshold($\tau$)', fontsize=22)
ax1.set_ylabel('Relative EPS Difference    ', fontsize=22)
ax1.set_xticklabels(x_labels, fontsize=18)
ax1.tick_params(axis='y', labelsize=18)
ax1.grid(True, linestyle='--', alpha=0.6)

# 绘制每个 graph 的折线 - refine.csv
for idx, graph in enumerate(graphs):
    # refine.csv 处理
    df_graph_refine = df_refine[df_refine['graph'] == graph]
    df_r0 = df_graph_refine[df_graph_refine['refine'] == 'r0']['eps'].mean()
    df_graph_refine['relative_eps'] = (df_graph_refine['eps'] - df_r0) / df_r0
    df_label_refine = df_graph_refine.groupby('refine')['relative_eps'].mean()
    df_label_refine = df_label_refine.reindex(refines)
    
    # 绘制 refine 图表
    ax2.plot(refines, df_label_refine.values, marker=markers[idx % len(markers)], markersize=8, 
             linewidth=2, label=graph, color=palette[idx])

# 设置 refine 图的参数
ax2.set_xlabel(r'Refinement Threshold($\mathbf{r}$)', fontsize=22)
ax2.set_xticklabels(x_labels, fontsize=18)
ax2.tick_params(axis='y', labelsize=18)
ax2.grid(True, linestyle='--', alpha=0.6)

# 共享图例
handles, labels = ax1.get_legend_handles_labels()
fig.legend(handles, labels, frameon=False, loc='upper center', fontsize=20, ncol=5)

# 调整布局
plt.tight_layout(rect=[0.005, 0, 1, 0.84])
plt.subplots_adjust(wspace=0.2)
# 保存图像
plt.savefig('flexible.pdf')

