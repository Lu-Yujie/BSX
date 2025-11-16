import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# 设定需要绘制的 graphs 和 methods
graphs = ["WordNet", "Youtube"]
methods = ["KSS", "GUP", "BICE", "VEQ", r"B$\circledS X$"]
method_numbers = {
    "KSS": 0,
    "GUP": 2,
    "BICE": 3,
    "VEQ": 4,
    r"B$\circledS X$": 5,
}
graph_num = 8  # 每个图的 query 数量，假设是 8
o_limit = 'o_10_5'

df = pd.read_csv('global_'+o_limit+'.csv')
# 使用.replace()方法进行数据替换
df['method'] = df['method'].replace({
    'KSS': 'KSS',
    'BSX': r'B$\circledS X$',
    'VEQ': 'VEQ',
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

fig = plt.figure(figsize=(15, 5))
ax = plt.gca()
sns.set_theme(style="whitegrid")

# 调色, 填充样式
palette = sns.color_palette("husl", 6, desat=0.6)
hatches = ['//', '\\\\', '||', 'oo', '--', 'xx']

total_queries = len(graphs) * graph_num
x = np.arange(total_queries)
width = 0.14  # 每个柱子的宽度

# 绘制每个 method 的数据，并分割不同的 graphs
for idx, graph in enumerate(graphs):
    # 过滤数据
    df_graph = df[(df["method"].isin(methods)) & (df["graph"] == graph)]

    # 计算每个 graph 中每个 method 的平均 time
    df_grouped = df_graph.groupby(['graph_id', 'method'])['time'].mean().unstack()

    # 绘制每个 method 的柱状图
    for i, method in enumerate(methods):
        values = df_grouped[method]
        method_number = method_numbers[method]
        ax.bar(x[graph_num*idx:graph_num*(idx+1)] + i * width, values + 0.01, width=width,
                color='none', hatch=hatches[method_number], edgecolor=palette[method_number], label=method if idx == 0 else "")

    # 在每个 graph 后添加分割线
    ax.axvline(x=graph_num * (idx + 1) - 0.22, color='gray', linestyle='--')

    # 在对应的区域上方添加 graph 标签
    mid_position = graph_num * idx + (graph_num / 2)  # 找到每个 graph 的中点位置
    ax.text(mid_position, -0.1 * ax.get_ylim()[0]-0.12, graph, ha='center', va='top', fontsize=28, transform=ax.get_xaxis_transform())


# 设置图形参数
ax.set_yscale('log')  # 使用对数刻度
ax.set_ylabel('QPT (s)', fontsize=24)
ax.set_xlabel('')
ax.tick_params(axis='y', labelsize=22)  # 修改 y 轴刻度的字体大小为 16

# 设置 x 轴标签为不同的 graph 名称和 query 编号
queries_per_graph = [f"q{i+1}" for i in range(graph_num)]
xtick_labels = []
for graph in graphs:
    xtick_labels += queries_per_graph

ax.set_xticks(np.arange(total_queries) + width * len(methods) / 2)
ax.set_xticklabels(xtick_labels, fontsize=22)
ax.set_xlim(-0.15, total_queries - 0.25)
# ax.set_ylim(0.01, 0.0112)
# ax.set_yticks(np.arange(0.01, 0.0113, 0.001))
ax.grid(True, which='major', linestyle='--', color='gray', alpha=0.3)  # 设置网格为灰色虚线

# 调整布局
plt.tight_layout()

plt.legend(frameon=False, fontsize=26, loc='upper center', ncol=6, bbox_to_anchor=(0.5, 1.19))
plt.minorticks_off()
plt.subplots_adjust(top=0.9, bottom=0.16, left=0.09)

# 保存图像
plt.savefig('global_'+o_limit+'.pdf')
