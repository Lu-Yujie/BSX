import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import seaborn as sns

# 创建数据
df = pd.read_csv('../eps.csv')
# graphs = ["Figeys", "YeastS", "Citeseer", "HPRD", "Human", "WordNet", "DBLP", "Stanford", "Youtube", "Twitch"]
# methods = ["gup", "bice", "DPiso-DPiso-KSS", "RM-RM-RM", "DPiso-nan-BSX", "NLF-nan-BSX", "CFL-nan-BSX", "VEQ-VEQ-VEQ"]
graphs = ["Figeys", "YeastS", "Citeseer", "HPRD", "Human", "WordNet", "DBLP", "Stanford", "Youtube", "Twitch"]
methods = [r"CFL-B$\circledS X$", r"DPiso-B$\circledS X$", r"NLF-B$\circledS X$"]

# 使用.replace()方法进行替换
df['method'] = df['method'].replace({
    'CFL-nan-BSX': r'CFL-B$\circledS X$',
    "DPiso-nan-BSX": r'DPiso-B$\circledS X$',
    "NLF-nan-BSX": r'NLF-B$\circledS X$',
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
hatches = ['xx', 'oo', '--']

# filter all data
df = df[(df["method"].isin(methods)) &(df["graph"].isin(graphs))]
# 计算每个graph中每个method的平均eps
df = df.groupby(['graph', 'method'])['eps'].mean().unstack().reset_index()

# 将数据按照 graphs 的顺序排列
sorted_df = df[df["graph"]==graphs[0]]
for i in range(1,len(graphs)):
    sorted_df = pd.concat([sorted_df, df[df["graph"]==graphs[i]]])
df = sorted_df
df[methods] = df[methods].div(df[r"DPiso-B$\circledS X$"], axis=0)
print(df)
df[r'NLF-B$\circledS X$'] = np.log10(df[r'NLF-B$\circledS X$'])
df[r"CFL-B$\circledS X$"] = np.log10(df[r"CFL-B$\circledS X$"])
df[r"DPiso-B$\circledS X$"] = 0
print(df)

# 设置柱状图的宽度和柱子的位置
bar_width = 0.3
index = np.arange(len(graphs))

def plot_ax(ax):
    # 绘制每个 graph 的柱状图
    for i, method in enumerate(methods):
        values = df[method]
        ax.bar(index + bar_width*i, values, bar_width, color="none",
                hatch=hatches[i % len(hatches)], edgecolor=palette[i], label=method,
                )
    for i in range(len(graphs)-1):
        ax.axvline(x=(i + 1) - 0.2, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    # 在每个 graph 后添加分割线
    return ax

fig = plt.figure(figsize = (8,4))
ax2 = fig.add_axes([0.12,0.14,0.87,0.4])
ax1 = fig.add_axes([0.12,0.54,0.87,0.4])

plot_ax( ax1)
plot_ax( ax2)
# 设置 y轴 范围
ax1.set_ylim(0,0.4)
ax2.set_ylim(-9,0)

# 设置 label
ax1.set_yticks(np.arange(0,9.1,3))
ax2.set_yticks(np.arange(-13,0.1,4))
y1ticklabels = ['1', r'$10^{3}$', r'$10^{6}$', r'$10^{9}$']
ax1.set_yticklabels(y1ticklabels, fontsize=10)
y2ticklabels = [ r'$10^{-13}$', r'$10^{-9}$', r'$10^{-5}$', '']
ax2.set_yticklabels(y2ticklabels, fontsize=10)

# 添加标签和标题
plt.xlabel('', fontsize=20)
plt.ylabel('Relative EPS Ratio                ', fontsize=15.5,labelpad=6)
ax2.set_xticks(index + 1.5 * bar_width - 0.12, graphs, fontsize=12)
ax1.set_xticks([])
plt.legend(loc='upper center', frameon=False, ncol=3, fontsize=13, bbox_to_anchor=(0.5, 1.05))
ax1.grid(False)#, which='both', linestyle='--', color='gray', alpha=0.3)  # 设置网格为灰色虚线
ax2.grid(False)#, which='both', linestyle='--', color='gray', alpha=0.3)  # 设置网格为灰色虚线

ax1.spines['top'].set_linewidth(0.9)
ax1.spines['top'].set_color('black')
ax1.spines['left'].set_linewidth(0.9)
ax1.spines['left'].set_color('black')
ax1.spines['right'].set_linewidth(0.9)
ax1.spines['right'].set_color('black')

ax2.spines['bottom'].set_linewidth(0.9)
ax2.spines['bottom'].set_color('black')
ax2.spines['left'].set_linewidth(0.9)
ax2.spines['left'].set_color('black')
ax2.spines['right'].set_linewidth(0.9)
ax2.spines['right'].set_color('black')

ax1.set_xlim(-0.25, len(graphs) - 0.1)
ax2.set_xlim(-0.25, len(graphs) - 0.1)
# 显示图形
plt.savefig('filter_base.pdf')
