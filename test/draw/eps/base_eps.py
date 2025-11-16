# here put the import lib
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import seaborn as sns

def draw_ax(file, ax):
    df = pd.read_csv(file)

    graphs = ["Figeys", "YeastS", "Citeseer", "HPRD", "Human", "WordNet", "DBLP", "Stanford", "Youtube", "Twitch"]
    methods = ["KSS", "RM", "GUP", "BICE", "VEQ", r"B$\circledS X$"]

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
    df = df[(df["method"].isin(methods)) &(df["graph"].isin(graphs))]

    # 计算每个graph中每个method的平均eps
    mean_eps = df.groupby(['graph', 'method'])['eps'].mean().unstack()
    mean_eps = mean_eps[methods].div(mean_eps["VEQ"], axis=0)
    # exit()
    mean_eps = np.log10(mean_eps)
    print(mean_eps)

    # 对柱子顺序进行排序
    # mean_eps = mean_eps.reindex(graphs)
    mean_eps = mean_eps.sort_values(by=r"B$\circledS X$")

    # 调色, 填充样式
    sns.set_style("whitegrid")
    palette = sns.color_palette("husl", len(methods), desat=0.6)
    hatches = ['//', '\\\\', '||', 'oo', '--', 'xx']

    # 计算每个method的柱状图位置
    bar_width = 0.14
    bar_positions = np.arange(len(mean_eps.index))

    # 绘制每个method的柱状图
    for i, method in enumerate(methods):
        ax.bar(bar_positions + i * bar_width, mean_eps[method], width=bar_width, color='none', hatch=hatches[i % len(hatches)], edgecolor=palette[i], label=method)

    # 添加标签和标题
    ax.set_xticks(bar_positions + bar_width * (len(mean_eps.columns) - 1) / 2)
    ax.set_xticklabels(mean_eps.index, rotation=0, fontsize=12)
    # ax.set_yscale('log')
    ax.grid(False)
    ax.set_xlim(-0.2, len(methods) + 4)

    return ax

def plot_broken(ax1,ax2):
    #绘制断裂处的标记
    # d = .85  #设置倾斜度    
    # kwargs = dict(marker=[(-1, -d), (1, d)], markersize=5,
    # linestyle='none', color='k', mec='k', mew=1, clip_on=False)
    # ax2.plot([0, 1], [0, 0],transform=ax2.transAxes, **kwargs)
    # ax1.plot([0, 1], [1, 1], transform=ax1.transAxes, **kwargs)
    ax2.spines['bottom'].set_visible(False)#关闭子图2中底部脊
    ax1.spines['top'].set_visible(False)##关闭子图1中顶部脊
    ax2.set_xticks([])

if __name__ == "__main__":
    fig = plt.figure(figsize = (10,3))
    ax1 = fig.add_axes([0.07,0.54,0.92,0.4])
    ax2 = fig.add_axes([0.07,0.14,0.92,0.4])
    draw_ax('../eps_old.csv', ax1)
    draw_ax('../eps_old.csv', ax2)

    # 设置 y轴 范围
    ax1.set_ylim(0,10)
    ax2.set_ylim(-87,0)

    # 设置 label
    ax1.set_yticks(np.arange(0,10,3))
    ax2.set_yticks(np.arange(-90,1,30))
    y1ticklabels = ['1', r'$10^{3}$', r'$10^{6}$', r'$10^{9}$']
    ax1.set_yticklabels(y1ticklabels, fontsize=8)
    y2ticklabels = [ r'$10^{-90}$', r'$10^{-60}$', r'$10^{-30}$', '1']
    ax2.set_yticklabels(y2ticklabels, fontsize=8)

    # 画截断点
    plot_broken(ax2,ax1)

    # 字体设置
    plt.rcParams.update({
        'font.size': 12,
        'axes.titlesize': 16,
        'axes.labelsize': 16,
        'legend.fontsize': 12,
        'xtick.labelsize': 16,
    })

    # 显示图表
    ax1.set_ylabel('Relative EPS Ratio                        ', fontsize=12, labelpad=10)#空格调节令ylabel居中
    ax1.legend(loc='upper center', frameon=False, ncol=6)
    plt.savefig('overall_base.pdf')
