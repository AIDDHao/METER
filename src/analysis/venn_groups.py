"""
自定义模块：绘制多组数据的交集venn图
"""

import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from venn import venn
import venn
import os
from matplotlib.colors import ListedColormap

# 两组数据绘制venn
def venn_2g(data1, data2, name_data, path_dir, title_name='OTTER-AI Target Intersections'):
    """
    data1,data2:2组数据
    name_data:各组的名称
    title_name='Intersection' ：venn图标题
    """

    os.makedirs(path_dir, exist_ok=True)
    # venn数据
    data1 = set(data1)
    data2 = set(data2)
    # 各组名称
    name1 = list(name_data)[0]
    name2 = list(name_data)[1]
    # 设置字体
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams['axes.unicode_minus'] = False
    # subsets参数
    subset = [[data1, data2]]
    for i in subset:
        my_dpi = 100
        plt.figure(figsize=(1800/my_dpi, 750/my_dpi), dpi=300)
        g = venn2(subsets=i,
                  set_labels=(name1, name2),
                  set_colors=("#1a4fb2", "#f8d90f"),
                  alpha=0.6,
                  normalize_to=1.0
                  )
        g.get_label_by_id('10').set_fontsize(20)
        g.get_label_by_id('11').set_fontsize(20)
        g.get_label_by_id('01').set_fontsize(20)
        for label in g.set_labels:
            label.set_fontsize(12)
            label.set_weight('bold')
            label.set_y(-0.6)
        plt.rcParams.update({'font.size': 20, 'font.weight': 'bold',
                             'axes.labelweight': 'bold',
                             'axes.titleweight': 'bold',
                             'xtick.labelsize': 'large',
                             'ytick.labelsize': 'large'
                             })
        plt.title(title_name, fontsize=20),
        plt.rcParams['font.sans-serif'] = ['Arial']
        plt.savefig(f'{path_dir}/{name1}_{name2}_OTTER2_Venn.svg', dpi=300)


# 3组数据venn图
def venn_3g(data1, data2, data3, name_data, path_dir, title_name='OTTER-AI Target Intersections'):
    """
    data1,data2,data3:3组数据
    name_data:各组的名称
    title_name='Intersection' ：venn图标题
    """
    os.makedirs(path_dir, exist_ok=True)
    # venn数据
    data1 = set(data1)
    data2 = set(data2)
    data3 = set(data3)
    # 各组名称
    name1 = list(name_data)[0]
    name2 = list(name_data)[1]
    name3 = list(name_data)[2]
    # 设置字体
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams['axes.unicode_minus'] = False
    # subsets参数
    subset = [[data1, data2, data3]]
    for i in subset:
        my_dpi = 100
        plt.figure(figsize=(1800 / my_dpi, 750 / my_dpi), dpi=300)
        g = venn3(subsets=i,
                  set_colors=("#509579", "#8d9fbb", "#d89c7c"),
                  alpha=0.6,
                  normalize_to=1.0,
                  )
        # 设置每部分的字体大小
        for i0 in list('100, 010, 110, 001, 101, 011, 111'.split(', ')):
            try:
                g.get_label_by_id('%s' % i0).set_fontsize(20)
            except (AttributeError, TypeError) as e:
                pass
        for label in g.set_labels:
            label.set_visible(False)

        plt.rcParams.update({'font.size': 20})
        plt.title(title_name, fontsize=24),
        plt.legend(list(name_data), loc='upper left', bbox_to_anchor=(0.95, 1), fontsize=16,

                   )

        plt.rcParams['font.sans-serif'] = ['Arial']
        plt.savefig(f'{path_dir}/{name1}_{name2}_{name3}_OTTER2_enn.svg', dpi=300)



# 4组数据venn图
def venn_4g(data1, data2, data3, data4, name_data, path_dir, title_name='OTTER-AI Target Counts', colors=None,
            output_name='OTTER2_Target_Counts_4_Groups.svg', dpi=300, alpha=0.6, font_size=20):
    """
    data1,data2, data3, data4:4组数据
    name_data:各组的名称
    title_name='Intersection' ：venn图标题
    """
    if len(name_data) != 4:
        raise ValueError("必须包含4个集合名。")
    # 颜色
    if colors is None:
        colors = ['#8fbedc', '#7fbbb2', '#cdc6e7', '#ebca8e']

    if len(colors) != 4:
        raise ValueError("colors必须包含4个颜色。")
    c_cmap = ListedColormap(colors)


    # 画venn图
    data1 = set(data1)
    data2 = set(data2)
    data3 = set(data3)
    data4 = set(data4)

    # 合并数据
    all_sets = [data1, data2, data3, data4]
    data_dict = {name: s for name, s in zip(name_data, all_sets)}

    # 设置字体
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams['axes.unicode_minus'] = False
    # 创建画布
    fig, ax = plt.subplots(figsize=(18, 7.5), dpi=dpi)
    legend_loc = 'upper left'
    venn.venn(
        data_dict,
        cmap=c_cmap,
        alpha=alpha,
        fontsize=font_size,
        ax=ax,
        legend_loc=legend_loc
    )

    # 标题
    ax.set_title(title_name, fontsize=font_size + 4)
    # 图例位置
    ax.legend(
        list(name_data),
        loc=legend_loc,
        bbox_to_anchor=(1.02, 1),
        fontsize=font_size - 4,
        borderaxespad=0
    )
    os.makedirs(path_dir, exist_ok=True)
    out_path = os.path.join(path_dir, output_name)
    fig.savefig(out_path, dpi=dpi, bbox_inches='tight')


# 5组数据venn图
def venn_5g(data1, data2, data3, data4, data5, name_data, path_dir, title_name='OTTER-AI Target Counts', colors=None,
            output_name='OTTER2_Target_Counts_5_Group.svg', dpi=300, alpha=0.6, font_size=20):
    """
    data1,data2, data3, data4, data5:5组数据
    name_data:各组的名称
    title_name='Intersection'：venn图标题
    colors=None:自定义颜色
    """

    if len(name_data) != 5:
        raise ValueError("必须包含5个集合名。")
    # 颜色
    if colors is None:
        colors = ['#FEDAA4', '#F6B7B4', '#ACD8C6', '#CBC355', '#B1C4D9']

    if len(colors) != 5:
        raise ValueError("colors必须包含5个颜色。")
    # 构建自定义颜色映射
    c_cmap = ListedColormap(colors)


    # 画venn图
    data1 = set(data1)
    data2 = set(data2)
    data3 = set(data3)
    data4 = set(data4)
    data5 = set(data5)

    # 合并所有数据
    all_sets = [data1, data2, data3, data4, data5]
    data_dict = {name: s for name, s in zip(name_data, all_sets)}
    # 设置字体
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams['axes.unicode_minus'] = False
    # 创建画布
    fig, ax = plt.subplots(figsize=(18, 7.5), dpi=dpi)
    legend_loc = 'upper left'
    venn.venn(
        data_dict,
        cmap=c_cmap,
        alpha=alpha,
        fontsize=font_size,
        ax=ax,
        legend_loc=legend_loc
    )

    # 标题
    ax.set_title(title_name, fontsize=font_size + 4)
    # 图例位置
    ax.legend(
        list(name_data),
        loc=legend_loc,
        bbox_to_anchor=(1.02, 1),
        fontsize=font_size - 4,
        borderaxespad=0
    )
    os.makedirs(path_dir, exist_ok=True)
    out_path = os.path.join(path_dir, output_name)
    fig.savefig(out_path, dpi=dpi, bbox_inches='tight')



if __name__ == '__main__':
    list1 = [1, 2, 3, 4, 5]
    list2 = (1, 2, 5, 6, 7)
    list3 = [1, 3, 5, 9, 8, 10, 6, 4]
    list4 = [1, 3, 5, 9, 7, 10, 6, 11]
    list5 = [1, 2, 6, 5, 11, 10]
    list6 = [1, 2, 6, 5, 3, 12, 9]
    venn_2g(list1, list2, ['list1', 'list2'], "test_venn")
    venn_3g(list1, list2, list3, ['list1', 'list2', 'list3'], "test_venn")
    venn_4g(list1, list2, list3, list4, ['list1', 'list2', 'list3', 'list4'], "test_venn")
    venn_5g(list1, list2, list3, list4, list5, ['list1', 'list2', 'list3', 'list4', 'list5'], "test_venn")


