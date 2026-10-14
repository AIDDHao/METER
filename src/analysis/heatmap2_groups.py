"""
自定义模块：绘制多组数据的热图
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm, Normalize

def create_three_color_colormap(colors, max_score):
    """
    使用三种指定颜色创建线性渐变颜色映射，并根据max_score归一化
    规则：
    0 ~ 1：使用第一段动态变化
    1 ~ max_score：使用第二段动态变化
    """
    n_bins = 256
    three_color_cmap = LinearSegmentedColormap.from_list('three_color', colors, N=n_bins)

    if max_score <= 1:
        norm = Normalize(vmin=0, vmax=1)
    else:
        norm = TwoSlopeNorm(vmin=0, vcenter=1, vmax=max_score)

    return three_color_cmap, norm

# 定义函数
def heatmap_scores(path_csv, col_data, index_data, path_dir, output_fn = " .png",
                   len_num = 50, min_number=0, max_number=2):
    """
    path_csv,
    header_data:标题
    data_l:data list
    col_data:columns=各列的名称
    index_data:index=gene列表
    path_dir:要保存的路径文件夹
    output_fn = " .png":输出的文件夹名，默认为空
    len_num:选择用多少个基因绘制热图
    """

    # 定义输出的图片文件名
    if output_fn == " .png":
        output_fn = "Heatmap_Vlag_Systemnum.png"
    else:
        pass

    path_o = os.path.join(path_dir, output_fn)

    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.weight'] = 'bold'

    # 读取文件获取数据
    # 读取数据文件
    data_read = pd.read_csv(path_csv)

    if len(index_data) > len_num:
        index_data = index_data[: len_num]

    else:
        index_data = index_data

    # 提取数据
    data = data_read.iloc[: len(index_data), 1:]
    gene_name_list = data_read.iloc[: len(index_data), 0].tolist()
    data_array = np.array(data)
    data_list = data_array.tolist()
    df = pd.DataFrame(data_list,
                      columns=col_data,
                      index=[gene_name_list])

    # 制作热图
    if len(index_data) <= 20:
        fig, ax = plt.subplots(figsize=(10, 10))
    elif len(index_data) <= 40:
        fig, ax = plt.subplots(figsize=(10, 14))
    elif len(index_data) <= 60:
        fig, ax = plt.subplots(figsize=(10, 16))
    elif len(index_data) <= 80:
        fig, ax = plt.subplots(figsize=(10, 20))
    else:
        fig, ax = plt.subplots(figsize=(10, 25))
    sns.set(font_scale=1.2)

    h_color = "vlag"
    # 字体大小
    if len(index_data) <= 20:
        data_size = 16
    elif len(index_data) <= 30:
        data_size = 14
    elif len(index_data) <= 50:
        data_size = 12
    else:
        data_size = 10

    fig.subplots_adjust(left=0.18, right=0.82, top=0.92, bottom=0.12)
    # 颜色棒
    cbar_ax = fig.add_axes((0.88, 0.15, 0.02, 0.22))
    heatmap = sns.heatmap(df, annot=True, cmap=h_color, fmt='.2f', vmin=min_number, vmax=max_number,
                          linewidths=0.004, linecolor='white',
                          annot_kws={"size": data_size},
                          ax=ax,
                          cbar_ax=cbar_ax
                          )
    cbar = heatmap.collections[0].colorbar
    mid_number = (min_number + max_number) / 2
    cbar.set_ticks([min_number, mid_number, max_number])
    cbar.set_ticklabels([f"{min_number:g}", f"{mid_number:g}", f"{max_number:g}"])
    cbar.ax.tick_params(labelsize=max(data_size - 2, 8))

    # 标签
    cbar.ax.set_title("Score", fontsize=max(data_size - 2, 8), pad=10, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")

    # 字体大小
    if len(index_data) <= 20:
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
    elif len(index_data) <= 50:
        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=12)
    elif len(index_data) <= 80:
        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=12)
    else:
        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=10)

    # 设置x轴标签
    x_labels = []
    for col in df.columns:
        if isinstance(col, tuple):
            col = col[0]
        x_labels.append(str(col).replace(" ", "\n"))
    ax.set_xticklabels(x_labels)
    for label in ax.get_xticklabels():
        label.set_rotation(0)
        label.set_ha('center')
        label.set_va('top')
        label.set_multialignment('center')

    # 标题
    ax.set_title('OTTER-AI Total Scores with Ranking', fontsize=20,
                 fontdict={'family': 'Arial', 'weight': 'bold'},
              pad=10)
    plt.savefig(path_o, dpi=800)
    plt.close()



def heatmap_groups(path_csv, col_data, index_data, path_dir, output_fn = " .png",
                   len_num = 50, min_number=0, max_number=6):
    """
    path_csv,
    header_data:标题
    data_l:data list
    col_data:columns=各列的名称
    index_data:index=gene列表
    path_dir:要保存的路径文件夹
    output_fn = " .png":输出的文件夹名，默认为空
    len_num:选择用多少个基因绘制热图
    """

    if output_fn == " .png":
        output_fn = "Heatmap_Vlag_Systemnum.png"
    else:
        pass

    path_o = os.path.join(path_dir, output_fn)
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.weight'] = 'bold'

    # 读取数据
    data_read = pd.read_csv(path_csv)
    if len(index_data) > len_num:
        index_data = index_data[: len_num]

    else:
        index_data = index_data


    data = data_read.iloc[: len(index_data), 2:]
    gene_name_list = data_read.iloc[: len(index_data), 1].tolist()
    data_array = np.array(data)
    data_list = data_array.tolist()

    df = pd.DataFrame(data_list,
                      columns=col_data,
                      index=[gene_name_list])

    # 制作热图
    if len(index_data) <= 20:
        fig, ax = plt.subplots(figsize=(10, 10))
    elif len(index_data) <= 40:
        fig, ax = plt.subplots(figsize=(10, 14))
    elif len(index_data) <= 60:
        fig, ax = plt.subplots(figsize=(10, 16))
    elif len(index_data) <= 80:
        fig, ax = plt.subplots(figsize=(10, 20))
    else:
        fig, ax = plt.subplots(figsize=(10, 25))
    sns.set(font_scale=1.2)


    # 根据基因数量调整数据字体大小
    if len(index_data) <= 20:
        data_size = 16
    elif len(index_data) <= 30:
        data_size = 14
    elif len(index_data) <= 50:
        data_size = 12
    else:
        data_size = 10

    # 自定义颜色列表
    colors = ["#eef0f3", "#72a0c1", "#2076a9"]  # 用这个 蓝色系
    df_raw = df.copy()
    df_color = df.copy().astype(float)

    # 根据列数确定理论总分最大值
    # 最后一列为总分，前面每列单项得分理论最大值为1
    max_score = len(df.columns) - 1
    if max_score < 1:
        max_score = 1

    cmap, norm = create_three_color_colormap(colors, max_score)

    # 注释
    annot_df = df_raw.copy()
    for col in annot_df.columns:
        annot_df[col] = annot_df[col].map(lambda x: f"{x:.2f}")

    # cbar
    fig.subplots_adjust(left=0.18, right=0.82, top=0.92, bottom=0.12)
    cbar_ax = fig.add_axes((0.88, 0.15, 0.02, 0.22))

    heatmap = sns.heatmap(df_color, annot=annot_df, cmap=cmap, fmt='',
                          norm=norm,
                          linewidths=0.4, linecolor='white',
                          annot_kws={"size": data_size},
                          ax=ax,
                          cbar_ax=cbar_ax
                          )

    cbar = heatmap.collections[0].colorbar
    if max_score <= 1:
        cbar.set_ticks([0, 1])
        cbar.set_ticklabels(["0", "1"])
    else:
        cbar.set_ticks([0, max_score])
        cbar.set_ticklabels([f"{0:g}",f"{max_score:g}"])

    cbar.ax.tick_params(labelsize=max(data_size - 2, 8))
    for label in cbar.ax.get_yticklabels():
        label.set_fontweight('bold')

    cbar.ax.set_title("Score",
                      fontsize=max(data_size - 2, 8),
                      pad=10,
                      fontweight='bold')

    # 标签为空
    ax.set_xlabel("")
    ax.set_ylabel("")

    # 字体大小
    if len(index_data) <= 20:
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
    elif len(index_data) <= 50:
        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=12)
    elif len(index_data) <= 80:
        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=12)
    else:
        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=10)

    x_labels = []
    for col in df_color.columns:
        if isinstance(col, tuple):
            col = col[0]
        x_labels.append(str(col).replace(" ", "\n"))
    ax.set_xticklabels(x_labels)
    for label in ax.get_xticklabels():
        label.set_rotation(0)
        label.set_ha('center')
        label.set_va('top')
        label.set_multialignment('center')


    # 标题
    ax.set_title('OTTER-AI Keywords Scores with Ranking', fontsize=20,
              fontdict={'family': 'Arial', 'weight': 'bold'},
              pad=10)
    plt.savefig(path_o, dpi=800)
    plt.close()


def create_three_color_colormap_o(colors):
    """
    使用三种指定颜色创建线性渐变颜色映射
    """

    n_bins = 256
    three_color_cmap = LinearSegmentedColormap.from_list('three_color', colors, N=n_bins)

    return three_color_cmap

def heatmap_ppis_kw(path_data, path_dir, output_fn=" .png",
                   len_num=50, min_number=0, max_number=1):

    if output_fn == " .png":
        output_fn = "Heatmap_PPI_Mining_Score.png"
    else:
        pass

    path_o = os.path.join(path_dir, output_fn)

    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.weight'] = 'bold'

    # 读取数据

    if isinstance(path_data, pd.DataFrame):
        data_read = path_data.copy()
    elif str(path_data).endswith(".csv"):
        data_read = pd.read_csv(path_data)
    else:
        data_read = pd.read_excel(path_data, sheet_name="rank")


    cols = [col for col in data_read.columns.tolist()[1:] if col != "PM_Num"]
    if len(cols) == 0:
        raise ValueError("除Symbol_Name和PM_Num外，没有可用于绘制热图的列。")

    # 提取数据
    gene_col = data_read.iloc[:, 0]
    gene_col_list = gene_col.tolist()
    if len(gene_col_list) > len_num:
        gene_data = gene_col_list[: len_num]
    else:
        gene_data = gene_col_list
    data = data_read.loc[:len(gene_data) - 1, cols]

    gene_name_list = data_read.iloc[: len(gene_data), 0].tolist()
    data_array = np.array(data)
    data_list = data_array.tolist()
    df = pd.DataFrame(data_list,
                      columns=cols,
                      index=gene_name_list)

    # 制作热图
    if len(gene_data) <= 20:
        fig, ax = plt.subplots(figsize=(10, 10))
    elif len(gene_data) <= 40:
        fig, ax = plt.subplots(figsize=(10, 14))
    elif len(gene_data) <= 60:
        fig, ax = plt.subplots(figsize=(10, 16))
    elif len(gene_data) <= 80:
        fig, ax = plt.subplots(figsize=(10, 20))
    else:
        fig, ax = plt.subplots(figsize=(10, 25))
    sns.set(font_scale=1.2)

    # 自定义颜色
    colors = ["#EDF6F9", "#D7D5E3", "#FAD2DC"]
    cmap = create_three_color_colormap_o(colors)
    df_raw = df.copy()

    df_color = df.copy().astype(float)

    for col in df_color.columns:
        col_min = df_color[col].min()
        col_max = df_color[col].max()
        if col_max > col_min:
            df_color[col] = (df_color[col] - col_min) / (col_max - col_min)
        else:
            df_color[col] = 0.5

    # 生成注释
    annot_df = df_raw.copy()
    for col in annot_df.columns:
        if col in ["PPI_Total", "PPI_Report"]:
            annot_df[col] = annot_df[col].map(lambda x: f"{int(x)}")
        elif col == "Score":
            annot_df[col] = annot_df[col].map(lambda x: f"{x:.2f}")
        else:
            annot_df[col] = annot_df[col].map(lambda x: f"{x}")

    # 字体大小
    if len(gene_data) <= 20:
        data_size = 16
    elif len(gene_data) <= 30:
        data_size = 14
    elif len(gene_data) <= 50:
        data_size = 12
    else:
        data_size = 10
    # cbar
    fig.subplots_adjust(left=0.18, right=0.82, top=0.92, bottom=0.12)
    cbar_ax = fig.add_axes((0.88, 0.15, 0.02, 0.22))


    heatmap = sns.heatmap(df_color, annot=annot_df, cmap=cmap, fmt='',
                          vmin=min_number, vmax=max_number,
                          linewidths=0.4, linecolor='white',
                          annot_kws={"size": data_size},
                          ax=ax,
                          cbar_ax=cbar_ax
                          )

    # 设置color bar
    cbar = heatmap.collections[0].colorbar
    mid_number = (min_number + max_number) / 2
    cbar.set_ticks([min_number, mid_number, max_number])
    cbar.set_ticklabels([f"{min_number:g}", f"{mid_number:g}", f"{max_number:g}"])
    cbar.ax.tick_params(labelsize=max(data_size - 2, 8))
    cbar.ax.set_title("Normalized\nRatio", fontsize=max(data_size - 2, 8), pad=6, fontweight="bold")

    ax.set_xlabel("")
    ax.set_ylabel("")

    # 设置字体大小
    if len(gene_data) <= 20:
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
    elif len(gene_data) <= 50:
        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=12)
    elif len(gene_data) <= 80:
        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=12)
    else:
        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=10)

    # 设置标签

    x_labels = []
    for col in df_color.columns:
        if isinstance(col, tuple):
            col = col[0]
        x_labels.append(str(col).replace(" ", "\n"))
    ax.set_xticklabels(x_labels)

    for label in ax.get_xticklabels():
        label.set_rotation(0)
        label.set_ha('center')
        label.set_va('top')
        label.set_multialignment('center')

    # 标题
    ax.set_title('OTTER-AI Total Scores for Interacting Proteins', fontsize=20,
                 fontdict={'family': 'Arial', 'weight': 'bold'},
                 pad=10)

    plt.savefig(path_o, dpi=800)
    plt.close()


def heatmap_ppis(path_data, path_dir, output_fn = " .png",
                   len_num = 50, min_number=0):
    """
    header_data:标题
    data_l:data list
    col_data:columns=,各列的名称
    index_data:index=,gene列表
    path_dir:要保存的路径文件夹
    output_fn = " .png":输出的文件夹名
    len_num:选择用多少个基因绘制热图
    """


    # 定义输出的图片文件名
    if output_fn == " .png":
        output_fn = "Heatmap_PPI_Score.png"
    else:
        pass

    path_o = os.path.join(path_dir, output_fn)

    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.weight'] = 'bold'

    # 读取数据
    if str(path_data).endswith(".csv"):
        data_read = pd.read_csv(path_data)
    else:
        data_read = pd.read_excel(path_data)
    cols = data_read.columns.tolist()[1:]
    # 提取数据
    gene_col = data_read.iloc[:, 0]
    gene_col_list = gene_col.tolist()
    if len(gene_col_list) > len_num:
        gene_data = gene_col_list[: len_num]
    else:
        gene_data = gene_col_list
    data = data_read.iloc[: len(gene_data), 1:]

    gene_name_list = data_read.iloc[: len(gene_data), 0].tolist()
    data_array = np.array(data)
    data_list = data_array.tolist()

    df = pd.DataFrame(data_list,
                      columns=cols,
                      index=gene_name_list)

    # 制作热图
    if len(gene_data) <= 20:
        fig, ax = plt.subplots(figsize=(10, 10))
    elif len(gene_data) <= 40:
        fig, ax = plt.subplots(figsize=(10, 14))
    elif len(gene_data) <= 60:
        fig, ax = plt.subplots(figsize=(10, 16))
    elif len(gene_data) <= 80:
        fig, ax = plt.subplots(figsize=(10, 20))
    else:
        fig, ax = plt.subplots(figsize=(10, 25))
    sns.set(font_scale=1.2)


    # 自定义颜色
    colors = ["#f8f7fc", "#ffc0cb", "#e18e96"]
    df_raw = df.copy()
    df_color = df.copy().astype(float)

    # 根据最大值
    max_score = len(df.columns) - 1
    if max_score < 1:
        max_score = 1

    cmap, norm = create_three_color_colormap(colors, max_score)

    # 生成注释
    annot_df = df_raw.copy()
    for col in annot_df.columns:
        annot_df[col] = annot_df[col].map(lambda x: f"{x:.2f}")

    # 字体大小
    if len(gene_data) <= 20:
        data_size = 16
    elif len(gene_data) <= 30:
        data_size = 14
    elif len(gene_data) <= 50:
        data_size = 12
    else:
        data_size = 10

    # cbar
    fig.subplots_adjust(left=0.18, right=0.82, top=0.92, bottom=0.12)
    cbar_ax = fig.add_axes((0.88, 0.15, 0.02, 0.22))
    heatmap = sns.heatmap(df_color, annot=annot_df, cmap=cmap, fmt='',
                          norm=norm,
                          linewidths=0.4, linecolor='white',
                          annot_kws={"size": data_size},
                          ax=ax,
                          cbar_ax=cbar_ax
                          )

    cbar = heatmap.collections[0].colorbar
    if max_score <= 1:
        cbar.set_ticks([0, 1])
        cbar.set_ticklabels(["0", "1"])
    else:
        cbar.set_ticks([0, max_score])
        cbar.set_ticklabels([f"{0:g}",f"{max_score:g}"])

    cbar.ax.tick_params(labelsize=max(data_size - 2, 8))
    for label in cbar.ax.get_yticklabels():
        label.set_fontweight('bold')

    cbar.ax.set_title("Score",
                      fontsize=max(data_size - 2, 8),
                      pad=10,
                      fontweight='bold')

    ax.set_xlabel("")
    ax.set_ylabel("")

    if len(gene_data) <= 20:
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
    elif len(gene_data) <= 50:
        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=12)
    elif len(gene_data) <= 80:
        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=12)
    else:
        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=10)

    x_labels = []
    for col in df_color.columns:
        if isinstance(col, tuple):
            col = col[0]
        x_labels.append(str(col).replace(" ", "\n"))
    ax.set_xticklabels(x_labels)

    for label in ax.get_xticklabels():
        label.set_rotation(0)
        label.set_ha('center')
        label.set_va('top')
        label.set_multialignment('center')

    # 标题
    ax.set_title('OTTER-AI Keywords Scores for Interacting Proteins', fontsize=20,
              fontdict={'family': 'Arial', 'weight': 'bold'},
              pad=10)

    plt.savefig(path_o, dpi=800)
    plt.close()

