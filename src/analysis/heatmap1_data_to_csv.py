"""
结果取交集
提取交集结果对应的第二列数据，生成新的排名，然后输出csv文件
对数据进行归一化处理

"""
import pandas as pd
from pathlib import Path
import csv
from collections import defaultdict


# 定义函数读取csv文件,筛选otter结果有效数据
def get_csv_file(keyword, otter_path_name):
    df_file = pd.read_csv(fr"{otter_path_name}\Output_{keyword}\Test_Score_Table_Top_50.csv")
    total_d = dict(zip(df_file['Symbol'], df_file['Total_Score']))
    # 读取结果的genes转化为列表
    gene_list = df_file['Symbol'].values.tolist()

    # 筛选有效数据
    for gene in gene_list:
        if total_d[gene] == 0:
            gene_list.pop(gene_list.index(gene))

    return gene_list


# 定义函数读取excel文件
def get_excel_data(keyword, otter_path_name, name_list):
    excel_name = fr'{otter_path_name}\Freq_List_Cela_{keyword}.xlsx'
    df_e = pd.read_excel(excel_name, usecols=[0, 1], header=None, names=['Name', 'Value'])
    dict_all = dict(zip(df_e['Name'], df_e['Value']))

    # 生成新字典
    dict_match = {gene: dict_all[gene] for gene in name_list if gene in dict_all}
    return dict_match


# 定义函数读取txt文件
def read_txt_file(txt_name):
    path_txt = Path(txt_name)
    contents_txt = path_txt.read_text(encoding='UTF-8')
    txt_list = contents_txt.splitlines()
    return txt_list


# 定义函数取交集
def get_intersection(gene1, gene2, gene3, gene4, gene5):
    path_i = Path(f"Inter_Gene_5Groups.txt")
    inters = set.intersection(set(gene1), set(gene2), set(gene3), set(gene4), set(gene5))
    write_i = ''

    for gene_name in inters:
        write_i += f"{gene_name}\n"
    path_i.write_text(write_i, encoding="UTF-8")
    return write_i


# 定义函数写入csv
def text_data_to_csv(message, i_l, kw_l, k_g_mess, path_name='All_Kws_Gene_Ratio_Data.csv'):
    """
    message:keyword:{gene:[report_num,total_num_ratio]}
    i_l:交集
    kw_l:keywords_list
    k_g_mess:keyword:gene_name
    """

    i_l = list(i_l)
    # 1. 获取需要的所有关键词与他们对应的基因名及其附加信息
    data = {}
    gene_totals = defaultdict(float)
    for keyword in kw_l:
        if keyword in k_g_mess:
            gene_names = k_g_mess[keyword]
            if keyword in message:
                gene_info = message[keyword]
                keyword_data = {}
                for gene in i_l:
                    if gene in gene_info:
                        keyword_data[gene] = gene_info[gene]
                        gene_frequency = float(gene_info[gene][-1])
                        gene_totals[gene] += gene_frequency
                    else:
                        keyword_data[gene] = [0, 0, 0]
                data[keyword] = keyword_data
    # 2. 计算基因频率总和并排序
    sorted_genes = sorted(gene_totals.items(), key=lambda x: x[1], reverse=True)
    sorted_genes_d = dict(sorted_genes)

    # 3.输出csv文件
    headers = ['Rank', 'Gene Name'] + kw_l + ['Total Score']
    with open(path_name, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([header.capitalize() for header in headers])
        for idx, (gene, total_sum) in enumerate(sorted_genes, start=1):
            row = []
            row.append(idx)
            row.append(gene)
            ratio_row = [data[keyword].get(gene, [0, 0])[-1] for keyword in kw_l]
            row.extend(ratio_row)
            row.append(total_sum)
            writer.writerow(row)

    return headers, list(sorted_genes_d.keys())


# 定义函数计算总分并写入csv文件
def total_ranked_csv(data_dict, output_filename="All_Kws_Gene_Score_Data.csv"):
    """
    处理嵌套字典并生成排序后的CSV文件
    data_dict : 嵌套字典 {keyword：{gene：{'Text_Score':值, 'PPI_Score':值, 'Total_Score':值}}}
    output_filename : 输出文件名
    """
    keywords = list(data_dict.keys())
    if not keywords:
        raise ValueError("输入字典不能为空")

    base_genes = set(data_dict[keywords[0]].keys())
    for key in keywords[1:]:
        if set(data_dict[key].keys()) != base_genes:
            raise ValueError(f"关键词 '{key}' 的基因列表不一致")

    data_l = []
    for gene in base_genes:
        gene_scores = {'gene': gene}
        total = 0.0

        for keyword in keywords:
            try:
                score = data_dict[keyword][gene]['Total_Score']
                total += float(score)
                gene_scores[keyword] = score
            except KeyError:
                raise KeyError(f"基因{gene}在关键词{keyword}中缺少Total_Score")

        gene_scores['total'] = round(total, 4)
        data_l.append(gene_scores)

    # 排序
    data_l.sort(key=lambda x: (-x['total'], x['gene']))

    with open(output_filename, 'w', newline='', encoding='UTF-8') as f:
        writer = csv.writer(f)
        header = ['Rank', 'Gene Name'] + keywords + ['Total Score']
        writer.writerow([header_name.capitalize() for header_name in header])

        for rank, item in enumerate(data_l, 1):
            row = [
                rank,
                item['gene'],
                *[format(item[key], '.4f') for key in keywords],
                format(item['total'], '.4f')
            ]
            writer.writerow(row)

    return header



# 定义归一化算法函数
def normalize_column(column):
    max_value = float(column.max())
    if max_value > 0:
        return column / max_value
    return column


# 定义归一化操作函数
def data_normalization(path_name_i='All_Kws_Gene_Score_Data.csv',
                       path_name_o='Heatmap_Data_Normalized.csv'):
    """
    path_name_i:input归一化的数据源文件名
    path_name_o:output归一化的数据源文件名
    """
    df = pd.read_csv(path_name_i)
    gene_names = df.iloc[:, :2]
    data_to_normalize = df.iloc[:, 2:-1]

    # 对数据应用归一化处理
    data_normalized = data_to_normalize.apply(normalize_column)
    data_normalized = data_normalized.round(4)
    # 计算新的总和
    new_total_sum = data_normalized.sum(axis=1)

    # 合并
    df_normalized = pd.concat([gene_names, data_normalized, new_total_sum], axis=1)
    max_value = df_normalized[df_normalized.columns[2:]].max().max()
    min_value = df_normalized[df_normalized.columns[2:]].min().min()

    # 更新列名
    df_normalized.columns = df.columns
    df_normalized = df_normalized.sort_values(by=df_normalized.columns[-1], ascending=False).reset_index(drop=True)

    # 重新生成排名列
    df_normalized[df.columns[0]] = df_normalized.index + 1
    df_normalized.to_csv(path_name_o, index=False)

    return max_value, min_value




