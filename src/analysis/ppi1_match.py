"""
获取每个基因的ppi相信
"""

import json
from pathlib import Path
import os
import time
import glob
import math
import pandas as pd
import numpy as np
import operator
import csv
from collections import OrderedDict
from collections import defaultdict
from itertools import combinations


# 定义函数读取文件
def read_txt_file(txt_name):
    """
    txt_name:file name
    return:将文件中的内容按行读进一个列表，返回该列表
    """
    path_txt = Path(txt_name)
    contents_txt = path_txt.read_text(encoding='UTF-8')
    txt_list = contents_txt.splitlines()
    return txt_list

def read_json_file(file_name_json):
    """
    file_name_json:file name
    return:返回字典内容
    """
    with open(file_name_json, 'r', encoding='UTF-8', errors='ignore') as f_read_json:
        dict_json = json.load(f_read_json)

    return dict_json

# 定义函数，输出JSON文件
def json_output(file_name_output, dict_for_json):
    """
    file_name_output:要输出的文件名
    dict_for_json:要输出的内容
    """
    path = Path(file_name_output)
    contents = json.dumps(dict_for_json, ensure_ascii=False, indent=4)
    path.write_text(contents, encoding='UTF-8')


# 定义函数，输出txt文件
def txt_output(txt_name_output, txt_contents):
    """
    txt_contents文件内容
    txt_name_output：待输出的文件名
    txt_contents：待输出的内容
    """

    with open(txt_name_output, 'w', encoding='UTF-8') as f_txt_output:
        if isinstance(txt_contents, list):
            for txt_content in txt_contents:
                f_txt_output.write(f"{txt_content}\n")
        if isinstance(txt_contents, set or tuple):
            for txt_content in txt_contents:
                f_txt_output.write(f"{txt_content}\n")
        if isinstance(txt_contents, dict):
            for k, (num, ratio) in txt_contents.items():
                ratio = float(ratio)
                f_txt_output.write(f"{k} \t{num} \t{ratio:.4f}\n")

# 定义函数，确保指定的文件夹存在
def ensure_dir_exists(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)

    else:
        pass


# 定义函数，输出一个空文件
def emtpy_file_output(file_name_empty):
    """
    新建一个空文件
    file_name_empty：文件名
    """
    with open(file_name_empty, 'w'):
        pass

# 定义函数分割txt文件中的内容
def split_txt1_d(path_file):
    """
    path_file:待处理文件名
    返回值：字典{gene:[ppi数量,gene总数,ratio}
    """

    file_l = read_txt_file(path_file)
    s_d = {}
    for line_info in file_l:
        line_info = line_info.split(" ")
        s_d[line_info[0]] = line_info[-1]

    return s_d


# 定义函数分割内容
def split_txt2_l(path_file):
    """
    path_file:待处理文件名
    返回值：列表
    """

    # 将文件读取为列表
    file_l = read_txt_file(path_file)
    s_l = []
    for line_info in file_l:
        line_info = line_info.split(" ")
        s_l.append(line_info)

    return s_l


# 定义函数，获取gene symbol到ensp编号的转换字典，过滤基因
def symbol_clean(symbol_input):
    """
    symbol_input:输入基因
    return：数据清洗后的字典{symbol:ENSP_num}
    """

    path_eg = "Symbol_To_EG.txt"
    path_ensp = "EG_To_ENSP.txt"

    eg_d = split_txt1_d(path_eg)
    ensp_d = split_txt1_d(path_ensp)

    # 创建输出信息
    s_ensp_d = {}
    for symbol in symbol_input:
        # 得到gene1的ensp信息
        if symbol in eg_d.keys():
            eg_symbol = eg_d[symbol]
            if eg_symbol in ensp_d.keys():
                ensp_info = ensp_d[eg_symbol]

                # 记录输出信息
                s_ensp_d[symbol] = ensp_info

            else:
                pass


        else:
            pass

    return s_ensp_d


# 定义函数，根据统计数量计算ppi得分
def ppi_score_ratio(num_data_dict):
    """
    num_data_dict : {gene:ppi_num_count}
    统计最大值，以ppi_num除以最大值，得到ppi分数
    return:{gene:[ppi_num,ratio]}
    """

    # 统计ppi计数的最大值
    num_counts = [float(num_c) for num_c in num_data_dict.values()]
    max_num = max(num_counts, default=0)
    # 计算得分
    if max_num != 0:
        score_dict = {}
        for gene, num_count in num_data_dict.items():
            gene_score = num_count / max_num
            score_dict[gene] = [num_count, gene_score]

    else:
        score_dict = {
            gene: [num_count, 0]
            for gene, num_count in num_data_dict.items()
        }

    # 排序
    score_data = dict(sorted(score_dict.items(), key=lambda item: item[1][-1], reverse=True))

    return score_data


def ppi_gene_match(gene_input_list):
    """
    gene_input_list：输入需要查询ppi数量的基因
    return：分数结果
    """

    if gene_input_list:
        # 读取输入文件
        path_ppi = "IDT/New_Data/PPI_String_Human_Link_M400_New.txt"
        ppi_l = split_txt2_l(path_ppi)

        # 创建输出数据
        genes_ppi = {}
        genes_input_d = symbol_clean(gene_input_list)
        # 字典：每个基因与哪些基因存在ppi
        ppi_gg_d = {}
        ppi_gg_n_d = {}

        ppi_gg = defaultdict(set)
        ppi_cs_d = defaultdict(dict)

        for ppi_line in ppi_l:
            try:
                gene_a, gene_b = ppi_line[0], ppi_line[1]
                c_score = float(ppi_line[-1]) if "." in ppi_line[-1] else int(ppi_line[-1])
                ppi_gg[gene_a].add(gene_b)
                ppi_gg[gene_b].add(gene_a)
                ppi_cs_d[gene_a][gene_b] = c_score
                ppi_cs_d[gene_b][gene_a] = c_score
            except (IndexError, ValueError) as e:
                print(f"Invalid PPI line: {ppi_line}, error: {str(e)}")
                continue

        # 预转换基因列表
        gene_e_list = list(genes_input_d.items())
        ensp_to_gene = {v: k for k, v in genes_input_d.items()}

        for idx, (gene_input1, ensp1) in enumerate(gene_e_list):
            # 获取所有与gene_input1存在交互的ENSP
            interacting_ensps = ppi_gg.get(ensp1, set())

            if interacting_ensps:
                # 转换为基因名并过滤
                genes_gg = []
                for ensp in interacting_ensps:
                    if ensp in ensp_to_gene and ensp_to_gene[ensp] != gene_input1:
                        genes_gg.append(ensp_to_gene[ensp])

                # 构建PPI关系结构
                ppi_entries = []
                for gene_gg in genes_gg:
                    ppi_entries.append({
                        'name': "",
                        'source': gene_input1,
                        'target': gene_gg
                    })

            else:
                genes_gg = []
                ppi_entries = []

            # 写入结果
            ppi_gg_n_d[gene_input1] = len(genes_gg)
            ppi_gg_d[gene_input1] = ppi_entries

        # 计算ppi得分
        genes_ppi_score_d = ppi_score_ratio(ppi_gg_n_d)
        genes_ppi.update(genes_ppi_score_d)

        # 得到差集
        gene_di_l = list(set(gene_input_list).difference(set(genes_input_d)))
        # 赋值为空列表
        ppi_gene_di_d = {gene_di: [0, 0] for gene_di in gene_di_l}
        ppi_gg_di_d = {gene_di: [] for gene_di in gene_di_l}
        # 将字典中的所有元素更新写入总输出字典
        ppi_gg_d.update(ppi_gg_di_d)
        genes_ppi.update(ppi_gene_di_d)
        print(f"len(ppi_gg_d, genes_ppi) final search:{len(genes_ppi)}, {len(ppi_gg_d)}")
        # 对基因排序
        sorted_dict = sorted(genes_ppi.items(), key=lambda item: item[1][-1], reverse=True)
        genes_ppi_new = dict(sorted_dict)

        return genes_ppi_new, ppi_gg_d

    else:
        return {}, {}


# 定义函数进行DEGs的ppi筛选
def gene_ppi_mining(gene_input_list):
    """
    PPI板块函数
    对输入的gene，扫描ppi列表，匹配互作蛋白（不与输入gene重复）
    return：{gene：interaction gene list}
    """
    if gene_input_list:
        path_ppi = "PPI_String_Human_Link_M400_New.txt"
        # 将文件内容每一行按空格分割
        ppi_l = split_txt2_l(path_ppi)
        # 读取输入文件
        path_eg = "Symbol_To_EG.txt"
        path_ensp = "EG_To_ENSP.txt"

        # 将文件内容每一行按空格分割
        eg_d = split_txt1_d(path_eg)
        ensp_d = split_txt1_d(path_ensp)

        # 数据清洗，得到genes_input dict
        genes_input_d = symbol_clean(gene_input_list)
        ppi_gene_d = {}
        #————————————————————————————————————————————————————————————————————

        # 预处理：构建PPI互作原始信息的字典
        # g、ig用ensp表示
        ppi_dict = defaultdict(set)
        for ppi in ppi_l:
            if len(ppi) >= 2:
                ensp1, ensp2 = ppi[0], ppi[1]
                ppi_dict[ensp1].add(ensp2)
                ppi_dict[ensp2].add(ensp1)

        # 构建反向映射：ensp-eg-gene
        ensp_to_eg = {ensp: eg for eg, ensp in ensp_d.items()}
        eg_to_gene = {eg: gene for gene, eg in eg_d.items()}

        # 输入基因列表转集合
        input_gene_names = set(gene_input_list)


        for gene_input1, ensp_info1 in genes_input_d.items():
            try:
                # 1.查询PPI互作ensp
                ig_ensps = ppi_dict.get(ensp_info1, set())

                # 2. 转换ensp到gene并过滤
                gene1_ig_genes = set()
                for ig_ensp in ig_ensps:
                    # ensp-eg-gene
                    ig_eg = ensp_to_eg.get(ig_ensp)
                    if ig_eg is None:
                        continue

                    ig_gene = eg_to_gene.get(ig_eg)
                    if ig_gene is None:
                        continue

                    if ig_gene != gene_input1 :
                        gene1_ig_genes.add(ig_gene)

                ppi_gene_d[gene_input1] = list(gene1_ig_genes)

            except Exception as e:
                ppi_gene_d[gene_input1] = []
                print(f"处理 {gene_input1} 时出错: {str(e)},结果为空：{ppi_gene_d}")


        # symbol input中的PPI检索结束，为不在里面的gene赋值0
        gene_di_l = list(set(gene_input_list).difference(set(genes_input_d)))
        gene_di_d = {gene_di: [] for gene_di in gene_di_l}
        # 写入总输出字典
        ppi_gene_d.update(gene_di_d)

        # 对基因进行排序
        sorted_dict = dict(sorted(ppi_gene_d.items(), key=lambda item: len(item[1]), reverse=True))

        return sorted_dict

    else:
        return {}

# 定义函数，在对DEG的互作蛋白进行文本挖掘之后，对DEG进行评分
def ppi_m_score(m_data_dict):
    """
    PPI板块函数
    # {k: [v[0], v[-1]]},k:i_gene;v[0]:有效报道数，v[-1]：该i_gene在pubmed中的报道占比(ratio)
    m_data_dict:{g1:{ig1:[v[0], v[-1]]},{ig2: [v[0], v[-1]]}...,g2：{...}}
    # m_data_dict:{g1:{ig1:num},{ig2:num}...,g2：{...}}
    根据所有ratio的加和分数，对DEG进行打分排序
    return： {g1:[ppi_num, ppi_report-num,score],g2:[...]...}
    """
    m_data_tn_d = {}
    for mg, ig_data in m_data_dict.items() if m_data_dict else []:
        # 计算该基因的互作蛋白数量
        mg_ppi_num = len(ig_data.keys())
        # 计算该基因的互作蛋白中有报道的数量
        ppi_r_l = [k for k, v in ig_data.items() if v[0] > 0]
        ppi_n_r_l = [k for k , v in ig_data.items() if v[0] == 0]
        mg_ppi_re_num = len(ppi_r_l)
        # 为g计算score
        if mg_ppi_re_num > 0:
            ig_sum = sum([v[-1] for k, v in ig_data.items()])
            mg_score = (ig_sum / mg_ppi_re_num) / mg_ppi_re_num

        else:
            mg_score = 0

        try:
            mg_tn = sum(v[0] for v in ig_data.values())
        except ValueError as e:
            print(f"{e}:warning. mg_tn will be zero.")
            mg_tn = 0

        m_data_tn_d[mg] = {"PPI_Total": mg_ppi_num, "PPI_Report":mg_ppi_re_num, "PM_Num": mg_tn, "Score": mg_score}

    # 排序
    m_g_score_d = dict(sorted(m_data_tn_d.items(), key=lambda item: item[1]["Score"], reverse=True))

    # 返回结果
    return m_g_score_d


# 定义函数，对csv文件添加条件筛选工作表
def data_rank_sheet(data_path):
    """
    data_path: 输入数据路径
    """
    # 读取CSV文件
    df = pd.read_csv(data_path)
    df_name = str(data_path).split(".csv")[0]

    # 获取PPI列的最大值
    max_ppi = df['PPI_Total'].max()

    # 创建rank数据框，筛选满足条件的行
    if max_ppi >= 100:
        rank_data = (df['PPI_Total'] >= 100) & (df['PPI_Report'] >= 50) & (df['PM_Num'] >= 300)
    elif max_ppi <= 50:
        rank_data = (df['PPI_Total'] >= 25) & (df['PPI_Report'] >= 1) & (df['PM_Num'] >= 10)
    else:
        rank_data = (df['PPI_Total'] >= 50) & (df['PPI_Report'] >= 25) & (df['PM_Num'] >= 150)

    # 应用筛选条件
    rank_df = df[rank_data]
    mask = rank_data.astype(bool)
    rank_df = df.loc[mask, :].copy()

    score_data = rank_df["Score"]

    max_score = score_data.max()
    max_score = max_score if max_score != 0 else 1.0

    rank_df.loc[:, "Score"] = np.sqrt(np.sqrt(score_data / max_score))

    # 结果写入
    with pd.ExcelWriter(f'{df_name}.xlsx', engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name="data", index=False)
        rank_df.to_excel(writer, sheet_name='rank', index=False)

    return rank_df



# 定义函数，读取csv文件，对数据进行归一化打分
def csv_data_normalize(csv_data_path, normalize_path, col_name="Score"):
    df = pd.read_csv(csv_data_path)
    df[col_name] = df[col_name].fillna(0)
    score_data = df[col_name]
    # 获取Score列最大值
    max_score = score_data.max()
    if max_score == 0:
        df[col_name] = 0.0
    else:
        # 执行归一化公式,替换原Score列
        df[col_name] = np.sqrt(np.sqrt(score_data / max_score))
    col_last = df.columns[-1]
    col_second = df.columns[1]
    # 按归一化后的值排序
    df.sort_values(by=[col_name, col_last, col_second], ascending=[False, False, False], inplace=True)

    # 保存
    df.to_csv(normalize_path, index=False)

# 定义函数进行归一化打分处理
def normalize_scores(gene_dict):
    """
    对基因得分进行归一化处理（最高分设为1）
    参数：
    gene_dict (dict): 原始基因数据字典，格式为 {基因名: [..., 最终得分]}

    返回值：
    dict: 归一化后的字典 {基因名: 归一化得分}
    """

    if gene_dict:
        try:
            scores = [values[-1] for values in gene_dict.values() if len(values) >= 1]

            if not scores:
                raise ValueError("输入字典中没有有效得分数据")
            if len(gene_dict) == 1:
                for gene, score_data in gene_dict.items():
                    score = score_data[-1]
                    return {gene: score}, {gene: score_data}

            else:
                max_score = max(scores)

                if max_score == 0:
                    return {gene0: 0.0 for gene0 in gene_dict}, gene_dict

                # 生成归一化字典
                normalize_dict = {
                    gene: math.sqrt(math.sqrt(values[-1] / max_score))
                    for gene, values in gene_dict.items()
                    if len(values) >= 1
                }

                # 对基因进行排序
                normalize_dict1 = sorted(normalize_dict.items(), key=lambda item: item[1], reverse=True)
                normalize_dict_s = dict(normalize_dict1)
                gene_out_d = {
                    gene: gene_dict[gene][: -1] + [n_score]
                    for gene, n_score in normalize_dict_s.items()
                }
                return normalize_dict_s, gene_out_d

        except (TypeError, IndexError) as e:
            print(f"数据格式错误：{str(e)}")
            return {}, gene_dict
        except Exception as e:
            print(f"未知错误：{str(e)}")
            return {}, gene_dict
    else:

        return {}, gene_dict


# 定义函数，对不同关键词对应的PPI模块取交集
def ppi_info_i(input_path, out_path="PPI_Data.xlsx"):
    """
    input_path:文件路径
    """

    # 读取所有CSV文件并存储为字典列表
    dict_list = []
    common_headers = None

    csv_files = []
    for item in os.listdir(input_path):
        item_path = os.path.join(input_path, item)
        if os.path.isdir(item_path) and item.endswith("_output_info"):
            ppi_info_path = os.path.join(item_path, "PPI_DE_Score.csv")
            csv_files.append(ppi_info_path)

    if not csv_files:
        return []

    for file_path in csv_files:
        with open(file_path, 'r', encoding='UTF-8') as f:
            reader = csv.DictReader(f)
            headers = reader.fieldnames

            # 校验标题一致性
            if common_headers is None:
                common_headers = headers
            elif headers != common_headers:

                continue

            # 构建当前文件的字典
            data_dict = {}
            for row in reader:
                gene_name = row[headers[0]]
                gene_row_d = {}
                for header in headers:
                    if header == headers[0]:
                        continue
                    try:
                        gene_row_d[header] = float(row[header])
                    except (ValueError, TypeError):
                        gene_row_d[header] = row[header]
                data_dict[gene_name] = gene_row_d
            dict_list.append(data_dict)

    # 检查是否有效
    if not dict_list:
        return []

    # 获取所有文件的共同gene—name
    gene_inter = set(dict_list[0].keys())
    for d in dict_list[1:]:
        gene_inter.intersection_update(d.keys())

    if not gene_inter:
        return []

    # 计算各列总和
    column_sums = defaultdict(lambda: defaultdict(float))

    for name in gene_inter:
        for d in dict_list:
            for col in common_headers[1:]:
                if col in d[name]:
                    try:
                        if isinstance(d[name][col], (int, float)):
                            column_sums[name][col] +=  d[name][col]
                    except KeyError:
                        pass

    # 排序
    sorted_names = sorted(
        gene_inter,
        key=lambda x: column_sums[x].get('Score', 0),
        reverse=True
    )

    # 结果数据
    result = []
    for name in sorted_names:
        row_data = {common_headers[0]: name}
        for col in common_headers[1:]:
            value = column_sums[name].get(col, 0)
            if isinstance(value, (int, float)):
                row_data[col] = round(value, 4)
            else:
                row_data[col] = value
        result.append(row_data)

    # 创建DataFrame
    df = pd.DataFrame(result, columns=common_headers)
    if not df.empty and 'PPI_Total' in df.columns:
        # 找到df中Score列的最大绝对值
        max_ppi_total_num = df['PPI_Total'].abs().max()
    else:
        max_ppi_total_num = 30
    # 创建rank数据框，筛选满足条件的行
    if max_ppi_total_num >=200:
        rank_condition = (df['PPI_Total'] >= 200) & (df['PPI_Report'] >= 100) & (df['PM_Num'] >= 600)
    elif max_ppi_total_num <= 45:
        rank_condition = (df['PPI_Total'] >= 10) & (df['PPI_Report'] >= 1) & (df['PM_Num'] >= 5)
    else:
        rank_condition = (df['PPI_Total'] >= 45) & (df['PPI_Report'] >= 5) & (df['PM_Num'] >= 15)

    rank_df = df[rank_condition].copy()
    if not rank_df.empty and 'Score' in rank_df.columns:
        # 找到rank_df中Score列的最大绝对值
        max_score_abs_rank = rank_df['Score'].abs().max()
        scale_factor_rank = 1
        if max_score_abs_rank > 0:
            exponent = 0
            while max_score_abs_rank * (10 ** exponent) < 1:
                exponent += 1
            scale_factor_rank = 10 ** exponent
            rank_df['Score'] = rank_df['Score'] * scale_factor_rank
    else:
        scale_factor_rank = 1


    # 写入Excel文件
    with pd.ExcelWriter(
            os.path.join(input_path, out_path),
            engine='xlsxwriter'
    ) as excel_out:
        df.to_excel(excel_out, index=False, sheet_name='all_data')
        if not rank_df.empty:
            rank_df.to_excel(excel_out, index=False, sheet_name='rank')
        else:
            pass


        workbook = excel_out.book
        worksheet = excel_out.sheets['all_data']
        # 创建数字格式（保留4位小数）
        number_format = workbook.add_format({'num_format': '0.0000'})

        # 应用数字格式到所有数值列
        for i, col in enumerate(df.columns):
            if col != common_headers[0]:

                # 获取当前列数据
                col_data = df[col]
                for j, value in enumerate(col_data):
                    # 如果值是浮动数（小数），应用四位小数格式
                    if isinstance(value, float) and not value.is_integer():
                        col_idx = df.columns.get_loc(col)
                        worksheet.write(j + 1, col_idx, value, number_format)
                    else:
                        col_idx = df.columns.get_loc(col)
                        worksheet.write(j + 1, col_idx, value)

        # 只有在rank工作表存在时才应用格式
        if not rank_df.empty:
            try:
                rank_worksheet = excel_out.sheets['rank']
                # 应用数字格式到所有数值列
                for i, col in enumerate(rank_df.columns):
                    if col != common_headers[0]:
                        col_data = rank_df[col]
                        for j, value in enumerate(col_data):
                            if isinstance(value, float) and not value.is_integer():
                                col_idx = rank_df.columns.get_loc(col)
                                rank_worksheet.write(j + 1, col_idx, value, number_format)
                            else:
                                col_idx = rank_df.columns.get_loc(col)
                                rank_worksheet.write(j + 1, col_idx, value)
            except KeyError:
                pass

    return out_path



# 定义函数，对ppi模块的rank基因取交集
def header_title(name: str) -> str:
    """
    将列名转换为首字母大写。
    1.英文/数字用 '_' 作为分隔，分段首字母大写，其余小写；
    2.中文或非 ASCII 字符保持原样；
    3.去除多余空格。
    """
    s = str(name).strip().replace(" ", "_")
    parts = [p for p in s.split("_") if p != ""]
    new_parts = []
    for p in parts:
        if p.isascii():
            new_parts.append(p[:1].upper() + p[1:].lower())
        else:
            new_parts.append(p)
    return "_".join(new_parts)


def ppi_rank_data_i(input_dir, sheet_name='rank', suffix_str = "_output_info",
                    data_file = "PPI_DE_Score_N.xlsx", out_file="PPI_Rank_Data.xlsx", min_genes=10):
    """
    input_dir:输入路径
    sheet_name:需要读取的工作表名
    suffix_str :文件夹的后缀名
    data_file:单个关键词的数据文件
    output_file:输出文件
    min_genes: 单个关键词中有效基因最少数量阈值
    """

    output_path = os.path.join(input_dir, out_file)
    input_dir = Path(input_dir)
    # ————1.读取单关键词的子目录，生成kw_dict，并记录显示列名
    kw_data_list = []
    del_kw = []
    for subdir in input_dir.iterdir():
        if not subdir.is_dir():
            continue
        data_dir_name = subdir.name
        if data_dir_name.endswith(suffix_str):
            kw_name = data_dir_name.split(suffix_str)[0].strip()
            excel_path = subdir / data_file
            if not excel_path.exists():
                continue

            df = pd.read_excel(excel_path, sheet_name=sheet_name)
            df.columns = df.columns.astype(str).str.strip()
            df_lower = df.copy()
            df_lower.columns = df_lower.columns.str.lower()
            # 确认数据列名
            data_col = {col_name.lower(): col_name for col_name in["Symbol_Name", "Score"]}
            if not set(data_col.keys()).issubset(set(df_lower.columns)):
                raise KeyError(f"{excel_path}缺少列{data_col}，实际列：{list(df.columns)}")

            # 读取数据
            s = pd.to_numeric(df_lower["Score".lower()], errors="coerce").fillna(0.0).round(2)
            s.index = df_lower["Symbol_Name".lower()].astype(str)
            s.name = kw_name
            rank_genes_num = int(np.count_nonzero(s))
            if rank_genes_num == 0 or rank_genes_num < min_genes:
                del_kw.append((kw_name, rank_genes_num))
                continue

            kw_data_list.append(s)

    if not kw_data_list:
        raise RuntimeError("未收集到任何有效数据。")

    # ————2.仅保留所有文件夹的交集基因，并逐列相加得到总分
    if len(kw_data_list) == 1:
        df_scores = kw_data_list[0].to_frame()

    else:
        # 先尝试所有有效关键词的共同交集
        df_scores_all = pd.concat(kw_data_list, axis=1, join="inner")

        if not df_scores_all.empty:
            df_scores = df_scores_all

        else:
            # 共同交集为0，改为两两取交集，选则交集最多的一对

            best_df = None
            best_pair = None
            best_size = 0
            for s1, s2 in combinations(kw_data_list, 2):
                inter_genes = pd.concat([s1, s2], axis=1, join="inner")

                gene_size = inter_genes.shape[0]
                if gene_size > best_size:
                    best_size = gene_size
                    best_df = inter_genes
                    best_pair = (s1.name, s2.name)

            if best_df is None or best_size == 0:
                raise RuntimeError("无法生成有效数据。")
            df_scores = best_df


    for col in df_scores.columns:
        df_scores[col] = pd.to_numeric(df_scores[col], errors="coerce").fillna(0.0).round(2)
    df_scores["Total_Score"] = df_scores.sum(axis=1)
    # ————3.输出Excel结果
    df_out = (
        df_scores.reset_index()
        .rename(columns={"index": "Symbol_Name".lower()})
    )

    df_out.columns = [header_title(c) for c in df_out.columns]
    cols = ["Symbol_Name"] + [c for c in df_out.columns if c not in ("Symbol_Name", "Total_Score")] + ["Total_Score"]
    df_out = df_out.loc[:, cols].sort_values("Total_Score", ascending=False)  # 降序排序

    with pd.ExcelWriter(output_path, engine="openpyxl") as w:
        df_out.to_excel(w, sheet_name="Score", index=False)

    return output_path


# 定义函数，筛选结果基因对应的ppi具体信息
def ppi_data_for_nodes(gg_dict, path_csv="Heatmap_Data_Normalized.csv", len_num=50,
                       output_json="Final_ECharts_PPI_Data.json",
                       output_score="Final_ECharts_PPI_Freq.txt"):
    """
    处理交集数据
    gg_dict:ppi具体信息的字典：{gene1:[{'name':'', 'source':gene1, 'target': gene2}]}
    path_csv:含有各基因最终打分与排名信息的文件
    len_num:最多提取多少个基因
    output_json:输出文件的文件名
    """

    try:

        with open(path_csv, 'r', newline='') as csv_input:
            csv_reader = csv.reader(csv_input)
            next(csv_reader, None)

            with open(output_score, 'w') as txt_score:
                for row in csv_reader:
                    new_row = "   ".join(row[1:])
                    txt_score.write(new_row + "\n")


        # 读取csv数据文件
        df = pd.read_csv(path_csv)
        # 提取gene名列
        gene_name_c = df.iloc[:, 1]
        gene_name_all_l = gene_name_c.tolist()
        if len(gene_name_all_l) > len_num:
            gene_names_l = gene_name_all_l[:len_num]
        else:
            gene_names_l = gene_name_all_l

        gene_names = df.iloc[:, :2]
        # 提取打分数据

        csv_data = {}
        with open(path_csv, 'r', encoding='UTF-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                gene = row['Gene name']
                if gene in gene_names_l:
                    try:
                        rank = int(row['Rank'])
                        total = float(row['Total score'])

                        csv_data[gene] = [rank, gene, total]
                    except (ValueError, KeyError) as e:
                        print(f"CSV数据格式错误-{gene}：{str(e)}")
    except FileNotFoundError:
        return

    # 筛选有效基因
    valid_genes = [g for g in gene_names_l
                   if g in csv_data and g in gg_dict]

    # 按rank排序基因
    sorted_genes = sorted(valid_genes,
                          key=lambda x: csv_data[x][0])

    result = OrderedDict()
    for gene in sorted_genes:
        result[gene] = [
            csv_data[gene],
            gg_dict[gene]
        ]

    # 写入JSON文件
    with open(output_json, 'w', encoding='UTF-8') as f_json:
        json.dump(result, f_json,
                  indent=2,
                  ensure_ascii=False,
                  default=lambda x: x.__dict__)

    return result


# 定义函数，筛选结果基因对应的ppi具体信息
def ppi_data_top_nodes(gg_dict, path_csv="Heatmap_Data_Normalized.csv", len_num=50,
                       output_json="Final_ECharts_PPI_Data.json",
                       output_score="Final_ECharts_PPI_Freq.txt"):
    """
    处理单个关键词的数据，top100数据
    gg_dict:ppi具体信息的字典：{gene1:[{'name':'', 'source':gene1, 'target': gene2}]}
    path_csv:含有各基因最终打分与排名信息的文件
    len_num:最多提取多少个基因
    output_json:输出文件的文件名
    """

    try:

        with open(path_csv, 'r', newline='') as csv_input:
            csv_content = csv.reader(csv_input)

            next(csv_content, None)

            # 收集所有基因的总分用于排名
            gene_scores = []
            all_rows = []

            for row in csv_content:
                all_rows.append(row)
                gene_name = row[0]
                try:
                    total_score = float(row[-1])
                    gene_scores.append((gene_name, total_score))
                except (ValueError, IndexError) as e:
                    print(f"分数解析错误：{gene_name}   {str(e)}")

            # 排序
            gene_scores.sort(key=lambda x: x[1], reverse=True)

            ranking = {}
            current_rank = 1
            prev_score = None

            for i, (gene, score) in enumerate(gene_scores):
                if score == prev_score:
                    ranking[gene] = current_rank
                else:
                    current_rank = i + 1
                    ranking[gene] = current_rank
                    prev_score = score

        # 写入文件
        with open(output_score, 'w') as txt_score:
            for row in all_rows:
                txt_score.write("   ".join(row) + "\n")

        # 收集JSON数据
        csv_data = {}
        for row in all_rows:
            gene = row[0]
            try:
                rank = ranking.get(gene, 0)
                total = float(row[-1])
                csv_data[gene] = [rank, gene, total]
            except (ValueError, IndexError) as e:
                print(f"数据格式错误：{gene}   {str(e)}")

        # 筛选有效基因并排序
        valid_genes = sorted(
            [g for g in list(csv_data.keys())[:len_num] if g in gg_dict],
            key=lambda x: csv_data[x][0]
        )

        # 构建结果字典并写入JSON
        result = OrderedDict()
        for gene in valid_genes:
            result[gene] = [csv_data[gene], gg_dict[gene]]

        with open(output_json, 'w', encoding='UTF-8') as f_json:
            json.dump(result, f_json, indent=2, ensure_ascii=False)

    except FileNotFoundError:
        print(f"CSV文件不存在：{path_csv}")
        result = {}

    return result

# 定义函数，筛选结果基因对应的ppi具体信息
def ppi_data_dp_nodes(ppi_dir, len_data=""):
    """
    处理PPI板块的数据，为绘制节点图准备数据
    ppi_dir: PPI板块的文件夹路径
    len_data=""：数据长度，默认空即处理所有基因，可以输入数字(整数)
    """

    json_path = os.path.join(ppi_dir, 'PPI_DE_PM_Num.json')
    ppi_i_path = os.path.join(ppi_dir, 'PPI_Data.xlsx')

    # 读取JSON数据
    ppi_all = read_json_file(json_path)

    # 读取Excel文件中的交集基因数据
    df_excel = pd.read_excel(ppi_i_path)
    ppi_i_data = df_excel.sort_values(by='Score', ascending=False, inplace=False)

    # 在总文件夹下生成文件
    dp_output_files(ppi_i_data, ppi_all, ppi_dir, data_len=len_data, kw=None)

    # 处理每个关键词对应的子文件夹
    kws_dir = glob.glob(os.path.join(ppi_dir, '*_output_info'))
    for kw_dir in kws_dir:
        keyword = os.path.basename(kw_dir).split('_output_info')[0]

        # 查找关键词文件夹中的CSV文件
        kw_genes_files = glob.glob(os.path.join(kw_dir, 'PPI_DE_Score_N.csv'))

        if kw_genes_files:
            for kw_genes_file in kw_genes_files:
                # 处理CSV文件
                df_csv = pd.read_csv(kw_genes_file)
                if 'Score' not in df_csv.columns:
                    continue

                ppi_kw_data = df_csv.sort_values(by='Score', ascending=False, inplace=False)
                # 生成对应的输出文件
                dp_output_files(ppi_kw_data, ppi_all, kw_dir, data_len=len_data, kw=f"{keyword}")
        else:
            continue

def dp_output_files(df, ppi_data, output_dir,data_len="", kw=None):
    """
    dp:data ppi
    PPI板块：为节点图生成数据
    df：ppi板块的结果数据，一般是归一化之后的排序
    ppi_data:包含所有kws的元基因与互作基因信息
    output_dir:输出数据的文件夹路径
    data_len="":空表示读取所有基因，调用函数时可以指定数据长度，如50
    prefix="": 默认无前缀，前缀一般是处理单个关键词时需要
    all_kws=True: true表示针对交集基因，统计每个基因得分时应该总和所有关键词
    kw=None：None表示针对交集基因
    """

    if data_len == "":
        gene_de = df.iloc[:, 0]
    else:
        if len(df.iloc[:, 0]) >= int(data_len):
            gene_de = df.iloc[int(data_len):, 0]
        else:
            gene_de = df.iloc[:, 0]
    if kw:
        kw_out_name = "_" + kw
    else:
        kw_out_name = ""

    with open(os.path.join(output_dir, f"Nodes_DE_Pair{kw_out_name}_PPI.txt"), 'w') as f1, \
            open(os.path.join(output_dir, f"Freq_ECharts{kw_out_name}_PPI_Pair.txt"), 'w') as f2:

        # 遍历排序后的基因
        for gene in gene_de:
            gene = str(gene)
            found = False
            all_ppi_data = {}

            if not kw:
                # 合并所有关键词下的互作基因数据
                for keyword, gene_dict in ppi_data.items():
                    if gene in gene_dict:
                        for ppi_gene, scores in gene_dict[gene].items():
                            count = scores[0] if isinstance(scores, list) and len(scores) > 1 else 0
                            score_val = scores[1] if isinstance(scores, list) and len(scores) > 1 else 0

                            if ppi_gene not in all_ppi_data:
                                all_ppi_data[ppi_gene] = [count, score_val]
                            else:
                                all_ppi_data[ppi_gene][0] += count
                                all_ppi_data[ppi_gene][1] += score_val

                # 写入合并后的数据
                for ppi_gene, scores in all_ppi_data.items():
                    count, score_val = scores
                    f1.write(f"{gene}   {ppi_gene}\n")
                    f2.write(f"{gene}   {ppi_gene}   {count}   {score_val}\n")

            elif kw and kw in ppi_data:
                # 只处理当前关键词
                gene_dict = ppi_data[kw]
                if gene in gene_dict:
                    for ppi_gene, scores in gene_dict[gene].items():
                        count = scores[0] if isinstance(scores, list) and len(scores) > 1 else 0
                        score_val = scores[1] if isinstance(scores, list) and len(scores) > 1 else 0

                        f1.write(f"{gene}   {ppi_gene}\n")
                        f2.write(f"{gene}   {ppi_gene}   {count}   {score_val}\n")


        with open(os.path.join(output_dir, f"Freq_ECharts{kw_out_name}_PPI.txt"), 'w') as f3:
            for _, row in df.iterrows():

                gene = str(row["Symbol_Name"])

                if 'Score' in row:
                    score = row['Score']

                else:
                    score = 0

                f3.write(f"{gene}   {score}\n")


def ppi_nodes_pair(input_dict_data, output_file):
    """
    将ppi对的节点数据输出为txt文件
    input_dict_data:ppi_data_for_nodes函数的输出数据
    """

    interactions = []

    for gene, gene_data in input_dict_data.items():
        # 提取基因信息
        rank = gene_data[0][0]
        score = gene_data[0][2]
        interaction_list = gene_data[1]

        # 收集所有互作对
        for interaction in interaction_list:
            source = interaction["source"]
            target = interaction["target"]
            interactions.append((rank, source, target))

    # 排序
    interactions.sort(key=operator.itemgetter(0))

    # 写入文件
    with open(output_file, 'w') as f:
        for rank, source, target in interactions:
            f.write(f"{source}   {target}\n")

