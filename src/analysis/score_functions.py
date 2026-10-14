"""
定义一系列打分函数
"""

import csv
from collections import defaultdict
import math

# 定义函数从csv文件中获取各关键词的基因文本分数数据
def get_csv_data_score(csv_path, kws):
    """
    csv文件格式：排名列，基因名列，各关键词得分列，总分列
    return：
    {
        'keyword1': {'gene1': value, 'gene2': value, ...},
        'keyword2': {'gene1': value, 'gene2': value, ...},
        ...
    }
    """

    output_dict = defaultdict(dict)

    try:
        with open(csv_path, 'r', encoding='UTF-8') as csvfile:
            reader = csv.DictReader(csvfile)

            # 验证必须存在的列
            if 'Gene name' not in reader.fieldnames:
                raise ValueError("CSV文件中缺少必需的'Gene name'列")

            valid_keywords = [k for k in kws if k.capitalize() in reader.fieldnames]
            invalid_keywords = set(kws) - set(valid_keywords)

            if invalid_keywords:
                print(f"警告：以下关键词不存在于CSV列中：{', '.join(invalid_keywords)}")

            # 初始化字典结构
            for keyword in valid_keywords:
                output_dict[keyword.capitalize()] = {}


            for row_num, row in enumerate(reader, 1):
                gene = row['Gene name'].strip()
                if not gene:
                    continue

                for keyword in valid_keywords:
                    try:
                        keyword = keyword.capitalize()
                        value = row[keyword].strip()
                        if value.replace('.', '', 1).isdigit():
                            value = float(value) if '.' in value else int(value)
                        output_dict[keyword][gene] = value
                    except KeyError:
                        print(f"第{row_num}行缺少'{keyword}'列数据")
                    except Exception as e:
                        print(f"处理第{row_num}行时发生错误：{str(e)}")

        return dict(output_dict)

    except FileNotFoundError:
        return {}
    except Exception as e:
        print(f"处理文件时发生意外错误：{str(e)}")
        return {}


# 定义函数，根据统计数量计算得分
def score_ratio(num_data_dict):
    """
    num_data_dict : {gene:report_num_count}
    统计最大值，以re_num除以最大值，得到分数
    return:{gene:[re_num,ratio]}
    """

    # 统计报道计数的最大值
    num_counts = [float(num_c) for num_c in num_data_dict.values()]
    max_num = max(num_counts, default=0)
    # 计算各基因得分
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


# 定义函数进行归一化打分处理
def normalize_scores1(score_dict):
    """
    对基因得分进行归一化处理
    score_dict (dict): 原始基因数据字典，格式为 {基因名: [..., 最终得分]}

    返回值：
    dict: 归一化后的字典{基因名: 归一化得分}
    """
    try:
        # 提取所有基因的最终得分
        scores = [values[-1] for values in score_dict.values() if len(values) >= 1]
        if not scores:
            raise ValueError("输入字典中没有有效得分数据")

        max_score = max(scores, default=0)

        if max_score == 0:
            return {gene0: 0.0 for gene0 in score_dict}

        normalize_dict = {
            gene: math.sqrt(math.sqrt(values[-1] / max_score))
            for gene, values in score_dict.items()
            if len(values) >= 1
        }


        # 对基因进行排序
        normalize_dict1 = sorted(normalize_dict.items(), key=lambda item: item[1], reverse=True)
        normalize_dict_s = dict(normalize_dict1)
        score_out_d = {
            gene: score_dict[gene][: -1] + [n_score]
            for gene, n_score in normalize_dict_s.items()
        }

        return normalize_dict_s, score_out_d


    except (TypeError, IndexError) as e:
        print(f"数据格式错误：{str(e)}")
        return {}
    except Exception as e:
        print(f"未知错误：{str(e)}")
        return {}

# 定义函数进行归一化打分处理
def normalize_scores2(score_dict):
    """
    对基因得分进行归一化处理
    score_dict (dict): 原始基因数据字典，格式为 {基因名: 最终得分}

    返回值：
    dict: 归一化后的字典 {基因名: 归一化得分}
    """
    try:
        scores = [values for values in score_dict.values()]
        if not scores:
            raise ValueError("输入字典中没有有效得分数据")

        max_score = max(scores, default=0)

        if max_score == 0:
            return {gene0: 0.0 for gene0 in score_dict}

        normalize_dict = {
            gene: math.sqrt(math.sqrt(values / max_score))
            for gene, values in score_dict.items()
        }


        # 对基因进行排序
        normalize_dict1 = sorted(normalize_dict.items(), key=lambda item: item[1], reverse=True)
        normalize_dict_s = dict(normalize_dict1)
        return normalize_dict_s

    except (TypeError, IndexError) as e:
        print(f"数据格式错误：{str(e)}")
        return {}
    except Exception as e:
        print(f"未知错误：{str(e)}")
        return {}



# 定义函数，计算文本分数和ppi分数的总分
def scores_total(text_dict, ppi_dict):
    """
    text_dict: 文本分数的数据，字典{基因：得分}
    PPI_dict:PPI分数的数据，字典
    return:字典{基因：{'Text_Score':值, 'PPI_Score':值, 'Total_Score':值}}
    """
    total_d = {}
    for gene, text_num in text_dict.items():
        gene_d = {}
        gene_d["Text_Score"] = text_num
        ppi_num = ppi_dict[gene]
        gene_d["PPI_Score"] = ppi_num
        total_num = text_num + ppi_num
        gene_d["Total_Score"] = total_num

        total_d[gene] = gene_d

    return total_d


def scores_to_csv(data, filename='DE_Scores.csv'):
    """
    将基因评分数据写入CSV文件
    """
    # 统一数据格式为行字典列表
    processed_data = []
    headers = ['Symbol', 'Text_Score', 'PPI_Score', 'Total_Score']

    try:
        if isinstance(data, dict):
            if all(isinstance(v, (list, tuple)) for v in data.values()):
                processed_data = [
                    {col: values[i] for col, values in data.items()}
                    for i in range(max(len(v) for v in data.values()))
                ]
            else:
                processed_data = [
                    {'Symbol': k, **v} if isinstance(v, dict) else {'Symbol': k, 'scores': v}
                    for k, v in data.items()
                ]
        elif isinstance(data, list):
            if all(isinstance(item, dict) for item in data):
                processed_data = data
            else:
                processed_data = [
                    {headers[i]: value for i, value in enumerate(row)}
                    for row in data
                ]

        valid_data = []
        for row in processed_data:
            validated = {}
            for col in headers:
                value = row.get(col)
                if value is None:
                    raise ValueError(f"缺少必需列：{col}")
                validated[col] = f"{float(value):.4f}" if isinstance(value, (int, float)) else str(value)
            valid_data.append(validated)

        # 写入CSV文件
        with open(filename, 'w', newline='', encoding='UTF-8') as f_score:
            writer = csv.DictWriter(f_score, fieldnames=headers)
            writer.writeheader()
            writer.writerows(valid_data)

        return True

    except Exception as e:
        print(f"写入文件失败：{str(e)}，输出空文件")
        with open(filename, "w"):
            pass
        return False


