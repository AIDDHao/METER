"""
目的：用转录组数据缩小基因范围
输入关键词列表
输入基因名
进行文献搜索：
    判断关键词是否在词汇库里，若在，就调用对应的文件结果，不在就从头搜索
输出：
    OTTER-AI step1、2
"""


import time
import os
import re
import glob
import venn_groups
import shutil
from functools import reduce
import heatmap1_data_to_csv
import pandas as pd
import heatmap2_groups
import files_operation as fo
import text_mining_all_abs as tm
import ppi1_match as pm
import ppi2_echarts as pe
import score_functions as sf


# 定义函数，调用词汇库中的txt文件，并匹配gene信息
def match_gene_ratio(ratio_info, gene_user, kw, mode1=True):
    """
    处理ratio文件信息
    mode1=ture:模式一，对有报道的表型进行挖掘。mode1=false：模式二，靶标几乎无报道
    输出ratio的信息文件
    """

    gene_name_info = []
    gene_info_all_dict = {}
    # 读取含有所有gene报道总量的txt文件，处理信息
    gene_re_num_info_list = fo.read_txt_file('Symbol_PMID_Count_2026.txt')  # 26年数据
    gene_total_num_info_d = {}
    for line_info in gene_re_num_info_list:
        gene_name = line_info.split("   ")[0]
        gene_total_value = line_info.split('   ')[-1]
        gene_total_num_info_d[gene_name] = gene_total_value  # 写入字典

    for i in range(len(ratio_info)):
        parts = re.split("[ \t]", ratio_info[i])
        parts = list(filter(None, parts))
        ratio_info[i] = parts
        if parts:
            gene_name_info.append(parts[0])
            gene_info_all_dict[parts[0]] = parts[1:]

    # 匹配用户输入的DEGs和词汇库信息
    gene_name_info_set = set(gene.lower() for gene in gene_name_info)
    gene_user_set = set(gene.lower() for gene in gene_user)
    gene_name_match_lower = set.intersection(gene_user_set, gene_name_info_set)
    gene_name_match = [gene for gene in gene_user if gene.lower() in gene_name_match_lower]
    # 创建输出数据。
    ratio_results_d = {}
    for gene_name_u in gene_user:
        if gene_name_u in gene_name_match:
            ratio_results_d[gene_name_u] = list(map(float, gene_info_all_dict[gene_name_u]))

        else:
            if gene_name_u in gene_info_all_dict.keys():
                non_re_num = 0
                gene_total_r_num = int(gene_total_num_info_d[gene_name_u])
                if gene_total_r_num != 0:
                    non_re_ratio = non_re_num / gene_total_r_num
                else:
                    non_re_ratio = 0
                non_re_gene_ratio_info_l = [non_re_num, non_re_ratio, non_re_ratio]
                non_re_gene_ratio_info_l = list(map(float, non_re_gene_ratio_info_l))
                ratio_results_d[gene_name_u] = non_re_gene_ratio_info_l
            else:
                pass

    # 输出DEGs的相关信息
    ratio_results_d_sorted = dict(sorted(ratio_results_d.items(), key=lambda item: item[1][-1], reverse=True))
    dir_results = f"{dir_output}/all results"
    os.makedirs(dir_results, exist_ok=True)
    if mode1:
        ratio_out_filename = f"{dir_results}/{kw}_Freq_ratio.txt"
        fo.txt_output(ratio_out_filename, ratio_results_d_sorted)
    else:
        if ratio_results_d_sorted:
            ratio_out_filename = f"{dir_results}/{kw}_Freq_PPI_ratio.txt"
            fo.txt_output(ratio_out_filename, ratio_results_d_sorted, out_mode="a")
    # 返回已报道的gene信息字典,筛选报道值
    filtered_results = {}
    for f_name, message in ratio_results_d_sorted.items():
        if int(message[0]) >= 1:
            filtered_results[f_name] = message
    gene_filter_l = list(filtered_results.keys())
    gene_non_l = list(set(gene_user).difference(set(gene_filter_l)))
    return filtered_results, gene_non_l

# 定义函数，调用词汇库中的txt文件，并匹配gene对应的文献报道信息
def match_gene_report(ratio_info, gene_user, kw, out_mum=True):
    """
    处理ratio文件信息
    匹配用户输入的DEGs和词汇库信息
    匹配不到的DEGs即为与该关键词无报道，创建对应数据后也收集起来
    无报道数据包括：报道0，总报道数，ratio（应为0） 类型：列表
    输出ratio的信息文件
    """
    gene_name_info = []
    gene_info_all_dict = {}
    # 读取含有所有gene报道总量的txt文件，处理信息
    gene_re_num_info_list = fo.read_txt_file('Symbol_PMID_Count_2026.txt')  # 26年摘要
    gene_total_num_info_d = {}
    for line_info in gene_re_num_info_list:
        gene_name = line_info.split("   ")[0]
        gene_total_value = line_info.split('   ')[-1]
        gene_total_num_info_d[gene_name] = gene_total_value  # 写入字典

    for i in range(len(ratio_info)):
        parts = re.split("[ \t]", ratio_info[i])
        parts = list(filter(None, parts))
        ratio_info[i] = parts
        if parts:
            gene_name_info.append(parts[0])
            gene_info_all_dict[parts[0]] = parts[1:]

    # 匹配用户输入的DEGs和词汇库信息
    gene_name_info_set = set(gene.lower() for gene in gene_name_info)
    gene_user_set = set(gene.lower() for gene in gene_user)
    gene_name_match_lower = set.intersection(gene_user_set, gene_name_info_set)
    gene_name_match = [gene for gene in gene_user if gene.lower() in gene_name_match_lower]
    # 创建输出数据。
    report_num_d = {}
    non_re_g_l = []
    for gene_name_u in gene_user:
        if gene_name_u in gene_name_match:
            gene_report_num = int(gene_info_all_dict[gene_name_u][0])
            report_num_d[gene_name_u] = gene_report_num
            if gene_report_num == 0:
                non_re_g_l.append(gene_name_u)

        else:
            if gene_name_u in gene_info_all_dict.keys():
                non_re_num = 0
                report_num_d[gene_name_u] = non_re_num
                non_re_g_l.append(gene_name_u)
            else:
                non_re_num = 0
                report_num_d[gene_name_u] = non_re_num
                non_re_g_l.append(gene_name_u)

    # 输出DEGs的相关信息
    ratio_result_d = sf.score_ratio(report_num_d)
    sorted_dict = sorted(ratio_result_d.items(), key=lambda item: item[1][-1], reverse=True)
    # 将排序后的字典项转换回字典
    ratio_results_d_sorted = dict(sorted_dict)
    if out_mum:
        dir_results = f"{dir_output}/all results"
        fo.ensure_dir_exists(dir_results)
        fo.txt_output(f"{dir_results}/{kw}_Freq_ratio.txt", ratio_results_d_sorted)
    # 返回已报道的gene信息字典,筛选报道值
    filtered_results = {}
    for f_name, message in ratio_results_d_sorted.items():
        if int(message[0]) >= 3:
            filtered_results[f_name] = message

    gene_filter_l = list(filtered_results.keys())
    gene_non_l = list(set(gene_user).difference(set(gene_filter_l)))
    return filtered_results, gene_non_l


# 定义函数，输出筛选过的检索结果
def output_f_mining_results(kw, f_m_r_d):
    """
    kw:keyword
    f_m_r_d:筛选过的检索结果，字典——{gene：报道数量，pubmed数量，ratio}
    """

    sorted_dict = dict(sorted(f_m_r_d.items(), key=lambda item: item[1][-1], reverse=True))
    # 输出
    dir_results = f"{dir_output}/results filtered"
    fo.ensure_dir_exists(dir_results)
    fo.txt_output(f"{dir_results}/{kw}_Freq_ratio.txt", sorted_dict)

def keyword_output_name(kw_input):
    """
    将关键词中的空格全部转换成字符串，便于输出文件与信息处理
    """
    kw_input = str(kw_input)
    if " " in kw_input:
        kw_output = kw_input.replace(" ", "_")
    else:
        kw_output = kw_input

    return kw_output


start_time = time.time()
# 读取keyword文件：用户输入
keyword_list = fo.read_txt_file(input("the keyword file path(eg:keyword.txt):"))
keyword_list = [s.lower() for s in keyword_list]
keyword_name_d = {}
for keyword_input in keyword_list:
    keyword_output = keyword_output_name(keyword_input)
    keyword_name_d[keyword_input] = keyword_output
print(f"User input: {len(keyword_list)} keywords.")
if len(keyword_list) == 1:
    kw_out_put = keyword_output_name(keyword_list[0])
    dir_output = f"OTTER_V2_{kw_out_put}_DE"
else:
    dir_output = f"OTTER_V2_{len(keyword_list)}kws_DE"
fo.ensure_dir_exists(dir_output)
step2 = True
output_len_num = 35

# 读取DEGs信息：用户输入
gene_list = fo.read_txt_file(input("the DEGs file path(eg:genes_id.txt):"))
print(f"User input: {len(gene_list)} genes.")

mode1_begin_time = time.time()
print(f"mode1开始————————————————————————————————————")
# 读取词汇库数据
disease_list = fo.read_txt_file("lexical library/Disease_names_all250227.txt")
disease_list = [s.lower() for s in disease_list]
phenotype_list = fo.read_txt_file("lexical library/Phenotype names_all250214.txt")
phenotype_list = [s.lower() for s in phenotype_list]
d_keyword_list = []
p_keyword_list = []
nex_keyword_list = []

# 判断关键词是否在词汇库里
k_g_d = {}  # 新字典——{keyword:genes}
k_m_d = {}  # 新字典——{k：message}
k_n_g_d = {}  # 新字典——{k:not report gene}
for keyword in keyword_list:
    # 在词汇库里，直接调用相关信息
    if keyword in disease_list:  # 疾病术语库
        d_keyword_list.append(keyword)
        kw_dir_path = f"lexical library/genes_to_ID250226_ratio_all_d"
        kw_path = os.path.join(kw_dir_path, f"{keyword}.txt")
        if not os.path.exists(kw_path):
            nex_keyword_list.append(keyword)
            continue
        kw_genes_ratio_info = fo.read_txt_file(kw_path)
        kw_res, gene_non_list = match_gene_report(kw_genes_ratio_info, gene_list, keyword)
        res_genes_list = list(kw_res.keys())
        k_g_d[keyword] = res_genes_list
        k_m_d[keyword] = kw_res
        k_n_g_d[keyword] = gene_non_list
        print(f"{keyword}:Input DEGs:{len(gene_list)}, "
              f"valid report genes: {len(res_genes_list)}, no reported genes: {len(gene_non_list)}")
        # 对该关键词检索完毕，整理并输出信息
        output_f_mining_results(keyword, kw_res)
        if keyword in phenotype_list:
            p_keyword_list.append(keyword)

    elif keyword in phenotype_list:  # 若不在疾病术语中，看是否在表型术语中
        p_keyword_list.append(keyword)
        kw_dir_path = f"lexical library/genes_to_ID1012_ratio_all_p"
        fo.ensure_dir_exists(kw_dir_path)
        fo.rename_files_to_lowercase(kw_dir_path)
        kw_path = os.path.join(kw_dir_path, f"{keyword}.txt")
        if not os.path.exists(kw_path):
            nex_keyword_list.append(keyword)
            continue
        kw_genes_ratio_info = fo.read_txt_file(kw_path)
        kw_res, gene_non_list = match_gene_report(kw_genes_ratio_info, gene_list, keyword)
        res_genes_list = list(kw_res.keys())
        k_g_d[keyword] = res_genes_list
        k_m_d[keyword] = kw_res
        k_n_g_d[keyword] = gene_non_list
        print(f"{keyword}:Input DEGs:{len(gene_list)}, "
              f"valid report genes: {len(res_genes_list)}, no reported genes: {len(gene_non_list)}")
        # 对该关键词检索完毕
        output_f_mining_results(keyword, kw_res)

    else:  # 不在词汇库里
        nex_keyword_list.append(keyword)

# -----------------------------------------------------------------------分节

if nex_keyword_list:
    # 设置提示信息，对于不在词汇库里的关键词，让用户选择系统自动进行文本挖掘或自行输入或直接跳过
    print(f"Press '1' means the system will perform text mining automatically.\n"
          f"Press '2' means you will upload your own Pubmed ID associated with the keyword.\n"
          f"Press '3' means the system will pass the keyword.\n"
          f"Please make your decision on these keywords in turn:")
    # 新建一个列表，建立需要从头搜索的关键词列表
    mining_kw_list = []
    for nex_keyword in nex_keyword_list:
        while True:
            user_choice = input(f"Please input your choice for this keyword:{nex_keyword}(Press'1', '2', '3')：")
            if user_choice in ['1', '2', '3']:
                if user_choice == '1':
                    print(1)
                    # 将该关键词新增进列表
                    mining_kw_list.append(nex_keyword)
                elif user_choice == '2':
                    print(2)
                elif user_choice == '3':
                    print(3)
                    print(f"The following keyword will be passed: {nex_keyword}")
                break
            else:
                print("不满足条件，请重新输入您的选择。")

    # -----------------------------------------分节
    # 若有需要从头检索的关键词，执行从头检索的代码。
    if mining_kw_list:
        for m_keyword in mining_kw_list:
            print("...." * 10)
            kw_r_dir = f"{dir_output}/{m_keyword}_output_info"
            # 调用函数，得到该关键词文本挖掘相关的id
            file_kw_collection = f"lexical library/Keyword_PMID_User/{m_keyword}_mining_id.txt"
            if os.path.exists(file_kw_collection):
                print(f"Keyword:{m_keyword} already has id information.")
                kw_id_path = file_kw_collection
                kw_id_l = fo.read_txt_file(file_kw_collection)
            else:
                kw_id_path, kw_id_l = tm.mining_kw_id_nlp(m_keyword, f"{m_keyword}_mining_id.txt",
                                                      kw_r_dir)

                # 得到一个新的关键词报道文献检索结果，将其储存至内置数据库
                if kw_id_l:
                    kw_new_id_info_dir = "lexical library/Keyword_PMID_User"
                    if not os.path.exists(kw_new_id_info_dir):
                        os.makedirs(kw_new_id_info_dir)
                    kw_new_id_info_path = f"{kw_new_id_info_dir}/{m_keyword}_mining_id.txt"
                    shutil.copy2(kw_id_path, kw_new_id_info_path)

            # 调用函数，得到该关键词文本挖掘相关的id对应的gene信息
            file_kw_id_collection = f"lexical library/ID_to_DE_User/{m_keyword}_mining_gene.txt"
            if os.path.exists(file_kw_id_collection):
                print(f"Keyword:{m_keyword} already has id_to_gene information。")
                kw_gene_path = file_kw_id_collection
                kw_gene_l = fo.read_txt_file(kw_gene_path)
            else:
                kw_gene_l, kw_gene_path = tm.id_to_gene(kw_id_path, f"{m_keyword}_mining_gene.txt",
                                                        kw_r_dir)
                # 得到一个新的关键词报道文献检索对应的gene结果，将其储存至内置数据库
                if kw_gene_l:
                    collect_genes_info_dirs = "lexical library/ID_to_DE_User"
                    if not os.path.exists(collect_genes_info_dirs):
                        os.makedirs(collect_genes_info_dirs)
                    kw_gene_new_info_path = f"{collect_genes_info_dirs}/{m_keyword}_mining_gene.txt"
                    shutil.copy2(kw_gene_path, kw_gene_new_info_path)

            # _____
            # 调用函数，将该关键词文本挖掘相关的id对应的gene信息生成对应的ratio
            file_kw_ratio_data = f"lexical library/DE_Ratio_User/{m_keyword}_mining_ratio.txt"
            if os.path.exists(file_kw_ratio_data):
                print(f"Keyword:{m_keyword} already has gene ratio information。")
                kw_ratio_path = file_kw_ratio_data
                kw_ratio_list = fo.read_txt_file(kw_gene_path)
                kw_ratio_d = fo.ratio_l_to_d(kw_ratio_list)

            else:  # 挖掘新的信息，并将结果保存至内置数据库
                kw_ratio_d, kw_ratio_path = tm.gene_to_ratio(m_keyword, kw_gene_path,
                                                             kw_id_path,
                                                             f"{m_keyword}_mining_ratio.txt", kw_r_dir)
                if kw_ratio_d:
                    collect_genes_ratio_dirs = "lexical library/DE_Ratio_User"
                    if not os.path.exists(collect_genes_ratio_dirs):
                        os.makedirs(collect_genes_ratio_dirs)
                    kw_ratio_new_info_path = f"{collect_genes_ratio_dirs}/{m_keyword}_mining_ratio.txt"
                    shutil.copy2(kw_ratio_path, kw_ratio_new_info_path)
            kw_ratio_info = fo.read_txt_file(kw_ratio_path)
            # 匹配用户输入的DEGs和词汇库信息
            kw_res, gene_non_list = match_gene_report(kw_ratio_info, gene_list, m_keyword)
            res_genes_list = list(kw_res.keys())
            if res_genes_list:
                k_g_d[m_keyword] = res_genes_list
                k_m_d[m_keyword] = kw_res
                k_n_g_d[m_keyword] = gene_non_list
                print(f"{m_keyword}:Input DEGs:{len(gene_list)}, "
                      f"valid report genes: {len(res_genes_list)}, no reported genes: {len(gene_non_list)}")
            # 对该关键词检索完毕，整理并输出信息
            output_f_mining_results(m_keyword, kw_res)

            now_time = time.time()
            run_time_now = now_time - start_time
            print(f"{m_keyword}检索结束，目前运行时间：{run_time_now}s")


# -----------------------------------------分节
# -----------------------------------------分节
if k_g_d:
    path_dir_re = dir_output
    # 判断两个列表都不为空,绘制交集venn
    if d_keyword_list and p_keyword_list:
        d_keyword_set = set(d_keyword_list)
        p_keyword_set = set(p_keyword_list)
        d_p_kw_i = set.intersection(d_keyword_set, p_keyword_set)
        if d_p_kw_i:
            venn_groups.venn_2g(d_keyword_set, p_keyword_set, ["disease terms", "phenotype terms"], path_dir_re)
            # 将交集输出为txt文件
            fo.txt_output(f"{path_dir_re}/current terms intersection.txt", d_p_kw_i)

    # -----------------------------------------分节

    # 输出检索结果
    fo.json_output(f"{path_dir_re}/all keywords match results.json", k_g_d)
    result_kws = []
    re_values = []
    t_score_d = {}
    for k, v in k_g_d.items():
        k_o = keyword_name_d[k.lower()]

        path_dir_score = f"{path_dir_re}/Score"
        os.makedirs(path_dir_score, exist_ok=True)
        try:
            if len(v) > 100:
                g_for_ppi_l = v[:100]

            else:
                g_for_ppi_l = v

        except(TypeError, IndexError) as e:
            g_for_ppi_l = []

        except Exception as e:
            g_for_ppi_l = []

        # 如果该关键词的g for ppi list不为空才进行后续操作
        if g_for_ppi_l:
            # 为每个基因计算PPI分数
            gene_ppi_d, ppi_gg_info_d = pm.ppi_gene_match(g_for_ppi_l)
            normalized_ppi_dict, normalized_ppi_info_dict = pm.normalize_scores(gene_ppi_d)
            print(f"normalized_ppi_dict:{normalized_ppi_dict}.normalized_ppi_info_dict:{normalized_ppi_info_dict}")
            # 输出ppi得分文件和归一化分数文件
            fo.txt_output(f"{path_dir_score}/{k_o}_Top100_DE_PPI_Score.txt", gene_ppi_d)
            fo.json_output(f"{path_dir_score}/{k_o}_Top100_DE_PPI_Score.json", gene_ppi_d)
            fo.json_output(f"{path_dir_score}/{k_o}_Top100_DE_PPI_Score_N.json", normalized_ppi_info_dict)

            # 为该关键词的文本分数进行归一化处理
            text_score_info_d = k_m_d[k]
            gene_text_d = {}
            for gene, text_score in text_score_info_d.items():
                if gene in g_for_ppi_l:
                    gene_text_d[gene] = text_score
            # 检查文本分数与ppi分数的字典长度应一致
            if len(gene_ppi_d) == len(gene_text_d):
                # 对文本分数进行归一化处理
                normalized_text_dict, normalized_text_info_d = sf.normalize_scores1(gene_text_d)

            else:
                normalized_text_dict = {}
                normalized_text_info_d = {}

            # 输出文本分数和归一化分数文件
            fo.json_output(f"{path_dir_score}/{k_o}_Top100_DE_Text_Score.json", gene_text_d)
            fo.json_output(f"{path_dir_score}/{k_o}_Top100_DE_Text_Score_N.json", normalized_text_info_d)

            # 计算ppi分数与文本分数的总和
            g_t_score_d = {}
            for gene, score_t in normalized_text_dict.items():
                score_d = {}
                t_s_n = score_t
                score_d["Text_Score"] = t_s_n
                p_s_n = normalized_ppi_dict[gene]
                score_d["PPI_Score"] = p_s_n
                total_s = t_s_n + p_s_n
                score_d["Total_Score"] = total_s
                g_t_score_d[gene] = score_d

            # 得到所有基因的分数信息
            g_t_score_d_s = dict(sorted(g_t_score_d.items(), key=lambda item: item[1]["Total_Score"], reverse=True))
            # 输出csv文件
            path_kw_total_score = f"{path_dir_score}/{k_o}_Top100_DE_Score.csv"
            sf.scores_to_csv(g_t_score_d_s, path_kw_total_score)
            t_score_d[k] = g_t_score_d_s

            # 可视化处理分节————————————————————————————————————————————————————————————————————————
            # 将该关键词的分数结果绘制成热图
            gene_l_s = list(g_t_score_d_s.keys())
            headers_score = ['Text_Score', 'PPI_Score', 'Total_Score']
            heatmap2_groups.heatmap_scores(
                path_csv=path_kw_total_score, col_data=headers_score, index_data= gene_l_s,
                path_dir=path_dir_score,
                output_fn=f"{k_o}_DE_Score_Heatmap_Vlag.svg",
                len_num=output_len_num,
                min_number=0, max_number=2
            )

        # --------------------------------可视化分节
            # 绘制点图-为该关键词创建数据
            path_top_nodes = f"{path_dir_score}/{k_o}_Top100_DE_Nodes.json"
            path_top_ppi_score = f"{path_dir_score}/Freq_ECharts_{k_o}_Top100_DE.txt"
            ppi_nodes = pe.ppi_data_top_nodes(ppi_gg_info_d, path_csv=path_kw_total_score, output_json=path_top_nodes,
                                              output_score=path_top_ppi_score)
            path_nodes_pair = f"{path_dir_score}/Nodes_DE_Pair_{k_o}_Top100.txt"
            pe.ppi_nodes_pair(ppi_nodes, path_nodes_pair)


        # ------------------------------可视化分节
            if len(v) >= 10:
                result_kws.append(k)
                re_values.append(v)
            else:
                print(f"keyword：{k} 文本匹配的基因数量小于10，不做可视化处理。")
    # 输出关键词结果
    fo.txt_output(f"{path_dir_re}/keywords_results.txt", result_kws)
    num_keys = len(result_kws)

    # 计算交集
    if re_values:
        gene_intersection = list(reduce(set.intersection, map(lambda x: set(x), re_values)))

    else:
        gene_intersection = []

    # 判断如果有交集，准备数据输出可视化图形
    if gene_intersection:
        path_dir_i_score = f"{path_dir_re}/Score/Intersections"
        os.makedirs(path_dir_i_score, exist_ok=True)
        # 为交集基因计算ppi分数和文本分数总和
        path_text_score = f"{path_dir_i_score}/All_Kws_Gene_Ratio_Text_Data.csv"
        header_csv, genes_r_l = heatmap1_data_to_csv.text_data_to_csv(message=k_m_d, i_l=gene_intersection,
                                                                      k_g_mess=k_g_d, kw_l=result_kws,
                                                                      path_name=path_text_score)
        if len(genes_r_l) > 100:
            gene_inter = genes_r_l[:100]
            if len(gene_inter) != 100:
                print(f"Warning:gene inters number should be 100, but the actually number is:{len(gene_inter)}!")

        else:
            gene_inter = genes_r_l
        # 计算ppi分数
        i_ppi_d, i_ppi_gg_d = pm.ppi_gene_match(gene_inter)
        # 对ppi分数进行归一化
        i_ppi_n_dict, i_ppi_n_info_d = pm.normalize_scores(i_ppi_d)
        # 输出ppi得分文件和归一化分数文件
        fo.txt_output(f"{path_dir_i_score}/Inter_DE_PPI_Score.txt", i_ppi_d)
        fo.json_output(f"{path_dir_i_score}/Inter_DE_PPI_Score.json", i_ppi_d)
        fo.json_output(f"{path_dir_i_score}/Inter_DE_PPI_Score_N.json", i_ppi_n_info_d)


        # 为每个关键词计算文本分数——————————————————————————————————————
        kws_text_d = sf.get_csv_data_score(path_text_score, result_kws)
        k_i_total_d = {}
        for kw_name, text_info in kws_text_d.items():
            kw_name_o = keyword_name_d[kw_name.lower()]
            i_text_d = {
                gene: text_s
                for gene, text_s in text_info.items()
                if gene in gene_inter
            }

            # 对文本分数的数据进行归一化处理
            i_text_n_dict = sf.normalize_scores2(i_text_d)
            if len(i_text_n_dict) != len(i_ppi_n_dict):
                print(f"Warning: Len error! Please check 'i_text_n_dict'.")

            # 计算文本分数和PPI总分
            i_total_d = sf.scores_total(text_dict=i_text_n_dict, ppi_dict=i_ppi_n_dict)
            i_total_d1 = sorted(i_total_d.items(), key=lambda item: item[1]["Total_Score"], reverse=True)
            i_total_d_s = dict(i_total_d1)
            path_kw_total_s = f"{path_dir_i_score}/{kw_name_o}_Top100_DE_Score.csv"
            sf.scores_to_csv(i_total_d_s, path_kw_total_s)

            # 该关键词统计得分写入总字典
            k_i_total_d[kw_name] = i_total_d_s

        # 数据准备与可视化作图——————————————————————————————————————
        # 获取关键词的数量,判断有几个关键词的结果
        if num_keys == 2:
            venn_groups.venn_2g(data1=set(k_g_d[result_kws[0]]), data2=set(k_g_d[result_kws[1]]),
                                name_data=result_kws, path_dir=path_dir_re)
        elif num_keys == 3:
            venn_groups.venn_3g(set(k_g_d[result_kws[0]]), set(k_g_d[result_kws[1]]), set(k_g_d[result_kws[2]]),
                                name_data=result_kws,
                                path_dir=path_dir_re)
        elif num_keys == 4:
            venn_groups.venn_4g(set(k_g_d[result_kws[0]]), set(k_g_d[result_kws[1]]),
                                set(k_g_d[result_kws[2]]), set(k_g_d[result_kws[3]]),
                                name_data=result_kws,
                                path_dir=path_dir_re)
        elif num_keys == 5:
            venn_groups.venn_5g(set(k_g_d[result_kws[0]]), set(k_g_d[result_kws[1]]),
                                set(k_g_d[result_kws[2]]), set(k_g_d[result_kws[3]]),
                                set(k_g_d[result_kws[4]]),
                                name_data=result_kws,
                                path_dir=path_dir_re)
        elif num_keys == 6:
            venn_groups.venn_6g(set(k_g_d[result_kws[0]]), set(k_g_d[result_kws[1]]),
                                set(k_g_d[result_kws[2]]), set(k_g_d[result_kws[3]]),
                                set(k_g_d[result_kws[4]]), set(k_g_d[result_kws[5]]),
                                name_data=result_kws,
                                path_dir=path_dir_re)
        elif num_keys == 1:
            print("只有一组数据，不能绘制Venn图")
        else:
            print("Venn 图只支持最多6个集合")

        # -----------------------------------------分节

        # 准备数据绘制热图
        path_heatmap_data_csv = f'{dir_output}/All_Kws_Gene_Score_Data.csv'
        header_csv_i = heatmap1_data_to_csv.total_ranked_csv(k_i_total_d, path_heatmap_data_csv)
        # 输出新的csv文件，设置文件名
        path_heatmap_data_csv_n = f'{dir_output}/Heatmap_Data_Normalized.csv'
        max_num, min_num = heatmap1_data_to_csv.data_normalization(path_name_i=path_heatmap_data_csv,
                                                                   path_name_o=path_heatmap_data_csv_n)

        # 用归一化后的数据做热图
        len_number = output_len_num
        heatmap2_groups.heatmap_groups(
            path_csv=path_heatmap_data_csv_n, col_data=[header_csv_i[2:]], index_data=gene_intersection,
            path_dir=path_dir_re,
            max_number=1,
            min_number=0,
            len_num=len_number,
            output_fn=f"Heatmap_Echarts_Final_Top{len_number}.svg"
        )

        # 创建nodes图的数据
        path_i_nodes = f"{dir_output}/Nodes_Data_DE.json"
        path_i_ppi_score = f"{dir_output}/Freq_ECharts_Final_PPI.txt"
        ppi_nodes = pe.ppi_data_for_nodes(i_ppi_gg_d, path_csv=path_heatmap_data_csv_n, output_json=path_i_nodes,
                                          output_score=path_i_ppi_score)
        # 创建nodes对的txt输出
        path_nodes_pair = f"{dir_output}/Nodes_DE_Pair_Final.txt"
        pe.ppi_nodes_pair(ppi_nodes, path_nodes_pair)

        mode1_time = time.time()
        print(f"mode1运行结束，运行时长：{mode1_time - mode1_begin_time}s.")


# ——————————————————————————————————————————————分节。扩大候选靶标范围至ppi
if k_n_g_d and step2:
    k_gi_d = {}
    k_ni_d = {}
    ppi_ri_dir = f"{dir_output}/PPI"
    os.makedirs(ppi_ri_dir, exist_ok=True)
    for kw_n, gene_m_l in k_n_g_d.items():
        kw_no = keyword_name_d[kw_n.lower()]
        kw_ri_dir = f"{ppi_ri_dir}/{kw_n}_output_info"
        os.makedirs(kw_ri_dir, exist_ok=True)
        g_gi_d = {}
        g_ni_d = {}
        gene_i_d = pm.gene_ppi_mining(gene_m_l)  # {gene：interaction gene list}
        if not gene_i_d:
            print(f"{kw_n}:gene_i_d is empty!")
        for gene_m, gene_i_l in gene_i_d.items() if gene_i_d else []:
            if kw_n in disease_list:  # 疾病术语库
                kw_dir_path = f"lexical library/genes_to_ID250226_ratio_all_d"
                if not os.path.exists(kw_dir_path):
                    os.makedirs(kw_dir_path)
                kw_path = os.path.join(kw_dir_path, f"{kw_n}.txt")
                if not os.path.exists(kw_path):
                    kw_ratio_info_path = f"lexical library/DE_Ratio_User/{kw_n}_mining_ratio.txt"
                    if os.path.exists(kw_ratio_info_path):
                        shutil.copy2(kw_ratio_info_path, kw_path)
                kw_genes_ratio_info = fo.read_txt_file(kw_path)
            elif kw_n in phenotype_list:
                kw_dir_path = f"lexical library/genes_to_ID1012_ratio_all_p"
                if not os.path.exists(kw_dir_path):
                    os.makedirs(kw_dir_path)
                kw_path = os.path.join(kw_dir_path, f"{kw_n}.txt")
                if not os.path.exists(kw_path):
                    kw_ratio_info_path = f"lexical library/DE_Ratio_User/{kw_n}_mining_ratio.txt"
                    if os.path.exists(kw_ratio_info_path):
                        shutil.copy2(kw_ratio_info_path, kw_path)
                kw_genes_ratio_info = fo.read_txt_file(kw_path)


            else:
                kw_ratio_info_path = f"lexical library/DE_Ratio_User/{kw_n}_mining_ratio.txt"
                if os.path.exists(kw_ratio_info_path):
                    kw_genes_ratio_info = fo.read_txt_file(kw_ratio_info_path)
                else:
                    print(f"Warning: no gene ratio info related to the keyword:{kw_n}!")
                    kw_genes_ratio_info = []

            kw_res_i, gene_non_i_list = match_gene_ratio(kw_genes_ratio_info, gene_i_l, kw_n, mode1=False)
            res_genes_list = list(kw_res_i.keys())
            ig_renum_d = {k: v[0] for k, v in kw_res_i.items()}
            ig_renum_r = {k: [v[0], v[-1]] for k, v in kw_res_i.items()}
            ig_non_re_d = {gene_non_i : [0, 0] for gene_non_i in gene_non_i_list}
            g_ni_check_l = fo.check_duplicate_data(ig_renum_r.keys(), ig_non_re_d.keys())
            if not g_ni_check_l:
                ig_renum_r.update(ig_non_re_d)
            else:
                print(f"Warning:ig_renum_d与g_ni_non_d不应有交集，但出现交集：{g_ni_check_l}!")
                for gene, value in kw_res_i.items():
                    if gene in g_ni_check_l:
                        print({gene: value})

            g_gi_d[gene_m] = res_genes_list
            g_ni_d[gene_m] = ig_renum_r

        # 该关键词所有DEG的互作蛋白文献检索完毕，记录结果
        g_ni_score_d = pm.ppi_m_score(g_ni_d)
        g_ni_score_path = f"{kw_ri_dir}/PPI_DE_Score.csv"
        header_csv_ni = ["Symbol_Name", "PPI_Total", "PPI_Report", "PM_Num", "Score"]
        fo.d_to_csv(g_ni_score_d, g_ni_score_path, header_csv_ni)


        # 对得分进行归一化,并保存为新的csv文件
        g_ni_score_n_path = f"{kw_ri_dir}/PPI_DE_Score_N.csv"
        pm.csv_data_normalize(g_ni_score_path, g_ni_score_n_path, col_name="Score")
        # 添加条件筛选工作表
        pm.data_rank_sheet(g_ni_score_n_path)
        k_gi_d[kw_n] = g_gi_d
        k_ni_d[kw_n] = g_ni_d
    pi_mining_output_path = f"{ppi_ri_dir}/PPI_DE_PM_Num.json"
    fo.json_output(pi_mining_output_path, k_ni_d)

    # 将各关键词的DEG的rank前十互作蛋白信息按照是否报道可视化为图谱
    pe.render_all_keywords(root_dir=ppi_ri_dir, global_json=pi_mining_output_path,
                           cfg=pe.Cfg(dir_suffix="_output_info", top_partners=5, open_browser=False))

    pm.ppi_info_i(ppi_ri_dir)
    # 为单个关键词的PPI打分结果生成热图
    try:
        kws_dir = glob.glob(os.path.join(ppi_ri_dir, '*_output_info'))
        for kw_dir in kws_dir:
            keyword = os.path.basename(kw_dir).split('_output_info')[0]

            # 查找关键词文件夹中的CSV文件(归一化之后的结果）
            kw_genes_files = glob.glob(os.path.join(kw_dir, 'PPI_DE_Score_N.xlsx'))

            if kw_genes_files:
                for kw_genes_file in kw_genes_files:
                    df_excel_kw = pd.read_excel(kw_genes_file, sheet_name="rank")
                    if ('Score' not in df_excel_kw.columns or 'PPI_Total'not in df_excel_kw.columns or
                            'PPI_Report' not in df_excel_kw.columns):
                        print(f"数据文件缺少Score列或PPI_Total列或PPI_Report列: {kw_genes_file}")
                        continue
                    # 提取数据绘制热图
                    len_num = output_len_num
                    heatmap2_groups.heatmap_ppis_kw(path_data=kw_genes_file, path_dir=kw_dir,
                                     output_fn=f"Heatmap_PPI_{keyword}.svg",
                                     len_num=len_num, min_number=0, max_number= 1)


            else:
                continue

    except Exception as e:
        pass

    try:
        # 取交集方法优化：
        path_ppi_i = pm.ppi_rank_data_i(ppi_ri_dir)
        # 将交集结果绘制成热图
        len_num = output_len_num
        heatmap2_groups.heatmap_ppis(path_data=path_ppi_i, path_dir=ppi_ri_dir,
                                     output_fn="Heatmap_PPI_Intersetion.svg",
                                     len_num=len_num)

        pe.ppi_data_dp_nodes(ppi_ri_dir, len_data="")

    except Exception as e:
        print(f"PPI绘制热图中产生未知错误：{str(e)}")
        pass

# 对OTTER2.0的数据结果的互作蛋白进行可视化处理
nodes_cmd = pe.run_exe_cmds(dir_output, len_num=output_len_num)
# 将cmd输出为txt文件
fo.txt_output(f"{dir_output}/Nodes_Echarts_cmds.txt", nodes_cmd)


end_time = time.time()
run_time = end_time - start_time
print(f"运行时间：{run_time}s")


