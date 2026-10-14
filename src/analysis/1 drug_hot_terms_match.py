"""
为用户推荐关键词
输入：药物名
"""
from pathlib import Path
import os
import re
import pandas as pd
import text_mining_all_abs as tm
from concurrent.futures import ProcessPoolExecutor, as_completed
import files_operation as fo
abs_cache = []


def init_worker(abs_list):
    """
    子进程初始化函数：将摘要列表缓存到每个worker进程中
    """
    global abs_cache
    abs_cache = abs_list


def normalize_id_list(id_list):
    """
    对id列表进行去空、去重、排序。
    """
    new_id_list = []
    for item in id_list:
        item_str = str(item).strip()
        if item_str:
            new_id_list.append(item_str)

    new_id_list = list(set(new_id_list))

    try:
        new_id_list = sorted(new_id_list, key=lambda x: int(x))
    except ValueError:
        new_id_list = sorted(new_id_list)

    return new_id_list


def chunk_list(data_list, chunk_size):
    """
    将大列表按固定大小切分。
    """
    for i in range(0, len(data_list), chunk_size):
        yield data_list[i:i + chunk_size]


def text_mining_abs(abs_line, kw):
    id_l = []
    hit_abs_l = []  # {idabs}
    id_num = "000000"
    for line in abs_line:
        # 字符串分割
        if "," in kw:
            new_line = re.split("[|:.;»]", line)
        else:
            new_line = re.split("[|:,.;»]", line)

        if "|||" in line:
            # 定义变量id_num
            id_num = new_line[0]

        for element in new_line:
            # 确保匹配的keyword完整
            if "," in kw:
                new_element = " " + element
                new_keyword = " " + kw
            else:
                new_element = " " + element + " "
                new_keyword = " " + kw + " "
            if new_keyword.lower() in new_element.lower():
                id_l.append(f"{id_num}")
                hit_abs_l.append(line)

                break
            else:
                continue

    return id_l, hit_abs_l


def process_keyword_batch(keyword_batch):
    """
    处理单个关键词：检查缓存、挖掘并保存结果。
    """
    batch_result_l = []

    for keyword in keyword_batch:
        kw_id_list, kw_hit_abs_list = text_mining_abs(abs_cache, keyword)
        if kw_id_list:
            # 将得到的id_list转换为集合去重，再转换回列表）
            kw_mined_list = normalize_id_list(kw_id_list)
            print(f"{keyword}共得到{len(kw_mined_list)}条PubMed_ID。")
        else:
            kw_mined_list = []

        batch_result_l.append((keyword, kw_mined_list))

    return batch_result_l


def main():
    drug = input("Drug Name(eg. Berberine):")
    input_drug_id_dir = "lexical library/Keyword_PMID_User/26ABS"
    drug_id_file = f"{drug}_mining_id.txt"
    drug_abs_file = f"{drug}_PM_ABS.txt"
    input_drug_id_path = os.path.join(input_drug_id_dir, drug_id_file)
    input_drug_abs_path = os.path.join(input_drug_id_dir, drug_abs_file)
    input_kw_path = "lexical library/Phenotype_names_250925_18231.txt"
    input_kw_list = fo.read_txt_file(input_kw_path)
    kw_dir_output = f"lexical library/{drug}_Related_Phenotype_PMID"
    fo.ensure_dir_exists(kw_dir_output)
    kws_count_d = {}

    # 对药物名进行文献挖掘
    if os.path.exists(input_drug_abs_path) and os.path.exists(input_drug_id_path):
        drug_id_l = normalize_id_list(fo.read_txt_file(input_drug_id_path))
        drug_pm_abs_l = fo.read_txt_file(input_drug_abs_path)
    else:
        drug_id_path, drug_id_l = tm.mining_kw_id(drug, drug_id_file,
                                              input_drug_id_dir)

        fo.ensure_dir_exists(input_drug_id_dir)
        drug_id_list = []
        drug_hit_abs_list = []

        for num in range(1, 1335):
            num1 = f'{num}'.rjust(4, '0')
            path_ab_file = f"ABS/pubmed26n{num1}.hao"
            path_ab = Path(path_ab_file)
            contents_ab = path_ab.read_text(encoding='UTF-8')
            line_list_f = contents_ab.splitlines()
            drug_id, drug_hit_abs = text_mining_abs(line_list_f, drug)
            drug_id_list.extend(drug_id)
            drug_hit_abs_list.extend(drug_hit_abs)


        # 输出信息1：将该关键词得到的ID输出为txt文件，文件名为关键词
        id_file_path = os.path.join(input_drug_id_dir, drug_id_file)
        # 输出信息2：将该关键词得到的摘要输出为txt文件，文件名为关键词
        id_abs_path = os.path.join(input_drug_id_dir, drug_abs_file)

        if drug_id_list:
            id_list_new = normalize_id_list(drug_id_list)
            mined_list = id_list_new
            print(f"{drug}共得到{len(id_list_new)}条PubMed_ID。")
        else:
            mined_list = []
        fo.txt_output(id_file_path, mined_list)
        fo.txt_output(id_abs_path, drug_hit_abs_list)
        drug_pm_abs_l = drug_hit_abs_list

    if not input_kw_list:
        return

    pending_kw_list = []
    for keyword in input_kw_list:
        kw_output_file = f"{keyword}_PMID.txt"
        output_path = os.path.join(kw_dir_output, kw_output_file)

        if os.path.exists(output_path):
            kw_mined_list = normalize_id_list(fo.read_txt_file(output_path))
            kws_count_d[keyword] = kw_mined_list

        else:
            pending_kw_list.append(keyword)

    # 设置并行进程
    max_workers = min(os.cpu_count() or 1, 8)
    batch_size = 20

    if pending_kw_list:
        keyword_batch_l = list(chunk_list(pending_kw_list, batch_size))
        with ProcessPoolExecutor(
            max_workers=max_workers,
            initializer=init_worker,
            initargs=(drug_pm_abs_l,)
        ) as executor:
            # 提交所有任务
            future_to_batch = {
                executor.submit(process_keyword_batch, keyword_batch): keyword_batch
                for keyword_batch in keyword_batch_l
            }

            for i, future in enumerate(as_completed(future_to_batch), start=1):
                keyword_batch = future_to_batch[future]
                try:
                    batch_result_l = future.result()

                    for keyword, kw_mined_list in batch_result_l:
                        kw_output_file = f"{keyword}_PMID.txt"
                        output_path = os.path.join(kw_dir_output, kw_output_file)
                        fo.txt_output(output_path, kw_mined_list)
                        kws_count_d[keyword] = kw_mined_list  # {keyword: report_id}

                    if i % 1 == 0 or i == len(future_to_batch):
                        print(
                            f"已完成 {i}/{len(future_to_batch)} 个批次"
                        )
                except Exception as e:
                    print(f"处理关键词批次{keyword_batch}时发生异常:{e}")

    # 所有关键词检索结束
    keyword_hit_rank = sorted(
        [(keyword, len(report_id)) for keyword, report_id in kws_count_d.items()],
        key=lambda x: x[1],
        reverse=True
    )
    # 将排序结果输出至Excel文件
    keyword_hit_df = pd.DataFrame(keyword_hit_rank, columns=["Keyword", "Count"])
    keyword_hit_rank_path = os.path.join(kw_dir_output, f"{drug}_Report_Keyword.xlsx")
    keyword_hit_df.to_excel(keyword_hit_rank_path, index=False, engine="openpyxl")
    print(f"关键词报道数量排序结果已输出至Excel文件：{keyword_hit_rank_path}")


if __name__ == "__main__":
    main()