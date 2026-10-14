"""
自定义模块-定义各种文件操作的函数
"""

from pathlib import Path
import csv
import json
import re
import os

# 定义函数读取文件
def read_txt_file(txt_name):
    """
    txt_name:file name
    return:文件内容
    """
    path_txt = Path(txt_name)
    contents_txt = path_txt.read_text(encoding='UTF-8')
    txt_list = contents_txt.splitlines()
    return txt_list


# 定义函数，读取JSON文件
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
def txt_output(txt_name_output, txt_contents, out_mode = "w", out_encode= "UTF-8"):
    """
    txt_contents: 可以用for循环遍历的内容
    txt_name_output：待输出的文件名
    txt_contents：待输出的内容
    """
    with open(txt_name_output, out_mode, encoding=out_encode) as f_txt_output:
        if isinstance(txt_contents, list):
            for txt_content in txt_contents:
                f_txt_output.write(f"{txt_content}\n")
        if isinstance(txt_contents, set or tuple):
            for txt_content in txt_contents:
                f_txt_output.write(f"{txt_content}\n")
        if isinstance(txt_contents, dict):
            if not txt_contents:
                first_v_len = 0
            else:
                first_value = next(iter(txt_contents.values()))
                if not isinstance(first_value, (list, tuple, str)):
                    first_v_len = 0

                else:
                    first_v_len = len(first_value)

            if first_v_len == 3:
                for k, (num, total_num, ratio) in txt_contents.items():
                    ratio = float(ratio)
                    # 以制表符分隔，并写入文件
                    f_txt_output.write(f"{k} \t{num} \t{total_num} \t{ratio:.4f}\n")
            elif first_v_len == 2:
                for k, (num, ratio) in txt_contents.items():
                    ratio = float(ratio)
                    f_txt_output.write(f"{k} \t{num} \t{ratio:.4f}\n")
            else:
                pass


# 定义函数，确保指定的文件夹存在
def ensure_dir_exists(dir_path):
    Path(dir_path).mkdir(parents=True, exist_ok=True)


# 定义函数，输出一个空文件
def emtpy_file_output(file_name_empty):
    """
    新建一个空文件
    file_name_empty：文件名
    """

    with open(file_name_empty, 'w'):
        pass


# 定义函数，文件名改为小写
def rename_files_to_lowercase(folder_path):
    """
    folder_path:文件夹路径
    """
    for filename in os.listdir(folder_path):
        # 文件路径
        old_path = os.path.join(folder_path, filename)

        if os.path.isfile(old_path):
            # 检查文件名
            if filename != filename.lower():
                # 转换文件名
                new_filename = filename.lower()
                new_path = os.path.join(folder_path, new_filename)

                # 如果新路径与旧路径不同，则重命名文件
                if old_path != new_path:
                    os.rename(old_path, new_path)


# 定义函数，将json文件转化输出为csv文件
def json_to_csv(json_name, output_name, header_info):
    """
    将基因评分数据写入CSV文件，支持多种输入格式

    参数：
    json_name:输入的json文件名
    output_name:输出的csv文件名
    header_info:标题信息,可以是列表，元组这种有序的格式
    """
    # 统一数据格式
    processed_data = []
    headers = list(header_info)
    with open(json_name, "r", encoding="UTF-8")  as f_j:
        data = json.load(f_j)


    try:
        if isinstance(data, dict):
            if all(isinstance(v, (list, tuple)) for v in data.values()):
                processed_data = [
                    {col: values[i] for col, values in data.items()}
                    for i in range(max(len(v) for v in data.values()))
                ]
            else:
                processed_data = [
                    {headers[0]: k, **v} if isinstance(v, dict) else {headers[0]: k, headers[1]: v}
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

        # 验证数据完整性
        valid_data = []
        for row in processed_data:
            validated = {}
            for col in headers:
                value = row.get(col)
                if value is None:
                    raise ValueError(f"缺少必需列：{col}")
                # 统一数值格式
                validated[col] = f"{float(value):.4f}" if isinstance(value, (int, float)) else str(value)
            valid_data.append(validated)

        # 写入CSV文件
        with open(output_name, 'w', newline='', encoding='UTF-8') as f_c:
            writer = csv.DictWriter(f_c, fieldnames=headers)
            writer.writeheader()
            writer.writerows(valid_data)

        return True

    except Exception as e:
        with open(json_name, "w"):
            pass
        return False

# 定义函数，将输入信息输出为csv文件
def d_to_csv(input_data, output_name, header_info):
    """
    将基因评分数据写入CSV文件，支持多种输入格式
    参数：
    input_data:输入的数据信息
    output_name:输出的csv文件名
    header_info:标题信息,header_info=["Symbol", "Num", "Score"]
    """

    processed_data = []

    try:
        if isinstance(input_data, dict):
            processed_data = []
            for gene, values in input_data.items():
                if isinstance(values, (list, tuple)) or len(values) == 2:
                    row = {
                        header_info[0]: gene,
                        header_info[1]: values[0],
                        header_info[2]: values[1]
                    }
                elif isinstance(values, dict) or len(values) == 4:
                    row = {
                        header_info[0]: gene,
                        header_info[1]: values["PPI_Total"],
                        header_info[2]: values["PPI_Report"],
                        header_info[3]: values["PM_Num"],
                        header_info[4]: values["Score"],
                    }
                else:
                    raise ValueError(f"基因 {gene} 对应的值格式错误，应为长度为2或4的列表/元组")
                processed_data.append(row)
        elif isinstance(input_data, list):
            pass

        # 验证数据完整性
        valid_data = []
        for row in processed_data:
            validated = {}
            for col in header_info:
                value = row.get(col)
                if value is None:
                    raise ValueError(f"缺少必需列：{col}")
                # 统一数值格式
                if col == 'Score' and isinstance(value, (int, float)):
                    validated[col] = f"{float(value):.4f}"
                else:
                    validated[col] = str(value)
            valid_data.append(validated)

        # 写入CSV文件
        with open(output_name, 'w', newline='', encoding='UTF-8') as f_c:
            writer = csv.DictWriter(f_c, fieldnames=header_info)
            writer.writeheader()
            writer.writerows(valid_data)

    except Exception as e:
        print(f"CSV生成失败：{str(e)}")


# 定义函数，将从文件中读取的gene的ratio信息（列表），转换为字典格式。eg.[ZNF592 	1 	30 	0.0333]
def ratio_l_to_d(ratio_file_list):
    ratio_dict = {}
    for ratio_line in ratio_file_list:
        # 去除首尾空白，然后用正则分割所有连续空白（包括空格和制表符）
        parts = re.split(r'\s+', ratio_line.strip())
        if len(parts) < 1:
            continue
        key = parts[0]
        values = parts[1:] if len(parts) > 1 else []
        ratio_dict[key] = values
    return ratio_dict

# 定义函数，检查两个数据集中的元素有没有重复
def check_duplicate_data(data1, data2):
    """
    data1&2:可以是列表、元组、集合、字典的键
    return：重复元素的列表
    """

    data1 = set(data1)
    data2 = set(data2)
    check_data = list(set.intersection(data1, data2))

    return check_data




