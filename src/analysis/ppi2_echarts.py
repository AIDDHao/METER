"""
为互作蛋白的节点数据进行可视化处理
1. 通用函数部分
2. 数据准备部分
3. 可视化处理部分
"""
from __future__ import annotations
import os
import subprocess
import glob
import pandas as pd
from pathlib import Path
import operator
import csv
import json, math, webbrowser
from dataclasses import dataclass
from typing import Dict, Any, List, Tuple
from collections import OrderedDict



# ————————————————————————————————————————————函数部分——————————————————————————————————————————————————
# ————————————————————————————————————————————函数部分——————————————————————————————————————————————————

# 定义函数读取文件
def read_txt_file(txt_name):
    """
    txt_name:file name
    return:将文件中的内容按行读进一个列表
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


# ————————————————————————————数据准备部分————————————————————————————————————————
# ————————————————————————————数据准备部分————————————————————————————————————————
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
                        print(f"CSV数据格式错误：{gene}    {str(e)}")
    except FileNotFoundError:
        return

    # 筛选有效基因
    valid_genes = [g for g in gene_names_l
                   if g in csv_data and g in gg_dict]

    sorted_genes = sorted(valid_genes,
                          key=lambda x: csv_data[x][0])

    # 构建字典
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

            # 首先收集所有行和分数
            for row in csv_content:
                all_rows.append(row)
                gene_name = row[0]
                try:

                    total_score = float(row[-1])
                    gene_scores.append((gene_name, total_score))
                except (ValueError, IndexError) as e:
                    print(f"分数解析错误：{gene_name}   {str(e)}")

            # 按总分降序排序
            gene_scores.sort(key=lambda x: x[1], reverse=True)

            # 排名
            ranking = {}
            current_rank = 1
            prev_score = None

            for i, (gene, score) in enumerate(gene_scores):
                # 处理相同分数的情况
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

        # 收集JSON所需数据
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
    ppi_i_path = os.path.join(ppi_dir, 'PPI_Rank_Data.xlsx')
    ppi_all = read_json_file(json_path)
    df_excel = pd.read_excel(ppi_i_path)
    ppi_i_data = df_excel.sort_values(by='Total_Score', ascending=False, inplace=False)

    # 生成文件
    dp_output_files(ppi_i_data, ppi_all, ppi_dir, data_len=len_data, kw=None)

    # 处理每个关键词对应的子文件夹
    kws_dir = glob.glob(os.path.join(ppi_dir, '*_output_info'))
    for kw_dir in kws_dir:
        keyword = os.path.basename(kw_dir).split('_output_info')[0]

        # 查找关键词文件夹中的CSV文件
        kw_genes_files = glob.glob(os.path.join(kw_dir, 'PPI_DE_Score_N.xlsx'))

        if kw_genes_files:
            for kw_genes_file in kw_genes_files:
                # 处理CSV文件
                df_excel_kw = pd.read_excel(kw_genes_file, sheet_name="rank")
                if 'Score' not in df_excel_kw.columns:
                    continue

                ppi_kw_data = df_excel_kw.sort_values(by='Score', ascending=False, inplace=False)
                dp_output_files(ppi_kw_data, ppi_all, kw_dir, data_len=len_data, kw=f"{keyword}")
        else:
            continue

def dp_output_files(df, ppi_data, output_dir,data_len="", kw=None):
    """
    dp:data ppi
    对象：单个关键词数据
    PPI板块：为节点图生成数据
    df：ppi板块的结果数据
    ppi_data:包含所有kws的元基因与互作基因信息
    output_dir:输出数据的文件夹路径
    data_len="":空表示读取所有基因，调用函数时可以指定数据长度，如50
    prefix="": 默认无前缀
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

    with open(os.path.join(output_dir, f"Nodes_DE_Pair{kw_out_name}_PPI.txt"), 'w') as f1:
        data_lines = []
        for gene in gene_de:
            gene = str(gene)
            found = False
            all_ppi_data = {}

            if not kw:  # kw=None
                for keyword, gene_dict in ppi_data.items():
                    if gene in gene_dict:
                        for ppi_gene, scores in gene_dict[gene].items():
                            data_line = f"{gene}   {ppi_gene}"
                            if data_line not in data_lines:
                                data_lines.append(data_line)

                for data_line in data_lines:
                    f1.write(f"{data_line}\n")

            elif kw and kw in ppi_data:
                # 只处理当前关键词
                gene_dict = ppi_data[kw]
                if gene in gene_dict:
                    for ppi_gene, scores in gene_dict[gene].items():
                        f1.write(f"{gene}   {ppi_gene}\n")


        path_txt_ppi = os.path.join(output_dir, f"Freq_ECharts{kw_out_name}_PPI.txt")
        df = df.fillna("")
        gene_col = df.columns[0]
        df[gene_col] = df[gene_col].astype(str).str.strip()
        data_de = df[df[gene_col].isin(gene_de)].copy()
        data_de = data_de.replace({"\r\n": " ", "\n": " ", "\r": " "}, regex=True)
        with open(path_txt_ppi, "w", encoding="UTF-8", newline="\n") as f:
            for row in data_de.itertuples(index=False, name=None):
                f.write("   ".join(map(str, row)) + "\n")


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

        for interaction in interaction_list:
            source = interaction["source"]
            target = interaction["target"]
            interactions.append((rank, source, target))

    # 按排名排序
    interactions.sort(key=operator.itemgetter(0))

    # 写入文件
    with open(output_file, 'w') as f:
        for rank, source, target in interactions:
            f.write(f"{source}   {target}\n")


# ————————————————————————————————————————3.可视化处理部分————————————————————————————————————————
# ————————————————————————————————————————3.可视化处理部分————————————————————————————————————————

# 定义工具函数：把路径中的空格改为下划线，并重命名
def spaces_clean(s: str) -> str:
    """
    清洗字符串中的空格
    """
    s = s.replace("\u3000", " ").strip()
    while " " in s:
        s = s.replace(" ", "_")
    return s


# 定义函数，将路径中的空格重命名为下划线
def normalize_path_name(path: Path, output_file: bool = False) -> Path:
    """
    output_file=True：输出文件路径；
    output_file=False：现有路径。
    """
    path = Path(path)
    parts = list(path.parts)
    if len(parts) == 0:
        return path

    root = Path(parts[0])
    path_cur  = root
    for i in range(1, len(parts)):
        old_part = parts[i]
        new_part = spaces_clean(old_part)
        old_path = path_cur / old_part
        new_path = path_cur / new_part
        is_last = (i == len(parts) - 1)
        if is_last and output_file:
            path_cur = new_path
            continue

        if old_path.exists():
            if old_part != new_part:

                if new_path.exists() and old_path.resolve() != new_path.resolve():
                    continue
                old_path.rename(new_path)
            path_cur = new_path
        else:

            path_cur = new_path

    return path_cur
def len_num_judge(ori_len, data_len):
    if ori_len <= data_len:
        len_num = ori_len
    else:
        len_num = data_len
    return len_num
def run_exe_cmds(dir_path, len_num=50):
    """
    dir_path:数据结果路径
    """
    # 定义文件路径
    exe_path = "OTTER_V2_ECharts_Generate.exe"
    exe_name = exe_path.split('.exe')[0]
    dir_i = dir_path
    dir_score = os.path.join(dir_path, "Score")
    dir_ppi = os.path.join(dir_path, "PPI")

    # 存储所有生成的cmd命令
    commands = []

    # 1. 处理交集数据
    freq_file = os.path.join(dir_i, "Freq_ECharts_Final_PPI.txt")
    freq_data_list = read_txt_file(freq_file)
    len_num = len_num_judge(len_num, len(freq_data_list))
    final_nodes = os.path.join(dir_i, "Nodes_DE_Pair_Final.txt")
    if os.path.exists(freq_file) and os.path.exists(final_nodes):
        output_html = os.path.join(dir_i, f"PPI_Echarts_Top_{len_num}.html")
        commands.append(
            f'"{exe_name}" {len_num} "{freq_file}" "{final_nodes}" "{output_html}"'
        )

    # 2. 处理Score数据目录下的keyword数据
    for freq_file in glob.glob(os.path.join(dir_score, "Freq_ECharts_*_Top*.txt")):

        if " " in freq_file:
            freq_file = freq_file.replace(" ", "_")
        len_num = len_num_judge(len_num, len(freq_data_list))
        # 提取keyword
        base_name = os.path.basename(freq_file)
        keyword = (base_name.split("Freq_ECharts_")[1]).split("_Top")[0]

        # 查找对应的Nodes文件
        nodes_files = glob.glob(os.path.join(dir_score, f"Nodes_DE_Pair_{keyword}_Top*.txt"))

        if nodes_files and len(nodes_files) == 1:
            if " " in nodes_files[0]:
                nodes_input_path = nodes_files[0].replace(" ", "_")
            else:
                nodes_input_path = nodes_files[0]
            output_html = os.path.join(
                dir_score,
                f"PPI_Echarts_{keyword}_Top_{len_num}.html"
            )
            if " " in output_html:
                output_html = output_html.replace(" ", "_")

            commands.append(
                f'"{exe_name}" {len_num} "{freq_file}" "{nodes_input_path}" "{output_html}"'
            )

    # 3. 处理ppi板块的数据目录下的keyword数据
    prefix_f = "Freq_ECharts"
    prefix_n = "Nodes_DE_Pair"
    suffix_ppi = "_PPI.txt"

    paths_ori = []
    # 遍历ppi板块路径
    for dirpath, _, filenames in os.walk(dir_ppi):
        dirpath = Path(dirpath)
        freq_kw = {}
        nodes_kw = {}

        for fn in filenames:
            # 符合条件的文件名
            base_fn = fn.strip()
            kw = ""

            if not base_fn.endswith(suffix_ppi):
                continue

            # 去掉后缀
            pre_file = base_fn.split(suffix_ppi)[0]

            if base_fn.startswith(prefix_f):
                # 去掉前缀
                mid_file = pre_file.split(prefix_f)[-1]
                if mid_file.startswith("_"):
                    mid_file = mid_file[1:]
                kw = mid_file.replace("\u3000", " ").strip()
                while "  " in kw:
                    kw = kw.replace("  ", " ")
                freq_kw[kw] = dirpath / base_fn

            elif base_fn.startswith(prefix_n):
                mid_file = pre_file.split(prefix_n)[1]
                if mid_file.startswith("_"):
                    mid_file = mid_file[1:]
                kw = mid_file.replace("\u3000", " ").strip()
                while "  " in kw:
                    kw = kw.replace("  ", " ")
                nodes_kw[kw] = dirpath / base_fn

            else:
                continue


        for kw in sorted(set(freq_kw) & set(nodes_kw)):
            freq_file = freq_kw[kw]
            nodes_file = nodes_kw[kw]

            label = f"Echarts_{kw}" if kw else "Echarts"
            output_html = dirpath / f"PPI_{label}_Top_{len_num}.html"
            paths_ori.append((exe_name, len_num, freq_file, nodes_file, output_html))

    for exe_path, L, freq_file, nodes_file, output_html in paths_ori:
        freq_file_new = normalize_path_name(Path(freq_file), output_file=False)
        nodes_file_new = normalize_path_name(Path(nodes_file), output_file=False)
        output_html_new = normalize_path_name(Path(output_html), output_file=True)

        # 生成命令
        cmd = f'"{exe_path}" {L} "{freq_file_new}" "{nodes_file_new}" "{output_html_new}"'
        commands.append(cmd)


    # 执行所有生成的命令
    for cmd in commands:
        print(f"执行命令: {cmd}")
        try:
            process = subprocess.Popen(
                cmd,
                shell=True,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            # 等待程序执行
            try:
                stdout, stderr = process.communicate(timeout=5)
            except subprocess.TimeoutExpired:
                process.stdin.write('\n')
                process.stdin.flush()
                stdout, stderr = process.communicate(timeout=20)

            # 检查执行结果
            if process.returncode == 0:
                print("命令执行成功")
                print(stdout)
            else:
                print(f"命令执行失败，错误代码: {process.returncode}")
                print(f"错误信息: {stderr}")
        except Exception as e:
            print(f"执行命令时发生异常: {str(e)}")
        finally:
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=5)
                except:
                    process.kill()
            print("\n" + "-" * 50 + "\n")

    return commands


# 为PPI板块单个关键词预测结果生成HTML网页
# ---------------- 1.config ----------------
@dataclass
class Cfg:
    dir_suffix: str = "_output_info"
    rank_excel: str = "PPI_DE_Score_N.xlsx"
    sheet: str = "rank"
    gene_col: str = "Symbol_Name"
    score_col: str = "Score"
    top_genes: int = 20
    top_partners: int = 10  # 互作蛋白数量
    open_browser: bool = False


# ---------------- HTML----------------
TPL = """<!doctype html><html><head><meta charset="utf-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1"/>
<title>{title}</title>
<script src="https://unpkg.com/cytoscape/dist/cytoscape.min.js"></script>
<script src="https://unpkg.com/cytoscape-svg/cytoscape-svg.js"></script>
<style>
body{{margin:0;font-family:Arial,Helvetica,sans-serif;background:#fff}}
.bar{{display:flex;gap:10px;align-items:center;padding:10px 14px;border-bottom:1px solid #eee}}
.kw{{font-weight:800;font-size:18px}}
select,button{{font-family:Arial,Helvetica,sans-serif;font-size:14px;padding:6px 10px}}
#cy{{width:100%;height:calc(100vh - 52px)}}
.sp{{flex:1}}
</style></head><body>
<div class="bar"><div class="kw">Keyword: {keyword}</div>
<div>Center:<select id="g"></select></div>
<button id="fit">Fit</button><button id="svg">SVG</button><button id="png">PNG</button><div class="sp"></div></div>
<div id="cy"></div>
<script>
const P={payload_json}, G=Object.keys(P), sel=document.getElementById("g");
G.forEach(x=>{{let o=document.createElement("option");o.value=x;o.textContent=x;sel.appendChild(o);}});
let cy=null;
function draw(x){{
  if(cy)cy.destroy();
  cy=cytoscape({{
    container:document.getElementById("cy"),
    elements:P[x].elements,
    layout:{{name:"preset"}},
    style:[
      {{selector:"node",style:{{
        "label":"data(label)",
        "font-family":"Arial,Helvetica,sans-serif",
        "font-weight":"bold",
        "text-valign":"center",
        "text-halign":"center",
        "text-wrap":"wrap",
        "text-max-width":90,
        "font-size":16,
        "width":"data(size)","height":"data(size)",
        "background-color":"data(color)",
        "border-color":"data(border)","border-width":3,
        "color":"#222"
      }}}},
      {{selector:'node[group="center"]',style:{{
        "shape":"round-rectangle",
        "text-max-width":160,
        "font-size":22
      }}}},
      {{selector:"edge",style:{{"width":4,"curve-style":"straight","line-color":"data(edgeColor)"}}}}
    ],
    userZoomingEnabled:true,userPanningEnabled:true
  }});
  cy.fit(undefined,70);
}}
sel.onchange=e=>draw(e.target.value);
document.getElementById("fit").onclick=()=>cy&&cy.fit(undefined,70);
document.getElementById("svg").onclick=()=>{{
  if(!cy)return;
  const s=cy.svg({{full:true,scale:1,bg:"white",useComputedStyles:true}});
  const b=new Blob([s],{{type:"image/svg+xml;charset=utf-8"}});
  const u=URL.createObjectURL(b);const a=document.createElement("a");
  a.href=u;a.download=sel.value+"_network.svg";a.click();URL.revokeObjectURL(u);
}};
document.getElementById("png").onclick=()=>{{
  if(!cy)return;
  const p=cy.png({{full:true,scale:6,bg:"white"}});
  const a=document.createElement("a");a.href=p;a.download=sel.value+"_network.png";a.click();
}};
if(G.length){{sel.value=G[0];draw(G[0]);}}
</script></body></html>"""


# ---------------- renderer ----------------
def _write_html(out_path: Path, *, keyword: str, payload_by_gene: Dict[str, Any]) -> None:
    html = TPL.format(
        title=f"{keyword} - Interacting Proteins",
        keyword=keyword,
        payload_json=json.dumps(payload_by_gene, ensure_ascii=False)
    )
    out_path.write_text(html, encoding="utf-8")


# ---------------- color ----------------
def _clamp01(x: float) -> float:
    return max(0.0, min(1.0, float(x)))

def _hex_to_rgb01(h: str) -> Tuple[float, float, float]:
    h = h.strip().lstrip("#")
    r = int(h[0:2], 16) / 255.0
    g = int(h[2:4], 16) / 255.0
    b = int(h[4:6], 16) / 255.0
    return r, g, b

def _rgb01_to_hex(rgb: Tuple[float, float, float]) -> str:
    r, g, b = rgb
    r = int(round(_clamp01(r) * 255))
    g = int(round(_clamp01(g) * 255))
    b = int(round(_clamp01(b) * 255))
    return f"#{r:02x}{g:02x}{b:02x}"

def _lerp(a: float, b: float, t: float) -> float:
    return a + (b - a) * t

def _color_lerp_hex(score01: float, hex_min: str, hex_max: str, gamma: float = 1.0) -> str:
    t = _clamp01(score01) ** float(gamma)
    mn = _hex_to_rgb01(hex_min)
    mx = _hex_to_rgb01(hex_max)
    rgb = (_lerp(mn[0], mx[0], t), _lerp(mn[1], mx[1], t), _lerp(mn[2], mx[2], t))
    return _rgb01_to_hex(rgb)


# ---------------- geometry ----------------
def _linemap(vals: List[float], a: float, b: float) -> List[float]:
    if not vals:
        return []
    vmin, vmax = min(vals), max(vals)
    if abs(vmax - vmin) < 1e-12:
        return [(a + b) / 2.0] * len(vals)
    return [a + (v - vmin) / (vmax - vmin) * (b - a) for v in vals]

def _semi(names: List[str], r: float, upper: bool, margin: float = 0.30) -> Dict[str, Tuple[float, float]]:
    if not names:
        return {}
    if upper:
        start, end = math.pi - margin, 0.0 + margin
    else:
        start, end = math.pi + margin, 2 * math.pi - margin
    angs = [(start + end) / 2] if len(names) == 1 else [
        start + (end - start) * i / (len(names) - 1) for i in range(len(names))
    ]
    return {n: (r * math.cos(a), -r * math.sin(a)) for n, a in zip(names, angs)}

def _top_genes(xlsx: Path, cfg: Cfg) -> List[str]:
    df = pd.read_excel(xlsx, sheet_name=cfg.sheet, engine="openpyxl")[[cfg.gene_col, cfg.score_col]].copy()
    df[cfg.gene_col] = df[cfg.gene_col].astype(str).str.strip()
    df[cfg.score_col] = pd.to_numeric(df[cfg.score_col], errors="coerce")
    df = df.dropna().drop_duplicates(subset=[cfg.gene_col], keep="first").sort_values(cfg.score_col, ascending=False)
    return df[cfg.gene_col].head(cfg.top_genes).tolist()

def _pick(gene_dict: Dict[str, Any], k: int):
    """
      reported: PM_Num > 0
      unreported: PM_Num == 0
    """
    reported, unreported = [], []

    for p, v in gene_dict.items():
        if not isinstance(v, (list, tuple)) or len(v) < 2:
            continue
        pm, sc = int(v[0]), float(v[1])
        if pm > 0:
            reported.append((p, pm, sc))
        else:
            unreported.append((p, pm, sc))

    # 排序规则不变
    reported.sort(key=lambda t: (t[1], t[2]), reverse=True)
    unreported.sort(key=lambda t: (t[1], t[2]), reverse=True)

    if len(reported) <= k:
        reported_sel = reported
    else:
        reported_sel = reported[:k]

    if len(unreported) <= k:
        unreported_sel = unreported
    else:
        unreported_sel = unreported[:k]

    return reported_sel, unreported_sel


# ---------------- payload ----------------
def build_payload(interactions_kw: Dict[str, Any], genes: List[str], k: int, base_r: float = 320.0, gamma: float = 1.0) -> Dict[str, Any]:
    CENTER_FILL, CENTER_BORDER = "#E1C1B4", "#C78E77"

    # reported
    REP_FILL_MIN, REP_FILL_MAX = "#F4DDDC", "#E4ACAB"
    REP_BRD_MIN,  REP_BRD_MAX  = "#E2A9A8", "#D17371"

    # unreported
    UN_FILL_FIXED = "#AFCBD7"
    UN_BRD_FIXED  = "#5B97B5"

    out: Dict[str, Any] = {}

    for g in genes:
        gdict = interactions_kw.get(g)
        if not isinstance(gdict, dict):
            continue

        rep, un = _pick(gdict, k)

        n = max(len(rep), len(un))
        r = base_r + max(0, n - 6) * 18

        # reported：top, unreported：bottom
        rep_names = [x[0] for x in rep]
        un_names  = [x[0] for x in un]
        rep_pos = _semi(rep_names, r, upper=True,  margin=0.30)
        un_pos  = _semi(un_names,  r * 0.92, upper=False, margin=0.30)

        # 由PM_Num决定节点大小
        rep_sizes = _linemap([math.log1p(x[1]) for x in rep], 30, 86)
        un_sizes  = _linemap([math.log1p(x[1]) for x in un], 26, 74)

        els: List[Dict[str, Any]] = [{
            "data": {
                "id": g, "label": g, "group": "center",
                "size": 96, "color": CENTER_FILL, "border": CENTER_BORDER,
                "edgeColor": CENTER_BORDER
            },
            "position": {"x": 0, "y": 0}
        }]

        # reported
        for (p, pm, sc), sz in zip(rep, rep_sizes):
            x, y = rep_pos.get(p, (0.0, 0.0))
            fill   = _color_lerp_hex(sc, REP_FILL_MIN, REP_FILL_MAX, gamma=gamma)
            border = _color_lerp_hex(sc, REP_BRD_MIN,  REP_BRD_MAX,  gamma=gamma)
            els.append({
                "data": {"id": p, "label": p, "group": "reported",
                         "size": float(sz), "color": fill, "border": border, "edgeColor": border},
                "position": {"x": float(x), "y": float(y)}
            })
            els.append({"data": {"id": f"{g}__{p}", "source": g, "target": p, "edgeColor": border}})

        # unreported
        for (p, pm, sc), sz in zip(un, un_sizes):
            x, y = un_pos.get(p, (0.0, 0.0))
            els.append({
                "data": {"id": p, "label": p, "group": "unreported",
                         "size": float(sz), "color": UN_FILL_FIXED, "border": UN_BRD_FIXED, "edgeColor": UN_BRD_FIXED},
                "position": {"x": float(x), "y": float(y)}
            })
            els.append({"data": {"id": f"{g}__{p}", "source": g, "target": p, "edgeColor": UN_BRD_FIXED}})

        out[g] = {"center": g, "elements": els}

    return out

def _genes_all_from_sheet(xlsx: Path, sheet: str, gene_col: str) -> List[str]:
    df = pd.read_excel(xlsx, sheet_name=sheet, engine="openpyxl")[[gene_col]].copy()
    df[gene_col] = df[gene_col].astype(str).str.strip()
    genes = df[gene_col].replace("", pd.NA).dropna().drop_duplicates().tolist()
    return genes

def _genes_data_filtered(xlsx: Path, *, sheet: str, gene_col: str, score_col: str,
                         top_n: int = 10, low_thr: float = 0.5) -> List[str]:
    """
    数据清洗
    """
    df = pd.read_excel(
        xlsx, sheet_name=sheet, engine="openpyxl"
    )[[gene_col, score_col]].copy()

    df[gene_col] = df[gene_col].astype(str).str.strip()
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    df = df.dropna(subset=[gene_col, score_col]) \
           .drop_duplicates(subset=[gene_col], keep="first")

    # Score最高的前top_n个
    top = (
        df.sort_values(score_col, ascending=False)
          .head(top_n)[gene_col]
          .tolist()
    )

    low = (
        df[df[score_col] < low_thr]
        .sort_values(score_col, ascending=False)
        .head(top_n)[gene_col]
        .tolist()
    )

    # 合并去重
    seen = set()
    out: List[str] = []
    for g in top + low:
        if g not in seen:
            seen.add(g)
            out.append(g)

    return out
def render_all_keywords(root_dir: str, global_json: str, cfg: Cfg = Cfg()) -> List[str]:
    """
    可视化
    """
    root = Path(root_dir)
    interactions = json.loads(Path(global_json).read_text(encoding="utf-8"))
    if not isinstance(interactions, dict):
        raise ValueError("Global JSON must be a dict keyed by keyword.")

    out_paths: List[str] = []
    folders = sorted([p for p in root.iterdir() if p.is_dir() and p.name.endswith(cfg.dir_suffix)])

    for folder in folders:
        keyword = folder.name[:-len(cfg.dir_suffix)].strip()
        xlsx = folder / cfg.rank_excel
        if not xlsx.exists():
            continue

        kw_interactions = interactions.get(keyword, {})

        # ---------- rank: top 10 ----------
        genes_rank = _top_genes(xlsx, cfg)  # 你原来的函数：按 score 取前 cfg.top_genes
        payload_rank = build_payload(kw_interactions, genes_rank, cfg.top_partners)

        out_rank = folder / f"{keyword}_Interacting_Proteins_{cfg.top_partners}.html"
        _write_html(out_rank, keyword=keyword, payload_by_gene=payload_rank)
        out_paths.append(str(out_rank))
        if cfg.open_browser:
            webbrowser.open_new_tab("file://" + str(out_rank.resolve()))

        # ----------  data: filtered genes ----------
        try:
            genes_data = _genes_data_filtered(
                xlsx,
                sheet="data",
                gene_col=cfg.gene_col,
                score_col=cfg.score_col,
                top_n=20,
                low_thr=0.5
            )
        except Exception:
            continue

        payload_data = build_payload(kw_interactions, genes_data, cfg.top_partners)

        out_data = folder / f"{keyword}_Interacting_Proteins_{cfg.top_partners}_data.html"
        _write_html(out_data, keyword=keyword, payload_by_gene=payload_data)
        out_paths.append(str(out_data))
        if cfg.open_browser:
            webbrowser.open_new_tab("file://" + str(out_data.resolve()))

    return out_paths
