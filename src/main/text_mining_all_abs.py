"""
该模块定义函数，设计代码对一个关键词进行文本挖掘，挖掘范围为所有pum摘要

"""
from pathlib import Path
import json
import re
import os
import files_operation as fo
from collections import OrderedDict


# ————————————————————————————————————————————————————————————————————————————————————
# ————————————————————————————————————————————————————————————————————————————————————
# 一、（1）定义函数从头检索一个关键词:NLP文献挖掘
import spacy
from spacy.matcher import PhraseMatcher
from spacy.tokens import DocBin

# --------------------------------------------------------------
# 模块级缓存
_mesh_index_cache = {}
_cache_file_list_cache = {}


# --------------------------------------------------------------

def normalize_text(text):
    if text is None:
        return ""
    text = str(text).lower()
    text = re.sub(r"[-_/]+", " ", text)
    text = re.sub(r"[^\w\s]", " ", text)
    text = re.sub(r"\s+", " ", text).strip()
    return text


def unique_keep_order(items):
    return list(OrderedDict.fromkeys(x for x in items if x))

# --------------------------------------------------------------
# MeSH新格式词库读取与索引
def load_mesh_synonyms(file_paths, allowed_record_types=None):
    """
    读取新格式MeSH词库：
    record_type|||record_id|||preferred_name|||term1\tterm2\tterm3
    """
    global _mesh_index_cache

    if isinstance(file_paths, (str, Path)):
        file_paths = [file_paths]

    cache_key = (
        tuple(str(Path(p)) for p in file_paths),
        tuple(sorted(allowed_record_types)) if allowed_record_types is not None else None
    )
    if cache_key in _mesh_index_cache:
        return _mesh_index_cache[cache_key]

    record_rows = {}
    term_to_record_keys = {}

    for file_path in file_paths:
        path = Path(file_path)
        if not path.exists():
            continue

        for line in path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line or "|||" not in line:
                continue

            parts = line.split("|||", 3)
            if len(parts) != 4:
                continue

            record_type, record_id, preferred_name, terms_part = parts
            record_type = record_type.strip()
            record_id = record_id.strip()
            preferred_name = preferred_name.strip()

            if allowed_record_types is not None and record_type not in allowed_record_types:
                continue

            terms_raw = [x.strip() for x in terms_part.split("\t") if x.strip()]
            if preferred_name:
                terms_raw = [preferred_name] + terms_raw

            terms_raw = unique_keep_order(terms_raw)
            terms_norm = unique_keep_order([normalize_text(x) for x in terms_raw])

            record_key = f"{record_type}|||{record_id}"
            record_rows[record_key] = {
                "record_type": record_type,
                "record_id": record_id,
                "preferred_name": preferred_name,
                "terms_raw": terms_raw,
                "terms_norm": terms_norm,
            }

            for term in terms_norm:
                term_to_record_keys.setdefault(term, set()).add(record_key)

    _mesh_index_cache[cache_key] = (record_rows, term_to_record_keys)
    return _mesh_index_cache[cache_key]


def resolve_kw(kw, mesh_data_files=None, allowed_record_types=None):
    """
    根据输入关键词 kw，在新格式词库中定位唯一记录。
    找不到返回 None；若匹配到多个记录则报错。
    """
    if mesh_data_files is None:
        mesh_data_files = [
        "Mesh_Data/Mesh_Desc_Synonyms.txt",
        "Mesh_Data/Mesh_Supp_Synonyms.txt",

    ]
    if allowed_record_types is None:
        allowed_record_types = {"desc", "supp"}

    kw_norm = normalize_text(kw)
    record_rows, term_to_record_keys = load_mesh_synonyms(
        mesh_data_files,
        allowed_record_types=allowed_record_types
    )

    record_keys = term_to_record_keys.get(kw_norm, set())

    if len(record_keys) == 0:
        return None

    if len(record_keys) > 1:
        raise ValueError(
            f"关键词'{kw}'在词库中匹配到多个记录：{sorted(record_keys)}"
        )

    return record_rows[list(record_keys)[0]]


# --------------------------------------------------------------
# 缓存摘要读取
def get_cache_paths(cache_dir, file_stem):
    docbin_path = Path(cache_dir) / f"{file_stem}.spacy"
    pmid_path = Path(cache_dir) / f"{file_stem}_pmid.txt"
    return docbin_path, pmid_path


def get_available_cache_files(cache_dir=None, start_num=None, end_num=None, file_template=None):
    """
    扫描一次可用缓存文件列表，后续复用。
    """
    global _cache_file_list_cache

    if cache_dir is None:
        cache_dir = "PM_NLP0329"
    if start_num is None:
        start_num = 1
    if end_num is None:
        end_num = 1334
    if file_template is None:
        file_template = "pubmed26n{num1}.hao"

    cache_key = (str(cache_dir), int(start_num), int(end_num), str(file_template))
    if cache_key in _cache_file_list_cache:
        return _cache_file_list_cache[cache_key]

    cache_files = []
    for num in range(start_num, end_num + 1):
        num1 = str(num).rjust(4, "0")
        stem = Path(file_template.format(num1=num1)).stem
        docbin_path, pmid_path = get_cache_paths(cache_dir, stem)

        if docbin_path.exists() and pmid_path.exists():
            cache_files.append((docbin_path, pmid_path))

    _cache_file_list_cache[cache_key] = cache_files
    return cache_files


def load_cached_docs(docbin_path, pmid_path, nlp):
    with open(docbin_path, "rb") as f:
        doc_bin = DocBin().from_bytes(f.read())
    docs = list(doc_bin.get_docs(nlp.vocab))

    with open(pmid_path, "r", encoding="utf-8") as f:
        pmids = [line.strip() for line in f if line.strip()]

    if len(docs) != len(pmids):
        raise ValueError(f"缓存不一致：{docbin_path} 与 {pmid_path}")

    return pmids, docs


# --------------------------------------------------------------
# NLP 匹配核心
def build_matcher(nlp, terms_raw, terms_norm):
    """
    构建 spaCy PhraseMatcher。
    """
    matcher = PhraseMatcher(nlp.vocab, attr="LOWER")
    for raw_term, norm_term in zip(terms_raw, terms_norm):
        if norm_term:
            matcher.add(raw_term, [nlp.make_doc(norm_term)])
    return matcher

# --------------------------------------------------------------
# 文献挖掘
def mining_kw_id_nlp(kw, output_filename, output_dir):
    """
    文献挖掘：

    参数：
        kw              : 待检索关键词
        output_filename : PMID输出文件名
        output_dir      : 输出文件夹名

    返回：
        id_file_path, mined_list
    """

    fo.ensure_dir_exists(output_dir)
    mesh_data_files = [
        "Mesh_Data/Mesh_Desc_Synonyms.txt",
        "Mesh_Data/Mesh_Supp_Synonyms.txt",

    ]
    allowed_record_types = {"desc", "supp"}
    cache_dir = "PM_NLP0329"
    start_num = 1
    end_num = 1334
    file_template = "pubmed26n{num1}.hao"

    # 1. 定位关键词对应记录
    record_row = resolve_kw(kw, mesh_data_files, allowed_record_types)

    if record_row is None:
        record_type = "NA"
        record_id = "NA"
        preferred_name = kw
        terms_raw = [kw]
        terms_norm = [normalize_text(kw)]

    else:
        record_type = record_row["record_type"]
        record_id = record_row["record_id"]
        preferred_name = record_row["preferred_name"]
        terms_raw = record_row["terms_raw"]
        terms_norm = record_row["terms_norm"]
        print(
            f"关键词: {kw}，使用记录:"
            f"{record_type}|||{record_id}|||{preferred_name}，"
            f"同义词数量: {len(terms_raw)}"
        )

    # 2. 构建 NLP 匹配器
    nlp = spacy.blank("en")
    matcher = build_matcher(nlp, terms_raw, terms_norm)

    # 3. 遍历缓存摘要并匹配
    id_set = set()

    cache_files = get_available_cache_files(
        cache_dir=cache_dir,
        start_num=start_num,
        end_num=end_num,
        file_template=file_template
    )

    for docbin_path, pmid_path in cache_files:
        pmids, docs = load_cached_docs(docbin_path, pmid_path, nlp)

        for pmid, doc in zip(pmids, docs):
            if matcher(doc):
                id_set.add(str(pmid))

    # 4. 输出PMID列表
    id_file_path = os.path.join(output_dir, output_filename)

    if id_set:
        mined_list = sorted(map(int, id_set))
        show_list = mined_list[-1:-101:-1] if len(mined_list) > 100 else mined_list
        print(
            f"共得到{len(mined_list)}条PubMed_ID。"
        )
    else:
        mined_list = []
        print(
            f"未匹配到 PubMed_ID。"
        )

    fo.txt_output(id_file_path, mined_list)
    return id_file_path, mined_list


# 下一个函数
# ——————————————————————————————————————————————————————————————————————————
# ——————————————————————————————————————————————————————————————————————————
# 二、定义函数，从ID信息到gene ratio信息
_Symbol_To_ID_JSON = "Symbol_To_PubMed_ID2026.json"
_Symbol_To_ID_Hao= "Symbol_2_PMID_2026.hao"
_PMID_To_Symbol_JSON = "PubMed_ID_To_Symbol2026.json"

_pmid_to_gene_cache = None



def build_pmid_to_gene_index_from_json(
    gene_to_pmid_json=None,
    pmid_to_gene_json=None
):
    """
    根据gene-PMID的json构建PMID-gene反向索引。
    """
    if gene_to_pmid_json is None:
        gene_to_pmid_json = _Symbol_To_ID_JSON
    if pmid_to_gene_json is None:
        pmid_to_gene_json = _PMID_To_Symbol_JSON

    with open(gene_to_pmid_json, "r", encoding="UTF-8") as f:
        gene_to_pmid = json.load(f)

    pmid_to_gene = {}

    for gene, pmid_list in gene_to_pmid.items():
        for pmid in pmid_list:
            pmid = str(pmid).strip()
            if not pmid or pmid == "XXXXX":
                continue
            pmid_to_gene.setdefault(pmid, []).append(gene)

    for pmid, genes in pmid_to_gene.items():
        pmid_to_gene[pmid] = unique_keep_order(genes)

    with open(pmid_to_gene_json, "w", encoding="UTF-8") as f:
        json.dump(pmid_to_gene, f, ensure_ascii=False)

    return pmid_to_gene


def build_pmid_to_gene_index_from_ori(
    hao_file=None,
    pmid_to_gene_json=None
):
    """
    根据Symbol_2_PMID_2026.hao构建PMID-gene反向索引。
    """
    if hao_file is None:
        hao_file = _Symbol_To_ID_Hao
    if pmid_to_gene_json is None:
        pmid_to_gene_json = _PMID_To_Symbol_JSON

    pmid_to_gene = {}

    with open(hao_file, "r", encoding="UTF-8") as f:
        for line in f:
            line = line.strip()
            if not line or "|||" not in line:
                continue

            gene, pmid_part = line.split("|||", 1)
            gene = gene.strip()

            pmid_list = [x.strip() for x in pmid_part.split(",") if x.strip()]
            for pmid in pmid_list:
                if pmid == "XXXXX":
                    continue
                pmid_to_gene.setdefault(pmid, []).append(gene)

    for pmid, genes in pmid_to_gene.items():
        pmid_to_gene[pmid] = unique_keep_order(genes)

    with open(pmid_to_gene_json, "w", encoding="UTF-8") as f:
        json.dump(pmid_to_gene, f, ensure_ascii=False)
    return pmid_to_gene


def load_pmid_to_gene_index(
    gene_to_pmid_json=None,
    gene_to_pmid_hao=None,
    pmid_to_gene_json=None
):
    """
    优先加载已经存在的反向索引；
    若不存在，则构建。
    """
    global _pmid_to_gene_cache

    if _pmid_to_gene_cache is not None:
        return _pmid_to_gene_cache

    if gene_to_pmid_json is None:
        gene_to_pmid_json = _Symbol_To_ID_JSON
    if gene_to_pmid_hao is None:
        gene_to_pmid_hao = _Symbol_To_ID_Hao
    if pmid_to_gene_json is None:
        pmid_to_gene_json = _PMID_To_Symbol_JSON

    if os.path.exists(pmid_to_gene_json):
        with open(pmid_to_gene_json, "r", encoding="UTF-8") as f:
            _pmid_to_gene_cache = json.load(f)
    else:
        if os.path.exists(gene_to_pmid_hao):
            _pmid_to_gene_cache = build_pmid_to_gene_index_from_ori(
                hao_file=gene_to_pmid_hao,
                pmid_to_gene_json=pmid_to_gene_json
            )
        else:
            _pmid_to_gene_cache = build_pmid_to_gene_index_from_json(
                gene_to_pmid_json=gene_to_pmid_json,
                pmid_to_gene_json=pmid_to_gene_json
            )

    return _pmid_to_gene_cache


def id_to_gene(id_info_path, output_file="id_to_gene.txt", output_dir="output_info"):
    """
    PMID-gene映射函数
    id_info_path : PMID 列表文件路径
    """
    pmid_to_gene = load_pmid_to_gene_index(
        gene_to_pmid_json=_Symbol_To_ID_JSON,
        gene_to_pmid_hao=_Symbol_To_ID_Hao,
        pmid_to_gene_json=_PMID_To_Symbol_JSON
    )

    with open(id_info_path, "r", encoding="UTF-8") as f:
        query_pmids = [line.strip() for line in f if line.strip()]

    query_pmids = list(set(query_pmids))

    gene_set = set()
    for pmid in query_pmids:
        genes = pmid_to_gene.get(str(pmid), [])
        gene_set.update(genes)

    fo.ensure_dir_exists(output_dir)
    gene_file_path = os.path.join(output_dir, output_file)

    if gene_set:
        id_gene_list = sorted(gene_set)
        fo.txt_output(gene_file_path, id_gene_list)
    else:
        id_gene_list = []
        fo.txt_output(gene_file_path, id_gene_list)

    return id_gene_list, gene_file_path

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# 下一个函数
# 定义函数，根据id对应的gene信息，输出ratio信息
def gene_to_ratio(kw, gene_info_path, id_info_path, output_file="id_to_gene.txt", output_dir="output_info"):
    """
    gene_info_path:kw对应的gene报道信息的文件路径
    id_info_path:kw对应的报道id信息的文件路径
    output_file:输出的文件
    output_dir:输出的文件夹
    """

    # 读取含有所有gene与id的json文件
    gene_id_info_d = fo.read_json_file("Symbol_To_PubMed_ID2026.json")# 26年摘要
    # 读取含有所有gene报道总量的txt文件
    gene_re_num_info_list = fo.read_txt_file('Symbol_PMID_Count_2026.txt')  # 26年摘要
    gene_total_num_info_d = {}
    for line_info in gene_re_num_info_list:
        gene_name = line_info.split("   ")[0]
        gene_total_value = line_info.split('   ')[-1]
        gene_total_num_info_d[gene_name] = gene_total_value

    # 对应gene信息，计算ratio
    fo.ensure_dir_exists(output_dir)
    gene_file_path = os.path.join(output_dir, output_file)
    if os.path.getsize(gene_info_path) > 0:
        gene_info_l = fo.read_txt_file(gene_info_path)
        # 读取报道的id文件为列表
        id_info_l = fo.read_txt_file(id_info_path)
        # 各gene的原始信息
        id_info_set = set(id_info_l)
        match_gene_mess_d = {}
        for match_gene in gene_info_l:
            match_gene_ids_all = gene_id_info_d[match_gene]
            matching_ids = set.intersection(set(match_gene_ids_all), id_info_set)
            match_gene_re_num = len(matching_ids)
            match_gene_re_total_num = int(gene_total_num_info_d[match_gene])
            match_gene_ratio = int(match_gene_re_num) / match_gene_re_total_num
            match_gene_ratio_info_l = [match_gene_re_num, match_gene_re_total_num, match_gene_ratio]
            match_gene_mess_d[match_gene] = match_gene_ratio_info_l
        # 新建字典-输出信息
        match_gene_output = {}
        if match_gene_mess_d:
            match_gene_mess_d_s = sorted(match_gene_mess_d.items(), key=lambda item: item[0], reverse=True)
            for (k, v) in match_gene_mess_d_s:
                match_gene_output[k] = v
            fo.txt_output(gene_file_path, match_gene_output)

        else:
            fo.emtpy_file_output(gene_file_path)

    else:
        match_gene_output = {}
        fo.emtpy_file_output(gene_file_path)

    return match_gene_output, gene_file_path
