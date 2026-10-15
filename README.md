# METER NLP Procedure Guide

## Repository Structure

METER/
├── README.md
├── src/
│ ├── main/ # Main pipeline scripts
│ │ ├── 3_kws_mining_nlp.py # Main entry for keyword mining
│ │ ├── text_mining_all_abs.py # Core text mining logic
│ │ └── files_operation.py # File I/O utilities
│ └── analysis/ # Downstream analysis scripts
│ ├── 1_drug_hot_terms_match.py # Drug–hot term matching
│ ├── ppi1_match.py # PPI network matching
│ ├── ppi2_echarts.py # PPI network visualization
│ ├── heatmap1_data_to_csv.py # Heatmap data preparation
│ ├── heatmap2_groups.py # Heatmap grouping
│ ├── venn_groups.py # Venn diagram grouping
│ └── score_functions.py # Scoring functions
├── data/ # Small input and intermediate files
│ ├── Cela_DE.txt
│ ├── cela_keyword1002.txt
│ ├── cela_kw_250828_2.txt
│ ├── cela_kw_ros.txt
│ ├── EG_To_ENSP.txt
│ ├── ITCM_kws_250507.txt
│ ├── PPI_String_Human_Link_M400_New.txt
│ ├── PubMed_ID_To_Symbol2026.json
│ ├── Salvianolic acid B_DEGs.txt
│ ├── Symbol_2_PMID_2026.hao
│ ├── Symbol_PMID_Count_2026.txt
│ ├── Symbol_To_EG.txt
│ └── Symbol_To_PubMed_ID2026.json
└── METER_large_files/ # Large files (not included in this repo)
  ├── lexical library/
  ├── ABS.zip
  ├── Mesh_Data.zip
  ├── METER_V2_ECharts_Generate.exe
  └── PM_NLP0329.zip

### Data

Large files are not included in this repository due to GitHub size limits.
Please download them from Baidu Netdisk and place them in the `METER_large_files/` folder:
Files to download:

- Link: https://pan.baidu.com/s/1uEnPLWrEQQxenMnHfaSklw?pwd=7uh6 
- Extraction code: 7uh6

Files to download:

- `PM_NLP0329.zip`
- `Mesh_Data.zip`
- `ABS.zip`
- `METER_V2_ECharts_Generate.exe`
- `lexical library/` — cached mining results (keywords, gene IDs, DE ratios). 
  It will be automatically generated when running the code, but you can also 
  download the pre-built cache from Baidu Netdisk to save time.

## STEP 1: Prepare METER Folder

Go to the METER folder (`METER_main_procedure_NLP`) and unzip:
- `PM_NLP0329.zip`
- `Mesh_Data.zip`
- `ABS.zip`
- `METER_V2_ECharts_Generate.exe`
- `lexical library/` — cached mining results (keywords, gene IDs, DE ratios). 
  It will be automatically generated when running the code, but you can also 
  download the pre-built cache from Baidu Netdisk to save time.

## STEP 2: Install Dependencies

Install the required Python packages:

```bash
python -m pip install matplotlib pandas numpy seaborn venn spacy xlsxwriter openpyxl requests lxml beautifulsoup4 networkx pyecharts Pillow
```

## STEP 3: Input Files

Place the following files in the METER folder or another specified path:

### a. `keyword.txt`

- One keyword per line.
- Multiple keywords allowed.

### b. `genes_id.txt`

- One gene symbol per line.
- No header.
- No empty lines.

## STEP 4: Run the Program

Open a terminal in the METER folder and run:

```
python src/main/3_kws_mining_nlp.py
```

When prompted, enter:

```
the keyword file path(eg:keyword.txt): keyword.txt

the DEGs file path(eg:genes_id.txt): genes_id.txt
```

### If a Keyword Is Not in the Local Library

You will be asked to select an option:

```
Press '1' means the system will perform text mining automatically.

Press '2' means you will upload your own Pubmed ID associated with the keyword.

Press '3' means the system will pass the keyword.
```

### Note

- Option 1 mines the keyword from the local PubMed abstract cache (`PM_NLP0329`).
- Mining results are cached in `lexical library/` and reused next time.

## STEP 5: Output

The output folder is created automatically:

- Single keyword:

```
METER_<keyword>_DE
```

- Multiple keywords:

```
METER_<N>kws_DE
```
