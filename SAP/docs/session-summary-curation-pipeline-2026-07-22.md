# 会话总结：SAP 分析与 Custom Module 定量管线

**日期：** 2026-07-22（更新 2026-07-26）  
**SAP 版本：** `seq-dep-SAP-v0.2.md`（Yanshan Jin）  
**项目：** 枸杞相关肠道微生物宏基因组/宏转录组测序深度 benchmarking

---

## 1. SAP 要做什么

### 1.1 主目标

在 **4 名受试者 × day0/day180 × metaG/metaT** 的高深度数据上，通过 **seqtk 降采样**（2M–50M read pairs × 5 seeds → 约 560 个数据集），确定：

> **能充分保留枸杞相关肠道微生物功能/分类信号的最低测序深度**，为正式临床试验选深度和组学类型。

### 1.2 主要分析维度（SAP §5.3 / §7–8）

- 预定义 feature 的 **检出率、丰度 concordance（vs 50M 参考）、seed 稳定性、模块完整度、纵向信号保留**
- Feature 类别包括：glycan、polyphenol、SCFA、carotenoid、goji 相关 taxa、全局功能、bile acid、多样性等
- 决策阈值（默认）：检出 ≥90%、Spearman ≥0.90、CV ≤30%、模块/纵向保留 ≥80%

### 1.3 手工 curated 特征集（`SAP/` 内 xlsx）

| 文件 | 内容 | 行数（约） |
|---|---|---|
| `gene_curation_Carotenoids.xlsx` | 类胡萝卜素 KO/EC | 97 |
| `gene_curation_LBP.xlsx` | 枸杞多糖降解 CAZyme（EC+GH 家族） | 87 |
| `gene_curation_LPS_biosynthesis_BRITE_ko01005.xlsx` | LPS 合成 KO | 162 |
| `gene_curation_Polyphenol.xlsx` | 多酚 → CAMPER/KO 映射 | 140 |
| `gene_curation_SCFA_48_enzyme.xlsx` | 48 酶系统（Frolova 2022） | 48 |
| `gene_curation_microbial_bile_acid.xlsx` | 微生物胆汁酸 KO | 21 |
| `goji_genus_level_changes_1.xlsx` | 枸杞干预属水平变化文献证据 | 21 |

SAP 正文见：`SAP/seq-dep-SAP-v0.2.md`（由 docx 转换，便于版本管理）。

---

## 2. 目前做到哪里

### 2.1 已完成的数据与主 pipeline

- **DRAM 注释 + distill + abundance** 主流程：`scripts/annotation_contigs_abundance.smk`
- 目录结构：
  - `annotation/dram_contigs_T2T/{omic}/{sample}/{prefix}/` — annotate、distill、`distill_input_fixed/`
  - `abundance/dram_contigs_T2T/{omic}/{sample}/{prefix}/` — featureCounts、`id_tpm`、DRAM 原生模块表
- **555 个 subsample** 已同时具备：
  - `distill_input_fixed/annotations.tsv`
  - `gene_counts.tsv`
- 每个 subsample 深度 × seed 命名如：`MG0001_50M_seed11`

### 2.2 Custom Module 管线（独立于主 smk）

| 文件 | 作用 |
|---|---|
| `scripts/build_curation_references.py` | SAP xlsx → 归一化 reference TSV（`id_only` / `id_ec` 两版） |
| `scripts/curation_refs/{id_only,id_ec}/*.tsv` | 6 个 curated 模块 reference |
| `scripts/curation_module_abundance.py` | 单 subsample：gene TPM → id_tpm → curated 模块丰度 |
| `scripts/aggregate_curation_long.py` | 跨 subsample 汇总长表（depth/seed/omic/id_variant） |
| `scripts/curation_modules.smk` | **独立** Snakemake，不改原 `annotation_contigs_abundance.smk` |
| `scripts/run_curation_modules.sbatch` | HPC 提交 id_only 全量 |
| `scripts/run_curation_modules_id_ec.sbatch` | HPC 提交 id_ec 全量 + 独立长表 |

**输出位置（每 subsample）：**

```
annotation/dram_contigs_T2T/{omic}/{sample}/{prefix}/custom_modules/
  id_only/
    {carotenoids,lbp,lps,polyphenol,bile_acid,scfa48}_per_feature.tsv
    {carotenoids,lbp,lps,polyphenol,bile_acid,scfa48}_per_module.tsv
  id_ec/
    （同上，sensitivity 版）
```

**跨样本长表（两版独立文件，互不覆盖）：**

```
annotation/dram_contigs_T2T/curation_module_summary/
  all_modules_per_module_long_id_only.tsv   # 主分析（69,375 行）
  all_modules_per_module_long_id_ec.tsv     # sensitivity（69,375 行）
```

### 2.3 验证与全量运行状态

| 运行 | Job ID | 状态 | 说明 |
|---|---|---|---|
| id_only 全量 | 19075537 | **COMPLETED**（~19 min） | 555 subsample × `custom_modules/id_only/` |
| id_ec 全量 | 19109573 | **COMPLETED**（~22 min） | 555 subsample × `custom_modules/id_ec/` |
| id_ec 首次提交 | 19099870 | FAILED | Snakemake 不能以 wildcard rule 名作 target；已修 sbatch |

- Smoke test：`MG0001_50M_seed11` 上 id_only / id_ec 均跑通并做过对比
- `~/2026-ysjin-seqdep` 为 Quobyte 软链接；日志可能显示 `/quobyte/angelazgrp/...` 路径，物理上为同一目录

### 2.4 刻意保留的「双版本」

| 特征 | DRAM 原生 | Custom curated |
|---|---|---|
| SCFA | `functional_abundance` + `scfa_per_*`（旧 reference） | `scfa48_per_*` |
| Polyphenol | CAMPER heatmap（`functional_abundance`） | `polyphenol_per_*`（见 §3.3 Polyphenol 策略） |
| Glycan/CAZy | DRAM CAZy heatmap | `lbp_per_*`（按具体枸杞多糖结构 group） |

---

## 3. 怎么做的（技术要点）

### 3.1 用哪个 annotations 文件？

**必须用：**

```
annotation/.../{prefix}/distill_input_fixed/annotations.tsv
```

**不要用顶层 `annotations.tsv`**。二者差异：

- `cazy_best_hit` 去掉 `.hmm` 后缀（`GH127.hmm` → `GH127`）— 对 LBP CAZy 匹配至关重要
- 纯数字 `fasta` ID 加 `bin.` 前缀（distill 兼容）

主 pipeline 与 custom pipeline 均读 fixed 版。

### 3.2 定量流程

```
featureCounts gene_counts.tsv
    + distill_input_fixed/annotations.tsv
        ↓
gene TPM（长度 + library size 归一化）
        ↓
id_tpm（每个 KO / EC / Pfam / CAZy / CAMPER D-id 等的 community 汇总 TPM）
        ↓
curation reference 每行 ids → feature_tpm = sum(id_tpm[id])
        ↓
按 module 列 group → detection_fraction、module_sum_tpm、module_weighted_tpm
```

- ID 提取复用 DRAM：`get_ids_from_annotations_by_row()`（与 distill 一致）
- 基因 TPM **均分**到其所有 ID；curation 一行内多 ID 为 **OR**，community 水平求和

### 3.3 ID 匹配策略：id_only vs id_ec

两套 reference 由 `build_curation_references.py` 从同一 xlsx 生成；组件列（`ko_ids` / `ec_ids` / `cazy_ids` / `camper_ids`）相同，仅拼出的 `ids` 不同：

| 变体 | `ids` 包含 | 不含 |
|---|---|---|
| **`id_only`**（**SAP 主分析**） | KO + CAZy + CAMPER `D#####` | EC |
| **`id_ec`**（sensitivity） | KO + EC + CAZy + CAMPER `D#####` | — |

### 3.4 各模块 reference ID 类型

| 模块 | id_only 的 `ids` | id_ec 的 `ids` | 备注 |
|---|---|---|---|
| carotenoids, LPS, bile_acid, scfa48 | KO | KO + EC | scfa48 多 KO 酶从 `KO_Mapping_Expanded` 合并 |
| **LBP (glycan)** | **CAZy** | CAZy + EC | 基本无 KO |
| **polyphenol** | **KO + CAMPER D-id** | **KO + EC + CAMPER D-id** | 见下节 |

#### Polyphenol 策略（已决定）

**Custom curated polyphenol** = SAP xlsx 手工映射 + 上表 ID 规则；**不是**单独依赖 DRAM CAMPER heatmap（后者仍保留作对照）。

- **Reference 来源列**（`gene_curation_Polyphenol.xlsx`）：`CAMPER raw gene_id`、`CAMPER identifier`、`KEGG Orthology identifier`、`EC number from raw gene_id`
- **Module 分组**：`Raw compound`（如 `Myricetin`、`rutin`），不是 CAMPER module 名
- **示例（id_only）**：FCR → `K26178,D00001`；CHI → `K26177,D00003`
- **示例（id_ec）**：CHI → `K26177,EC:5.5.1.6,D00003`
- **局限**：community `id_tpm` 仅 **~15 种 D-id**，多数 feature 仍主要靠 KO 命中；**~15/140** 行 `n_ids=0`（如 `mul-berin` unresolved）

**与 session 初稿选项的对应：**

| 初稿选项 | 是否采用 |
|---|---|
| 仅 KO/EC | 否（主分析含 D-id） |
| 加入 CAMPER D-id | **是**（id_only / id_ec 均含） |
| 主要依赖 DRAM CAMPER 版 | 否（并行保留，不作 custom 主表） |

### 3.5 Module 分组含义（`*_per_module.tsv`）

| feature_class | `module` = |
|---|---|
| LBP_glycan | 一种多糖名（如 `LBGP70-OL`） |
| SCFA48 | `产物:通路变体`（如 `butyrate:P1 Acetyl-CoA`）— **不同丁酸通路分开** |
| Polyphenol | 化合物（如 `Myricetin`、`rutin`） |
| Carotenoids | 类胡萝卜素衍生物 |
| Bile_acid | 功能模块（如 `7α-dehydroxylation`） |
| LPS | lipid A / core / O-antigen 等 |

**未做：** 多糖降解 ↔ SCFA 产生的跨模块代谢链链接。

### 3.6 SCFA48 多 KO 处理

- Expanded mapping 中同一酶的多 KO **合并为一行** `ids`（不拆多行 feature）
- 丰度：**OR 检出**（任一 ID > 0 即 detected）+ **TPM 求和**
- 与 id_ec 的差异：id_only 不含 EC，检出更严（如 50M smoke test：26/48 vs 39/48 detected）

---

## 4. 遇到的问题

### 4.1 环境与基础设施

- **Quobyte 挂载偶发卡死**：Glob/`ls`/`find` 长时间无响应；IDE 可见文件但工具读不到
- **fnm 权限错误**：`~/.local/state/fnm_multishells/` Permission denied，每次 shell 启动 ~2 分钟
- **Snakemake cache**：需 `export XDG_CACHE_HOME=$PWD/.cache` 避免 `~/.cache` 写权限问题
- **项目路径**：`~/2026-ysjin-seqdep` → `/quobyte/angelazgrp/2026-ysjin-seqdep/` 软链接；`Path.resolve()` 日志显示后者

### 4.2 设计与理解上的纠偏

| 问题 | 结论 |
|---|---|
| 之前以为 `id_tpm` 无 CAMPER id | **有误**：样本内有 **15 种 `D#####`**（如 D00001），但远少于 CAMPER 全目录 |
| `metabolism_summary` 似乎全是 KO | 展示列偏 KO/CAZy；**`id_tpm` 含 KO(3486)、EC(1221)、Pfam、CAZy、15×D-id** |
| KO+EC「同一基因」 | **确认方式**：`annotations.tsv` **一行 = 一个 gene**；`ko_id` + `kegg_hit` 中 `[EC:...]` 同行提取 |
| KO+EC 相加是否 double count | **基因水平**：多 ID 均分 TPM，KO+EC 相加 ≤ 基因 TPM；**community/id_tpm 水平**：curation 求和是 OR 匹配，不是单基因逻辑 |
| CAMPER sheet 为何 K 和 D 都有 | **CAMPER_distillate** 本身混用：303 行 KO 开头、49 行 D 开头、312 行组合写法 |
| 只有 `id_only/` 文件夹 | **正常**：id_only 与 id_ec 分目录；全量 id_ec 需单独 job（19099870 失败后 19109573 已成功） |
| 长表被覆盖风险 | 已拆为 `*_id_only.tsv` / `*_id_ec.tsv` 两个独立文件 |

### 4.3 SAP / 数据小问题（待修）

- SAP md 若干拼写笔误（`inlcude`、`annoation`、`constrairs` 等）
- `goji_genus_level_changes_1.xlsx` 最后一行 DOI 粘连文字
- Custom polyphenol ~15/140 行 `n_ids=0`（化合物 unresolved / 无 KO）

### 4.4 SCFA「替换」决策

用户曾选「用 48 酶版替换旧 scfa_reference」，但实施时 **保留双版本**：

- 旧 `scfa_reference.tsv` 覆盖 acetate/lactate/ethanol 等 + **essential-step 完整度逻辑**
- `scfa48` 按酶/通路变体分组，**无 essential 权重**
- 48 酶表实际也含 acetate/formate 等模块，将来可评估合并

---

## 5. 仍存在的困惑 / 待决事项

1. **Custom module 的 OR + sum 匹配是否满足 SAP 统计意图？** id_only 已比 id_ec 更严；是否还需 gene-level 或 essential-step 权重？
2. **Module 完整度**：当前为简单 `detection_fraction`（检出 feature 数 / 总数），**无 essential-step 权重**；SAP §7 模块 80% 阈值如何映射到各 curation 类型？
3. **Module 内 TPM 相加**可能在共享 KO / 多亚基 / id_ec 下偏高；是否需 dedup 或 max/min 规则？
4. **LBP 用 CAZy 家族** — 深度 benchmarking 反映的是 breadth 而非精确基因？
5. **与 DRAM 原生表的对应关系** — 正式分析时以哪套为主、如何并排报告？
6. **SCFA 分析用哪套** — `scfa48_per_module` vs DRAM `scfa_per_pathway` vs 两者对照

---

## 6. 打算怎么解决 / 可选改进

| 方向 | 做法 |
|---|---|
| SAP 深度分析 | 基于 `all_modules_per_module_long_id_only.tsv` 按 depth 算 median detection_fraction、相对 50M 保留率 |
| id_ec sensitivity | 与 id_only 长表并排比较（同路径结构，不同 `id_variant` 列） |
| Essential steps | 在 reference 加 `essential` 列，module 完整度改 weighted logic（仿 `scfa_reference`） |
| 多亚基 / 同源 KO | `curation_module_abundance.py` 加 `match_mode`（`or_sum` / `and_min` / `ko_only`） |
| 文档 | SAP Table 4 可从 `curation_refs/` 自动生成 |

---

## 7. 下一步建议（优先级）

1. **【高】SAP 深度 benchmarking 统计分析** — 读 `all_modules_per_module_long_id_only.tsv`，按 SAP §8 做 depth × feature_class 检出率 / concordance（Figure 3–6 方向）

2. **【中】id_only vs id_ec sensitivity 摘要** — 对 carotenoids、scfa48 等 EC 敏感模块做差异汇总

3. **【中】明确 SCFA 正式报告用哪套** — `scfa48_per_module` vs DRAM `scfa_per_pathway`

4. **【低】修 SAP 笔误、goji xlsx DOI、补充 SAP §6.3 一句说明 CAMPER 经 DRAM `--use_camper` + custom 多酚 curation 并行

---

## 8. 关键文件索引

```
SAP/
  seq-dep-SAP-v0.2.md
  gene_curation_*.xlsx
  goji_genus_level_changes_1.xlsx
  docs/
    session-summary-curation-pipeline-2026-07-22.md   ← 本文件

scripts/
  annotation_contigs_abundance.smk    # 主 pipeline（未改 custom 逻辑）
  curation_modules.smk                # 独立 curated module pipeline
  build_curation_references.py
  curation_module_abundance.py
  aggregate_curation_long.py
  curation_refs/
    id_only/*.tsv
    id_ec/*.tsv
  run_curation_modules.sbatch         # id_only 全量
  run_curation_modules_id_ec.sbatch   # id_ec 全量
  scfa_reference.tsv                  # DRAM 原生 SCFA 用

annotation/dram_contigs_T2T/
  {omic}/{sample}/{prefix}/
    distill_input_fixed/annotations.tsv   # 定量必读
    custom_modules/
      id_only/                            # 主分析 per-subsample 输出
      id_ec/                              # sensitivity per-subsample 输出
  curation_module_summary/
    all_modules_per_module_long_id_only.tsv
    all_modules_per_module_long_id_ec.tsv

abundance/dram_contigs_T2T/
  {omic}/{sample}/{prefix}/
    gene_counts.tsv
    id_tpm.tsv
    functional_abundance.tsv              # DRAM 原生（含 CAMPER）
    scfa_per_*.tsv
```

---

## 9. 一句话状态

**SAP 与 6 套 gene curation 已就绪；DRAM 主注释/丰度与 custom module（id_only + id_ec，各 555 subsample）均已全量跑完；跨深度长表 `all_modules_per_module_long_id_only.tsv` 可进入 SAP 规定的 depth benchmarking 统计分析阶段。**
