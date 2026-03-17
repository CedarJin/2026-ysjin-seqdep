# Min Contig Length Analysis Report

## Dataset

Sample: `SRR10692699_50M_seed11`

This report combines two analyses:

1. Read alignment rate and assembly-size trends across contig length thresholds.
2. DRAM annotation behavior across contig length thresholds.

The alignment scan includes thresholds `200`, `500`, `1000`, and `2500` bp.
The DRAM annotation comparison also uses `200`, `500`, `1000`, and `2500` bp.

## Part 1: Alignment Scan

Source table: [alignment_scan.tsv](/group/datalabgrp/ctbrown/ysjin/2026-ysjin-seqdep/analysis/contig_len_alignment_scan/SRR10692699_50M_seed11/alignment_scan.tsv)

### Overall alignment rate

Figure: [alignment_scan.svg](/group/datalabgrp/ctbrown/ysjin/2026-ysjin-seqdep/analysis/contig_len_alignment_scan/SRR10692699_50M_seed11/alignment_scan.svg)

Interpretation:

- The overall alignment rate declines gradually as the minimum contig length increases: `98.36%` at `200 bp`, `97.89%` at `500 bp`, `96.68%` at `1000 bp`, and `93.04%` at `2500 bp`.
- This suggests that shorter contigs do contribute mappable sequence content, but the reduction is modest between `200` and `1000 bp` and becomes more pronounced at `2500 bp`.
- In practice, the alignment-rate penalty for moving from `200` to `500 bp` is small, while the penalty from `1000` to `2500 bp` is much larger.

### Contig count

Figure: [contig_count_scan.svg](/group/datalabgrp/ctbrown/ysjin/2026-ysjin-seqdep/analysis/contig_len_alignment_scan/SRR10692699_50M_seed11/contig_count_scan.svg)

Interpretation:

- Contig count drops sharply with stricter filtering: from `131,672` at `200 bp` to `13,617` at `2500 bp`.
- Most of the assembly fragmentation is concentrated in the short-contig range.
- This means that increasing the threshold removes a very large number of short contigs even when total retained sequence decreases much more slowly.

### Total assembly bp

Figure: [total_bp_scan.svg](/group/datalabgrp/ctbrown/ysjin/2026-ysjin-seqdep/analysis/contig_len_alignment_scan/SRR10692699_50M_seed11/total_bp_scan.svg)

Interpretation:

- Total retained assembly length decreases from `278.8 Mbp` at `200 bp` to `196.8 Mbp` at `2500 bp`.
- The assembly loses many contigs quickly, but not as much total sequence at first, indicating that many discarded contigs are individually short.
- The transition from `200` to `500 bp` removes a lot of contigs while preserving most total bp, which is often a favorable tradeoff.

### Average contig length

Figure: [avg_bp_scan.svg](/group/datalabgrp/ctbrown/ysjin/2026-ysjin-seqdep/analysis/contig_len_alignment_scan/SRR10692699_50M_seed11/avg_bp_scan.svg)

Interpretation:

- Average contig length increases monotonically with threshold, from about `2.1 kbp` at `200 bp` to `14.5 kbp` at `2500 bp`.
- This is expected, but it also highlights how strongly the retained assembly becomes enriched for longer contigs as the threshold rises.
- The jump is especially large above `1000 bp`, where the retained assembly begins to look substantially more contiguous.

### Maximum contig length

Figure: [max_bp_scan.svg](/group/datalabgrp/ctbrown/ysjin/2026-ysjin-seqdep/analysis/contig_len_alignment_scan/SRR10692699_50M_seed11/max_bp_scan.svg)

Interpretation:

- Maximum contig length stays constant at `907,670 bp` across all thresholds.
- This means the threshold only affects the lower tail of the contig distribution and does not change the longest assembled sequences.
- The main tradeoff is therefore between retaining small contigs and simplifying the assembly, not between retaining and losing the longest contigs.

### N50

Figure: [n50_bp_scan.svg](/group/datalabgrp/ctbrown/ysjin/2026-ysjin-seqdep/analysis/contig_len_alignment_scan/SRR10692699_50M_seed11/n50_bp_scan.svg)

Interpretation:

- N50 increases from `15,572 bp` at `200 bp` to `41,309 bp` at `2500 bp`.
- This reflects a strong improvement in the contiguity of the retained assembly subset as short contigs are filtered out.
- However, the N50 gain needs to be considered together with the alignment-rate loss and the annotation/runtime tradeoffs below.

### Alignment scan summary

- `200 bp` maximizes retained sequence and read alignment.
- `500 bp` appears to be a moderate filter that removes many short contigs with only a small alignment-rate penalty.
- `1000 bp` still preserves most alignment signal while substantially simplifying the assembly.
- `2500 bp` produces a much cleaner and more contiguous retained assembly, but at a noticeably larger cost in read alignment and total retained sequence.

## Part 2: DRAM Annotation Scan

Source table: [dram_minlen_summary.tsv](/group/datalabgrp/ctbrown/ysjin/2026-ysjin-seqdep/analysis/dram_contig_minlen_summary/SRR10692699_50M_seed11/dram_minlen_summary.tsv)

Note on metric definitions:

- `raw CAZy hits` are counted from non-empty `cazy_hits` rows in `annotations.tsv`.
- `raw CAMPER hits` are counted from non-empty `camper_hits` rows in `annotations.tsv`.
- `unique CAZy IDs` are counted from `cazy_ids`.
- `unique CAMPER IDs` are counted from `camper_id`.
- `functional module count` is defined here as the number of non-zero rows across functional sheets in `metabolism_summary.xlsx`, excluding `MISC`, `rRNA`, and `tRNA`.

### Raw hit counts

Figure: [dram_minlen_raw_hits.svg](/group/datalabgrp/ctbrown/ysjin/2026-ysjin-seqdep/analysis/dram_contig_minlen_summary/SRR10692699_50M_seed11/dram_minlen_raw_hits.svg)

Interpretation:

- Raw hit counts decline steadily with stricter thresholds.
- CAZy hits decrease from `7,486` at `200 bp` to `5,787` at `2500 bp`.
- CAMPER hits decrease from `1,077` at `200 bp` to `794` at `2500 bp`.
- This shows that shorter contigs do contribute additional annotated genes, especially when the threshold is reduced from `2500` to `1000` and below.

### Unique ID counts

Figure: [dram_minlen_unique_ids.svg](/group/datalabgrp/ctbrown/ysjin/2026-ysjin-seqdep/analysis/dram_contig_minlen_summary/SRR10692699_50M_seed11/dram_minlen_unique_ids.svg)

Interpretation:

- Unique ID counts are much more stable than raw hit counts.
- Unique CAZy IDs change only modestly: `481` at `200 bp` versus `471` at `2500 bp`.
- Unique CAMPER IDs also change, but less dramatically than raw hits: `93` at `200 bp` versus `80` at `2500 bp`.
- This suggests that lowering the threshold mostly adds more occurrences of already-observed functional categories, rather than introducing a large number of new functions.

### Runtime

Figure: [dram_minlen_runtime.svg](/group/datalabgrp/ctbrown/ysjin/2026-ysjin-seqdep/analysis/dram_contig_minlen_summary/SRR10692699_50M_seed11/dram_minlen_runtime.svg)

Interpretation:

- Runtime is highly sensitive to minimum contig length, but not perfectly monotonic because different runs used different cluster conditions and restart histories.
- The longest runtime is at `200 bp` (`59.78 h`), followed by `500 bp` (`30.52 h`).
- The `1000 bp` run is substantially shorter (`7.34 h`), while the `2500 bp` default run is `11.52 h`.
- Even allowing for scheduler effects, the dominant pattern is clear: including more short contigs greatly increases annotation cost.

### Functional module count

Figure: [dram_minlen_functional_modules.svg](/group/datalabgrp/ctbrown/ysjin/2026-ysjin-seqdep/analysis/dram_contig_minlen_summary/SRR10692699_50M_seed11/dram_minlen_functional_modules.svg)

Interpretation:

- Functional module count decreases from `1,629` at `200 bp` to `1,509` at `2500 bp`.
- The decline is real, but smaller than the drop in raw hit counts.
- This again suggests that stricter filtering removes many annotated genes while preserving much of the higher-level functional profile.
- In other words, module-level functional coverage is relatively robust to moderate threshold increases.

### DRAM scan summary

- Lower thresholds increase annotation yield, especially raw CAZy and CAMPER hits.
- Unique functional diversity changes more slowly than raw hit abundance.
- Module-level functional content is relatively stable compared with gene-level hit counts.
- The main cost of lowering the threshold is runtime, which rises sharply at `200 bp` and `500 bp`.

## Integrated interpretation

Taken together, these plots suggest that the main tradeoff is not simply “more information versus less information,” but rather “marginal annotation gain versus very large computational cost.”

Key points:

- `200 bp` maximizes read alignment, retained sequence, raw annotation hits, and module count, but it is by far the most expensive annotation setting.
- `500 bp` retains much of the biological signal seen at `200 bp`, but runtime is still very high.
- `1000 bp` appears to be a strong compromise: it retains most alignment signal and much of the functional diversity while reducing runtime substantially.
- `2500 bp` is computationally simpler and yields a more contiguous retained assembly, but it discards a noticeable amount of read-supported and annotation-supported content.

## Suggested interpretation for threshold choice

If the goal is to maximize annotation recovery, `200 bp` is the most inclusive setting.
If the goal is to balance annotation yield with computational efficiency, `1000 bp` looks like the most defensible compromise in this dataset.
If the goal is to prioritize assembly contiguity and a conservative contig set, `2500 bp` is reasonable, but it comes with a clear loss of read-supported and annotation-supported content.
