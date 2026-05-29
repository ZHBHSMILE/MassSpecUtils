
# MassSpecUtils

用于质谱蛋白质组学数据处理的 R 工具包，提供肽段位置查找、UniMod 修饰位点解析、UniProt 条目格式化等功能。

## 安装

```r
# 从 GitHub 安装
devtools::install_github("yourusername/MassSpecUtils")
```

## 依赖包

```r
install.packages(c("stringr", "purrr", "data.table"))

# Bioconductor
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("Biostrings")
```

## 函数列表

### `find_peptide_pos()` — 查找肽段在蛋白质中的位置

在单条蛋白质序列中批量查找多个肽段的起止位置，未匹配返回 `NA`，支持一个肽段多次出现。

```r
library(MassSpecUtils)

peptides <- c("AAPPL", "SKAEYLK", "ZZZZZ")
protein  <- "MAAAAPPLSKAEYLKGKR"

find_peptide_pos(peptides, protein)
#>    peptide start end
#> 1:   AAPPL     3   7
#> 2: SKAEYLK     9  15
#> 3:   ZZZZZ    NA  NA
#> 4:     GKR    16  18
```

大规模数据（如 10000 肽段 × 2500 蛋白）推荐并行处理：

```r
library(furrr)
library(data.table)

# 提前转换蛋白序列
prot_seqs        <- as.character(fa)
names(prot_seqs) <- names(fa)

plan(multisession, workers = parallel::detectCores() - 1)

results <- furrr::future_map(unique(df$Protein.Group), function(prot_id) {
  prot_seq     <- prot_seqs[[prot_id]]
  peptides_sub <- unique(df$Stripped.Sequence[df$Protein.Group == prot_id])
  pos_dt       <- MassSpecUtils::find_peptide_pos(peptides_sub, prot_seq)
  pos_dt$Protein.Group <- prot_id
  pos_dt
}, .options = furrr_options(seed = NULL))

plan(sequential)

pos_all <- data.table::rbindlist(results, fill = TRUE)

# join 回原始数据，自动处理一对多
df5 <- df |>
  dplyr::left_join(pos_all, by = c("Protein.Group", "Stripped.Sequence" = "peptide"))
```

---

### `get_mod_pos()` — 解析 UniMod 修饰位点

从含 UniMod 标注的肽段序列中提取指定修饰的氨基酸位置。位置 `0` 表示 N 端修饰，位置 `1` 表示第一个氨基酸。

```r
seq <- "(UniMod:8888)AAAPPL(UniMod:2114)SK(UniMod:2114)AEYLK"

# 只返回位置（默认）
get_mod_pos(seq, "2114")
#> [1] 6 9

# 同时返回氨基酸字母
get_mod_pos(seq, "2114", option = "with_aa")
#> $letter
#> [1] "L" "K"
#> $pos
#> [1] 6 9

# N 端修饰返回 0
get_mod_pos("(UniMod:2114)AAAPPLK", "2114")
#> [1] 0
```

批量处理：

```r
df$positions <- sapply(df$sequence, get_mod_pos, mod_id = "2114")
```

---

### `format_uni()` — 格式化 UniProt 条目

从 UniProt 标准格式字符串（`ID|Name|Gene`）中提取蛋白质名称，支持多条目分号分隔输入。

```r
format_uni("P12345|Protein A|GENE1;Q67890|Protein B|GENE2")
#> [1] "Protein A;Protein B"
```

## 典型工作流

```r
library(MassSpecUtils)
library(dplyr)

# 1. 格式化蛋白名称
df <- df |>
  mutate(Protein.Name = sapply(Protein.Group, format_uni))

# 2. 解析修饰位点
df <- df |>
  mutate(mod_pos = sapply(Modified.Sequence, get_mod_pos, mod_id = "2114"))

# 3. 查找肽段位置
prot_seqs        <- as.character(fa)
names(prot_seqs) <- names(fa)

pos_all <- ... # 见上方并行示例

df <- df |>
  left_join(pos_all, by = c("Protein.Group", "Stripped.Sequence" = "peptide"))
```

## 开源协议

MIT

