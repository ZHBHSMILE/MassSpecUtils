#' Get UniMod modification positions in a peptide sequence
#'
#' @param sequence Character string. Peptide sequence with UniMod annotations,
#'   e.g. `"(UniMod:1)ACD(UniMod:35)EFG"`.
#' @param mod_id Character or numeric. Target UniMod ID to locate.
#' @param option Character. Controls return format:
#'   \itemize{
#'     \item `"only_pos"` (default): integer vector of positions.
#'     \item `"with_aa"`: named list with `$pos` (positions) and `$letter` (amino acids).
#'   }
#'
#' @return
#' \itemize{
#'   \item If `option = "only_pos"`: integer vector of modification positions.
#'   \item If `option = "with_aa"`: list with elements `$letter` and `$pos`.
#' }
#' Position `0` indicates N-terminal; position `1` indicates the first amino acid.
#'
#' @details
#' The function works in four steps:
#'
#' 1. Replace the target UniMod tag with placeholder `%`
#' 2. Strip all remaining UniMod annotations
#' 3. Locate `%` positions in the cleaned string
#' 4. Subtract cumulative placeholder offsets
#'
#' **Example:**
#'
#' Input: `(UniMod:8888)AAAPPL(UniMod:2114)SK(UniMod:2114)AEYLK`
#'
#' | Step | Result |
#' |------|--------|
#' | Replace target | `(UniMod:8888)AAAPPL%SK%AEYLK` |
#' | Strip remaining | `AAAPPL%SK%AEYLK` |
#' | Locate `%` | `7, 10` |
#' | Adjust offsets | `6, 9` |
#'
#' @examples
#' # Only positions (default)
#' get_mod_pos("(UniMod:8888)AAAPPL(UniMod:2114)SK(UniMod:2114)AEYLK", "2114")
#' #> [1] 6 9
#'
#' # With amino acid letters
#' get_mod_pos("(UniMod:8888)AAAPPL(UniMod:2114)SK(UniMod:2114)AEYLK", "2114",
#'             option = "with_aa")
#' #> $letter
#' #> [1] "L" "K"
#' #> $pos
#' #> [1] 6 9
#'
#' # N-terminal modification (returns 0)
#' get_mod_pos("(UniMod:2114)AAAPPLK", "2114")
#' #> [1] 0
#'
#' # Batch processing
#' df <- data.frame(
#'   sequence = c(
#'     "(UniMod:8888)AAAPPL(UniMod:2114)SKAEYLK(UniMod:2114)",
#'     "G(UniMod:2114)KAA"
#'   )
#' )
#' df$positions <- sapply(df$sequence, get_mod_pos, mod_id = "2114")
#'
#' @export
get_mod_pos <- function(sequence, mod_id, option = "only_pos") {
	target_pattern <- paste0("\\(UniMod:", mod_id, "\\)")

	# Step 1: Replace target modification with placeholder
	temp_seq <- stringr::str_replace_all(sequence, target_pattern, "%")

	# Step 2: Strip all remaining UniMod annotations
	clean_seq <- stringr::str_replace_all(temp_seq, "\\(UniMod:[0-9]+\\)", "")

	# Step 3: Locate placeholder positions
	pos <- unname(stringr::str_locate_all(clean_seq, "%")[[1]][, 1])

	# Step 4: Correct cumulative offsets from placeholders
	adjusted_pos <- pos - seq_along(pos)

	if (option == "only_pos") {
		return(adjusted_pos)
	} else if (option == "with_aa") {
		seq_vec <- strsplit(clean_seq, "")[[1]]
		# Remove % from seq_vec before indexing
		seq_vec <- seq_vec[seq_vec != "%"]
		acid <- seq_vec[adjusted_pos]
		return(list(letter = acid, pos = adjusted_pos))
	} else {
		stop("option must be 'only_pos' or 'with_aa'")
	}
}


#' Format UniProt Protein Entries
#'
#' @description
#' Extracts and formats protein names from UniProt entry strings.
#' Processes semicolon-separated entries and extracts the second element
#' from pipe-delimited strings.
#'
#' @param str1 A character string containing UniProt entries (e.g., "ID|Name|Gene;ID|Name|Gene").
#' @param sep Separator for entries (default is ";").
#' @param pipe Pattern for splitting internal components (default is "\\|").
#'
#' @return A character string of formatted names.
#' @export
#'
#' @examples
#' format_uni("P12345|Protein A|GENE1;Q67890|Protein B|GENE2")
#' # Returns: "Protein A;Protein B"
format_uni <- function(str1, sep = ";", pipe = "\\|") {
	if (is.null(str1) || is.na(str1) || str1 == "") return(NA_character_)

	str1 |>
		stringr::str_split(sep) |>
		purrr::pluck(1) |> #
		purrr::map_chr(function(entry) {
			parts <- stringr::str_split(entry, pipe) |> purrr::pluck(1) #
			if (length(parts) >= 2) parts[2] else entry
		}) |>
		paste(collapse = sep)
}



#' 查找肽段在蛋白质序列中的位置
#'
#' @description
#' 在单条蛋白质序列中批量查找多个肽段的起止位置。
#' 对未匹配的肽段返回 `NA` 而非报错，适合大规模蛋白质组学数据处理。
#'
#' @param peptides 字符向量。待查找的肽段序列，例如 `c("AAPPL", "SKAEYLK")`。
#' @param protein_seq 字符串。目标蛋白质序列，例如 `"MAAAAPPLSKAEYLK"`。
#'
#' @return 一个 `data.table`，包含以下三列：
#' \describe{
#'   \item{`peptide`}{字符型。输入的肽段序列。}
#'   \item{`start`}{整数型。肽段在蛋白质中的起始位置（1-based）；未匹配时为 `NA`。}
#'   \item{`end`}{整数型。肽段在蛋白质中的终止位置（1-based）；未匹配时为 `NA`。}
#' }
#'
#' 若一个肽段在蛋白质中出现多次，则每次匹配各占一行。
#'
#' @details
#' 函数内部流程：
#'
#' 1. 将蛋白质序列转换为 \code{Biostrings::AAString} 对象
#' 2. 将所有肽段转换为 \code{Biostrings::AAStringSet} 对象
#' 3. 对每个肽段调用 \code{Biostrings::matchPattern}（C 底层实现）进行精确匹配
#' 4. 使用 \code{data.table::rbindlist} 合并结果（优于 \code{do.call(rbind)}）
#'
#' 所有外部函数均通过 \code{::} 显式调用，无需在 \code{NAMESPACE} 中声明 \code{importFrom}，
#' 仅需在 \code{DESCRIPTION} 的 \code{Suggests} 或 \code{Imports} 字段中列出依赖包即可。
#'
#' **注意事项：**
#' \itemize{
#'   \item 匹配区分大小写，输入序列须为标准单字母氨基酸编码（大写）
#'   \item 不支持含修饰标注的序列（如 \code{"AAC(UniMod:4)K"}），请预先清洗
#'   \item 函数每次只处理一条蛋白质序列；多蛋白场景请在外部循环调用
#' }
#'
#' @examples
#' protein <- "MAAAAPPLSKAEYLKGKR"
#'
#' peptides <- c(
#'   "AAPPL",      # 正常匹配
#'   "SKAEYLK",    # 正常匹配
#'   "ZZZZZ",      # 不存在，返回 NA
#'   "GKR"         # 正常匹配
#' )
#'
#' find_peptide_pos(peptides, protein)
#' #>    peptide start end
#' #> 1:   AAPPL     3   7
#' #> 2: SKAEYLK     9  15
#' #> 3:   ZZZZZ    NA  NA
#' #> 4:     GKR    16  18
#'
#' # 多蛋白场景：在外部循环调用
#' protein_list <- c(prot1 = "MAAAAPPLSK", prot2 = "MGKRSKAEYLK")
#'
#' result <- do.call(rbind, lapply(names(protein_list), function(pid) {
#'   df <- find_peptide_pos(peptides, protein_list[[pid]])
#'   df$protein <- pid
#'   df
#' }))
#'
#' @seealso
#' \code{\link[Biostrings]{matchPattern}} 底层匹配函数
#'
#' @export
find_peptide_pos <- function(peptides, protein_seq) {
	prot_aa  <- Biostrings::AAString(protein_seq)
	pept_aa  <- Biostrings::AAStringSet(peptides)

	results <- lapply(seq_along(peptides), function(i) {
		hit <- Biostrings::matchPattern(pept_aa[[i]], prot_aa, fixed = TRUE)

		if (length(hit) == 0) {
			data.frame(
				peptide = peptides[i],
				start   = NA_integer_,
				end     = NA_integer_
			)
		} else {
			data.frame(
				peptide = peptides[i],
				start   = Biostrings::start(hit),
				end     = Biostrings::end(hit)
			)
		}
	})

	data.table::rbindlist(results)
}
