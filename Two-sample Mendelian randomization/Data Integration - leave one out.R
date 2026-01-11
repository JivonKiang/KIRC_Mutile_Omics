rm(list = ls())
library(readxl)

# 获取当前目录的一级子文件夹路径
sub_dirs <- list.dirs(getwd(), recursive = FALSE)

# 初始化有效文件列表
valid_files <- character(0)

# 遍历每个子文件夹
for (dir in sub_dirs) {
  # 获取该子文件夹下的所有.xlsx文件
  dir_files <- list.files(
    path = dir,
    pattern = "\\.xlsx$",
    full.names = TRUE,
    recursive = FALSE
  )
  
  if (length(dir_files) == 0) next  # 跳过空文件夹
  
  # 检查文件夹中是否存在至少一个有效文件（含目标工作表）
  has_valid_file <- any(sapply(dir_files, function(file) {
    "Leave-One-Out" %in% excel_sheets(file)
  }))
  
  if (has_valid_file) {
    # 筛选该文件夹中包含目标工作表的文件
    valid_dir_files <- sapply(dir_files, function(file) {
      if ("Leave-One-Out" %in% excel_sheets(file)) file else NA
    })
    valid_dir_files <- na.omit(valid_dir_files)
    valid_files <- c(valid_files, valid_dir_files)
  }
}

# 读取所有有效文件的数据
data_list <- lapply(valid_files, function(file) {
  read_excel(file, sheet = "Leave-One-Out")###      Leave-One-Out
})

# 打印结果
print(data_list)

# 可选：为列表元素命名（方便后续操作）
names(data_list) <- basename(valid_files)

str(data_list)

# 提取所有MR_cell_Results.xlsx的数据
cell_data <- data_list[grepl("MR_cell_Results\\.xlsx", names(data_list))]

# 合并为单个数据框
combined_cell <- do.call(rbind, cell_data)

# 查看合并后的结构
str(combined_cell)

# 提取所有MR_gene_Results.xlsx的数据
gene_data <- data_list[grepl("MR_gene_Results\\.xlsx", names(data_list))]

# 合并为单个数据框
combined_gene <- do.call(rbind, gene_data)

# 查看合并后的结构
str(combined_gene)

library(openxlsx)
phenotype <- read.xlsx("R10_manifest_芬兰  表型下载地址.xlsx",sheet = "Sheet2")

phenotype

# 假设你的数据框名为df
phenotype$phenocode <- paste0("finngen_R10_", phenotype$phenocode)

# 加载包
library(dplyr)

# 合并 cell 数据
combined_cell <- combined_cell %>%
  left_join(phenotype[,1:5], 
            by = c("id.outcome" = "phenocode"))

# 合并 gene 数据
combined_gene <- combined_gene %>%
  left_join(phenotype[,1:5], 
            by = c("id.outcome" = "phenocode"))

library(TwoSampleMR)
combined_gene <- generate_odds_ratios(combined_gene)
combined_cell <- generate_odds_ratios(combined_cell)

### ----------------------------------------------------------------------------
###            多重校正
### ----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(purrr)

# 定义分列校正函数（针对每个变量单独分组校正）
apply_columnwise_correction <- function(data, group_columns) {
  # 遍历每个分组变量
  for (col in group_columns) {
    # 生成校正后的列名
    bh_col <- paste0("pval_BH_", col)
    bonf_col <- paste0("pval_Bonferroni_", col)
    
    # 按当前列分组并校正
    data <- data %>%
      group_by(across(all_of(col))) %>%
      mutate(
        !!bh_col := p.adjust(pval, method = "BH"),
        !!bonf_col := p.adjust(pval, method = "bonferroni")
      ) %>%
      ungroup()
  }
  return(data)
}

# 筛选 nsnp > 3 的行并分档 num_cases
cell_processed <- na.omit(combined_cell[,-5])

colnames(cell_processed)[8] <- c("pval")
# 定义需要独立校正的列
cell_group_columns <- c("id.outcome", "exposure", "SNP", "category")

# 应用分列校正
cell_corrected <- apply_columnwise_correction(cell_processed, cell_group_columns)

# 查看新增列
names(cell_corrected)
# 输出示例：
# [1] "pval_BH_outcome"          "pval_Bonferroni_outcome"  "pval_BH_exposure"        
# [4] "pval_Bonferroni_exposure" "pval_BH_method"           ...

library(openxlsx)

# ==================== 定义通用保存函数 ====================
save_multiple_sheets <- function(data, filename) {
  # 创建新工作簿
  wb <- createWorkbook()
  
  # 添加原始数据分页
  addWorksheet(wb, "Original_Data")
  writeData(wb, sheet = "Original_Data", x = data)
  
  # 获取所有p值列（包括原始pval和校正后的所有p值）
  pval_cols <- grep("^pval", names(data), value = TRUE)
  
  # 遍历每个p值列筛选显著结果
  for (col in pval_cols) {
    # 筛选当前p值列 <0.05 的行
    sig_data <- data[data[[col]] < 0.05, ]
    
    # 仅当存在显著结果时保存
    if (nrow(sig_data) > 0) {
      # 生成分页名称（确保符合Excel规范）
      sheet_name <- substr(col, 1, 31)  # Excel分页名最长31字符
      sheet_name <- gsub("\\.", "_", sheet_name)  # 替换特殊字符
      
      addWorksheet(wb, sheetName = sheet_name)
      writeData(wb, sheet = sheet_name, x = sig_data)
    }
  }
  
  # 保存工作簿
  saveWorkbook(wb, file = filename, overwrite = TRUE)
}

# ==================== 执行保存操作 ====================
# 保存 cell 数据
save_multiple_sheets(cell_corrected, "cell_results Leave-One-Out.xlsx")

### ----------------------------------------------------------------------------
###            upset
### ----------------------------------------------------------------------------

library(dplyr)
library(UpSetR)
library(openxlsx)
library(tibble)

# Cell数据
cell_pval_cols <- grep("^pval_", names(cell_corrected), value = TRUE)

calculate_effect <- function(or, lci, uci) {
  if (or > 1 & lci > 1) {
    return("Risk Factor")
  } else if (or < 1 & uci < 1) {
    return("Protective Factor")
  } else {
    return("Non-Significant")
  }
}

library(purrr)
library(dplyr)

cell_results <- map_dfr(cell_pval_cols, ~{
  cell_corrected %>%
    filter(!!sym(.x) < 0.05) %>%
    mutate(
      pval_type = .x,
      effect = pmap_chr(list(or, or_lci95, or_uci95), calculate_effect)
    ) %>%
    select(exposure, id.outcome, SNP, effect, pval_type)
}) %>%
  distinct()  # 去重



cell_results

library(dplyr)
library(tidyr)

# 处理cell_results
cell_final <- cell_results %>%
  mutate(combo = paste(exposure, id.outcome, SNP, effect,  sep = " || ")) %>%
  select(pval_type, combo) %>%
  group_by(pval_type) %>%
  mutate(row_id = row_number()) %>%
  ungroup() %>%
  pivot_wider(
    names_from = pval_type,
    values_from = combo,
    values_fill = list(combo = NA_character_)
  ) %>%
  select(-row_id) %>%
  slice(1:max(map_int(., ~sum(!is.na(.x)))))  # 动态计算最大有效行数

# 查看结果
print(cell_final, n = 3, width = 120)

# 转换为UpSetR兼容格式（以gene_final为例）
library(UpSetR)

# 将宽格式转换回长格式（恢复pval_type列）
cell_upset_data <- cell_final %>%
  pivot_longer(
    cols = everything(),
    names_to = "pval_type",
    values_to = "combo",
    values_drop_na = TRUE
  ) %>%
  distinct()  # 去重处理

cell_list <- split(cell_upset_data$combo, cell_upset_data$pval_type)

# 细胞数据处理 ---------------------------------------------------------------
# 创建输出目录（如果不存在）[3](@ref)
if(!dir.exists("upset")) dir.create("upset")
# 生成图表（重复基因流程）
pdf("upset/cell_upset  Leave-One-Out.pdf", width=8, height=6)
upset(fromList(cell_list),
      nsets = 99, nintersects = 20,
      mainbar.y.label = "Significant Combinations",
      sets.x.label = "p-value Types Count")
dev.off()

png("upset/cell_upset  Leave-One-Outp.png", width=2400, height=1800, res=300)
upset(fromList(cell_list),
      nsets = 99, nintersects = 20,
      mainbar.y.label = "Significant Combinations",
      sets.x.label = "p-value Types Count")
dev.off()

# 导出交集矩阵
cell_intersect_matrix <- as.data.frame.matrix(table(cell_upset_data$combo, cell_upset_data$pval_type))
write.xlsx(cell_intersect_matrix, "upset/cell_intersections  Leave-One-Out.xlsx",
           colNames = TRUE, rowNames = TRUE)



