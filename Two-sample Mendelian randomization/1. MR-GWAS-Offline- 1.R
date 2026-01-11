####----------------------前期准备----------------------------------
# 合并重复的库加载（避免重复加载）
library(TwoSampleMR)
library(data.table)
library(tidyverse)
#library(ieugwasr)
library(MRInstruments)
library(mr.raps)
library(MendelianRandomization)
library(vroom)
rm(list = ls())

target_outcome <- c("KIRC") ####根据需求替换

target_gene <- c("CRHBP|MFSD4|MPP7|UCN|CRH|HSP90AA1|PIK3CB|PIK3CG|PIK3C3|PIK3R4|PIK3C2B|PIK3CA|ATG14|MTOR|BECN1|AKT1|KIT|KDR|RICTOR|FOXO1|VEGFA|NOS3|MDM2") ####根据需求替换,注意必须要1个基因以上才行

#data(gwas_catalog)
#bmi <- subset(gwas_catalog, Phenotype=="Body mass index")
#bmi <- format_data(bmi)

# 设置固定工作目录路径（根据你的需求修改此处）
target_dir <- "E:/20241226 Yang Wang/20250808 IJS Revise/mr"  # 注意使用正斜杠或双反斜杠
setwd(target_dir)

### 预加载细胞的数据
load("01immune.cell.5e-08.RData")
exposure_dat <- immu.cell.5e.08

# 自动获取ZIP文件并处理路径
list.files(pattern = "\\.gz$")

####这个地方要改一下序号----------
zip_file <- list.files(pattern = "\\.gz$")[1]          # 取当前目录第一个ZIP文件
stopifnot("未找到gz文件" = !is.na(zip_file))          # 确保文件存在

library(tools)
zip_name <- file_path_sans_ext(zip_file)               # 自动提取文件夹名称

# 0. 清理或创建文件夹（保留文件夹结构）
if (dir.exists(zip_name)) {
  # 清空文件夹内容但保留结构（关键修改）
  unlink(file.path(zip_name, "*"), recursive = TRUE, force = TRUE)
  message("已清空文件夹内容：", zip_name)
} else {
  dir.create(zip_name)
  message("已创建新文件夹：", zip_name)
}

# 1. 创建新文件夹
dir.create(zip_name)                                   # 新建空文件夹

# 结局因素的处理
#读取本地结局
outcome_dat<-vroom(paste0("./",zip_name,".gz"))

# 2. 设置新工作目录
setwd(zip_name)                                        # 切换到解压文件夹
message("当前工作目录已设置为：", getwd())

#head(outcome_dat,10)

library(dplyr)

# 重命名核心列
outcome_dat <- outcome_dat %>%
  rename(
    SNP = rsids,          # SNP标识符
    effect_allele = alt,  # 效应等位基因（如ALT）
    other_allele = ref,   # 非效应等位基因（如REF）
    beta = beta,          # 效应值（原列名已正确）
    se = sebeta,          # 标准误（原列名已正确）
    pval = pval,          # P值（原列名已正确）
    eaf = af_alt          #可选：添加效应等位基因频率列（如果存在）
  )

library(dplyr)
# 方法2：直接使用distinct()（需确保SNP唯一性）
outcome_dat_cleaned <- outcome_dat %>%
  distinct(SNP, .keep_all = TRUE)  # 但无法控制p值/beta优先级

# 验证去重结果
sum(duplicated(outcome_dat_cleaned$SNP))  # 应为0

outcome_dat <- outcome_dat_cleaned %>% format_data(type = "outcome")
outcome_dat$id.outcome <- zip_name
saveRDS(outcome_dat,"outcome_dat.RDS")

####----------------------基因的分析----------------------------------
### outcome 处理好了
#rm(list = ls())
#outcome_dat <- readRDS("outcome_dat.RDS")
#data("drug_interactions")
data("gtex_eqtl")
summary(as.factor(gtex_eqtl$tissue))

# 筛选CRHBP、MFSD4、MPP7基因
exposure_dat <- subset(gtex_eqtl, 
                       grepl(target_gene, gene_name, ignore.case = TRUE))

summary(as.factor(gtex_eqtl$tissue))

# 在调用format_data前生成唯一ID
exposure_dat$unique_id <- paste(exposure_dat$gene_name, 
                                exposure_dat$tissue, 
                                sep = "_")

head(exposure_dat)

exposure_dat$exposure <- paste(exposure_dat$gene_name,"in",exposure_dat$tissue)
head(exposure_dat)

# 3. 检查每个基因-组织的SNP数量
snp_counts <- exposure_dat %>%
  group_by(gene_name, tissue) %>%
  summarise(N_SNPs = n())

exposure_dat <- format_data(
  exposure_dat,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  id_col = "unique_id",  # 新增的唯一ID列
  eaf_col = "af_alt",     # 假设原始数据包含效应等位基因频率
  gene_col = "gene_name",
  samplesize_col = "n"  # 样本量列（如有）
)

head(exposure_dat)

head(exposure_dat[, c("SNP", "id.exposure", "beta.exposure")])


### harmonise ------------------------------------------------------------------
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
saveRDS(dat,"dat gene.RDS")

summary(as.factor(dat$id.exposure))

library(tidyr)

# 步骤1：拆分列
dat <- dat %>%
  separate(
    col = id.exposure,         # 要拆分的列
    into = c("id.exposure", "exposure"), # 新列名：基因 + 组织
    sep = "_",                 # 以下划线作为分隔符
    extra = "merge",           # 将多余分隔符合并到最后一列（保留组织名称中的空格/下划线）
    remove = TRUE              # 删除原始列
  )

summary(as.factor(dat$id.exposure))
summary(as.factor(dat$exposure))

### MRPRESSO ------------------------------------------------------------------
library(MRPRESSO)

# 初始化结果容器（添加字符串处理参数）
egger_intercept_df <- data.frame(stringsAsFactors = FALSE)
presso_global_df <- data.frame(stringsAsFactors = FALSE)
presso_corrected_df <- data.frame(stringsAsFactors = FALSE)

# 按每个基因-组织组合分组分析
for (id_exp in unique(dat$id.exposure)) {
  dat_sub <- dat[dat$id.exposure == id_exp, ]
  num_snps <- nrow(dat_sub)
  
  # 跳过SNP<3的情况
  if (num_snps < 3) {
    warning(paste("跳过", id_exp, ": SNP < 3，无法进行MR-Egger和MR-PRESSO"))
    next
  }
  
  ### 1. MR-Egger截距检验（添加错误处理）
  pleiotropy_test <- tryCatch({
    mr_pleiotropy_test(dat_sub)
  }, error = function(e) {
    warning(paste("MR-Egger失败于", id_exp, ":", e$message))
    return(list(egger_intercept = NA, se = NA, pval = NA))
  })
  
  egger_intercept_df <- rbind(egger_intercept_df, 
                              data.frame(
                                id.exposure = id_exp,
                                egger_intercept = pleiotropy_test$egger_intercept,
                                se = pleiotropy_test$se,
                                pval = pleiotropy_test$pval,
                                stringsAsFactors = FALSE
                              ))
  
  ### 2. MR-PRESSO分析（仅当SNP≥10时运行）
  if (num_snps >= 10) {
    presso <- tryCatch({
      MRPRESSO::mr_presso(
        BetaOutcome = "beta.outcome",
        BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome",
        SdExposure = "se.exposure",
        data = dat_sub,
        NbDistribution = 1000
      )
    }, error = function(e) {
      warning(paste("MR-PRESSO失败于", id_exp, ":", e$message))
      return(NULL)
    })
    
    # 保存PRESSO结果
    if (!is.null(presso)) {
      # 全局检验结果
      global_test <- presso$`MR-PRESSO results`$`Global Test`
      presso_global_df <- rbind(presso_global_df, 
                                data.frame(
                                  id.exposure = id_exp,
                                  RSSobs = global_test$RSSobs,
                                  Pvalue = global_test$Pvalue,
                                  stringsAsFactors = FALSE
                                ))
      
      # 校正后的MR结果
      corrected_res <- presso$`Main MR results`[2, ]
      presso_corrected_df <- rbind(presso_corrected_df, 
                                   data.frame(
                                     id.exposure = id_exp,
                                     Causal_Estimate = ifelse(is.null(corrected_res$`Causal Estimate`), 
                                                              NA, corrected_res$`Causal Estimate`),
                                     Sd = ifelse(is.null(corrected_res$Sd), NA, corrected_res$Sd),
                                     T_stat = ifelse(is.null(corrected_res$`T-stat`), NA, corrected_res$`T-stat`),
                                     P_value = ifelse(is.null(corrected_res$`P-value`), NA, corrected_res$`P-value`),
                                     stringsAsFactors = FALSE
                                   ))
    } else {
      # 记录PRESSO失败的情况
      presso_global_df <- rbind(presso_global_df, 
                                data.frame(
                                  id.exposure = id_exp,
                                  RSSobs = NA,
                                  Pvalue = NA,
                                  stringsAsFactors = FALSE
                                ))
      presso_corrected_df <- rbind(presso_corrected_df, 
                                   data.frame(
                                     id.exposure = id_exp,
                                     Causal_Estimate = NA,
                                     Sd = NA,
                                     T_stat = NA,
                                     P_value = NA,
                                     stringsAsFactors = FALSE
                                   ))
    }
  }
}

# 打印结果
print("MR-Egger截距检验结果:")
print(egger_intercept_df)

print("MR-PRESSO全局检验结果:")
print(presso_global_df)

print("MR-PRESSO校正后结果:")
print(presso_corrected_df)

### MR 分析--------------------------------------------------------------------
mr_method_list <- mr_method_list()

print(mr_method_list$obj)

str(dat)
# MR 分析
result <- mr(dat)#默认的四种方法，mr_method_list可以查看所有方法

saveRDS(result,"gene_mr_result.RDS")

# 计算OR值
result <- generate_odds_ratios(result)

result$outcome <- zip_name

head(result)

# 异质性检验
heterogeneity <- mr_heterogeneity(dat)
#若Q检验P < 0.05：存在显著异质性，提示潜在水平多效性。

# 水平多效性检验
pleiotropy <- mr_pleiotropy_test(dat)

# 散点图
p1 <- mr_scatter_plot(result, dat)

# 森林图
result_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(result_single)

# 留一图
result_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(result_loo)

# 漏斗图
result_single <- mr_singlesnp(dat)
p4 <- mr_funnel_plot(result_single)

# 加载必要包
library(TwoSampleMR)
library(ggplot2)

# 创建总文件夹
if (!dir.exists("gene")) dir.create("gene")

# 获取所有唯一的基因ID
gene_ids <- unique(dat$id.exposure)

# 循环处理每个基因
for (gene_id in gene_ids) {
  # 创建基因专属文件夹
  gene_dir <- file.path("gene", gene_id)
  if (!dir.exists(gene_dir)) dir.create(gene_dir)
  
  # 筛选当前基因数据
  dat_gene <- dat[dat$id.exposure == gene_id, ]
  
  # 跳过SNP不足的情况
  if (nrow(dat_gene) < 3) {
    warning(paste("跳过", gene_id, ": SNP < 3，无法分析"))
    next
  }
  
  # 执行MR分析
  result_gene <- mr(dat_gene)
  
  # === 1. 散点图 ===
  p_scatter <- mr_scatter_plot(result_gene, dat_gene)[[1]] +
    labs(title = paste("Scatter Plot:", gene_id))
  
  # === 2. 森林图 ===
  result_single <- mr_singlesnp(dat_gene)
  p_forest <- mr_forest_plot(result_single)[[1]] +
    ggtitle(paste("Forest Plot:", gene_id))
  
  # === 3. 留一图 ===
  result_loo <- mr_leaveoneout(dat_gene)
  p_loo <- mr_leaveoneout_plot(result_loo)[[1]] +
    ggtitle(paste("Leave-One-Out Plot:", gene_id))
  
  # === 4. 漏斗图 ===
  p_funnel <- mr_funnel_plot(result_single)[[1]] +
    ggtitle(paste("Funnel Plot:", gene_id))
  
  # === 保存4种格式的图形 ===
  plot_list <- list(
    scatter = p_scatter,
    forest = p_forest,
    loo = p_loo,
    funnel = p_funnel
  )
  
  for (plot_name in names(plot_list)) {
    file_prefix <- file.path(gene_dir, paste0(gene_id, "_", plot_name))
    
    # JPG格式
    ggsave(paste0(file_prefix, ".jpg"), plot_list[[plot_name]], 
           width = 10, height = 8, dpi = 300)
    
    # PNG格式
    ggsave(paste0(file_prefix, ".png"), plot_list[[plot_name]],
           width = 10, height = 8, dpi = 300)
    
    # PDF格式（矢量图）
    ggsave(paste0(file_prefix, ".pdf"), plot_list[[plot_name]],
           width = 10, height = 8, device = cairo_pdf)
    
    # TIFF格式
    ggsave(paste0(file_prefix, ".tiff"), plot_list[[plot_name]],
           width = 10, height = 8, dpi = 300, compression = "lzw")
  }
  
  # 进度提示
  message(paste("已完成基因:", gene_id, "| 保存位置:", gene_dir))
}

####----------------------数据及结果保存---------------------------------
# 安装并加载writexl（无Java依赖，不生成临时文件）
if (!require("writexl")) install.packages("writexl")
library(writexl)

# 创建cell_result文件夹（如果不存在）
if (!dir.exists("gene_result")) {
  dir.create("gene_result")
  message("文件夹创建成功: gene_result/")
} else {
  message("文件夹已存在: gene_result/")
}

# 文件1：Egger+PRESSO结果
write_xlsx(
  list(
    Egger_Intercept = egger_intercept_df,
    PRESSO_Global = presso_global_df,
    PRESSO_Corrected = presso_corrected_df
  ),
  path = "gene_result/MR_PRESSO_AND_Egger.xlsx"
)

# 文件2：主分析结果
write_xlsx(
  list(
    `Main Results` = result,
    `Single SNP` = result_single,
    `Leave-One-Out` = result_loo
  ),
  path = "gene_result/MR_Main_Results.xlsx"
)

# 文件3：敏感性分析
write_xlsx(
  list(
    heterogeneity = heterogeneity,
    pleiotropy = pleiotropy
  ),
  path = "gene_result/Sensitivity_Analysis.xlsx"
)

# 确认输出完成
message("Results saved to: ", file.path(getwd(), "MR_gene_Results.xlsx"))


####----------------------细胞的分析----------------------------------

exposure_dat <- immu.cell.5e.08

### harmonise ------------------------------------------------------------------
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
saveRDS(dat,"dat cell.RDS")

rm(list = ls())

####这个地方要改一下序号----------
# 设置固定工作目录路径（根据你的需求修改此处）
target_dir <- "E:/20241226 Yang Wang/20250808 IJS Revise/mr"  # 注意使用正斜杠或双反斜杠
setwd(target_dir)

message("当前工作目录已设置为：", getwd())

# 自动获取ZIP文件并处理路径
list.files(pattern = "\\.gz$")
zip_file <- list.files(pattern = "\\.gz$")[1]          # 取当前目录第一个ZIP文件
stopifnot("未找到gz文件" = !is.na(zip_file))          # 确保文件存在

library(tools)
zip_name <- file_path_sans_ext(zip_file)               # 自动提取文件夹名称

# 2. 设置新工作目录
setwd(zip_name)                                        # 切换到解压文件夹

dat <- readRDS("dat cell.RDS")
summary(as.factor(dat$id.exposure))

library(tidyr)

### MRPRESSO ------------------------------------------------------------------
library(MRPRESSO)

# 初始化结果容器（添加字符串处理参数）
egger_intercept_df <- data.frame(stringsAsFactors = FALSE)
presso_global_df <- data.frame(stringsAsFactors = FALSE)
presso_corrected_df <- data.frame(stringsAsFactors = FALSE)

# 按每个基因-组织组合分组分析
for (id_exp in unique(dat$id.exposure)) {
  dat_sub <- dat[dat$id.exposure == id_exp, ]
  num_snps <- nrow(dat_sub)
  
  # 跳过SNP<3的情况
  if (num_snps < 3) {
    warning(paste("跳过", id_exp, ": SNP < 3，无法进行MR-Egger和MR-PRESSO"))
    next
  }
  
  ### 1. MR-Egger截距检验（添加错误处理）
  pleiotropy_test <- tryCatch({
    mr_pleiotropy_test(dat_sub)
  }, error = function(e) {
    warning(paste("MR-Egger失败于", id_exp, ":", e$message))
    return(list(egger_intercept = NA, se = NA, pval = NA))
  })
  
  egger_intercept_df <- rbind(egger_intercept_df, 
                              data.frame(
                                id.exposure = id_exp,
                                egger_intercept = pleiotropy_test$egger_intercept,
                                se = pleiotropy_test$se,
                                pval = pleiotropy_test$pval,
                                stringsAsFactors = FALSE
                              ))
  
  ### 2. MR-PRESSO分析（仅当SNP≥10时运行）
  if (num_snps >= 10) {
    presso <- tryCatch({
      MRPRESSO::mr_presso(
        BetaOutcome = "beta.outcome",
        BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome",
        SdExposure = "se.exposure",
        data = dat_sub,
        NbDistribution = 1000
      )
    }, error = function(e) {
      warning(paste("MR-PRESSO失败于", id_exp, ":", e$message))
      return(NULL)
    })
    
    # 保存PRESSO结果
    if (!is.null(presso)) {
      # 全局检验结果
      global_test <- presso$`MR-PRESSO results`$`Global Test`
      presso_global_df <- rbind(presso_global_df, 
                                data.frame(
                                  id.exposure = id_exp,
                                  RSSobs = global_test$RSSobs,
                                  Pvalue = global_test$Pvalue,
                                  stringsAsFactors = FALSE
                                ))
      
      # 校正后的MR结果
      corrected_res <- presso$`Main MR results`[2, ]
      presso_corrected_df <- rbind(presso_corrected_df, 
                                   data.frame(
                                     id.exposure = id_exp,
                                     Causal_Estimate = ifelse(is.null(corrected_res$`Causal Estimate`), 
                                                              NA, corrected_res$`Causal Estimate`),
                                     Sd = ifelse(is.null(corrected_res$Sd), NA, corrected_res$Sd),
                                     T_stat = ifelse(is.null(corrected_res$`T-stat`), NA, corrected_res$`T-stat`),
                                     P_value = ifelse(is.null(corrected_res$`P-value`), NA, corrected_res$`P-value`),
                                     stringsAsFactors = FALSE
                                   ))
    } else {
      # 记录PRESSO失败的情况
      presso_global_df <- rbind(presso_global_df, 
                                data.frame(
                                  id.exposure = id_exp,
                                  RSSobs = NA,
                                  Pvalue = NA,
                                  stringsAsFactors = FALSE
                                ))
      presso_corrected_df <- rbind(presso_corrected_df, 
                                   data.frame(
                                     id.exposure = id_exp,
                                     Causal_Estimate = NA,
                                     Sd = NA,
                                     T_stat = NA,
                                     P_value = NA,
                                     stringsAsFactors = FALSE
                                   ))
    }
  }
}

# 打印结果
print("MR-Egger截距检验结果:")
print(egger_intercept_df)

print("MR-PRESSO全局检验结果:")
print(presso_global_df)

print("MR-PRESSO校正后结果:")
print(presso_corrected_df)

### MR 分析--------------------------------------------------------------------
mr_method_list <- mr_method_list()

print(mr_method_list$obj)

str(dat)
# MR 分析
result <- mr(dat)#默认的四种方法，mr_method_list可以查看所有方法

saveRDS(result,"cell_mr_result.RDS")

# 计算OR值
result <- generate_odds_ratios(result)

result$outcome <- zip_name

head(result)

# 异质性检验
heterogeneity <- mr_heterogeneity(dat)
#若Q检验P < 0.05：存在显著异质性，提示潜在水平多效性。

# 水平多效性检验
pleiotropy <- mr_pleiotropy_test(dat)

# 散点图
p1 <- mr_scatter_plot(result, dat)

# 森林图
result_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(result_single)

# 留一图
result_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(result_loo)

# 漏斗图
p4 <- mr_funnel_plot(result_single)

# 加载必要包
library(TwoSampleMR)
library(ggplot2)

# 创建总文件夹
if (!dir.exists("cell")) dir.create("cell")

# 获取所有唯一的基因ID
gene_ids <- unique(dat$id.exposure)

# 循环处理每个基因
for (gene_id in gene_ids) {
  # 创建基因专属文件夹
  gene_dir <- file.path("cell", gene_id)
  if (!dir.exists(gene_dir)) dir.create(gene_dir)
  
  # 筛选当前基因数据
  dat_gene <- dat[dat$id.exposure == gene_id, ]
  
  # 跳过SNP不足的情况
  if (nrow(dat_gene) < 3) {
    warning(paste("跳过", gene_id, ": SNP < 3，无法分析"))
    next
  }
  
  # 执行MR分析
  result_gene <- mr(dat_gene)
  
  # === 1. 散点图 ===
  p_scatter <- mr_scatter_plot(result_gene, dat_gene)[[1]] +
    labs(title = paste("Scatter Plot:", gene_id))
  
  # === 2. 森林图 ===
  result_single <- mr_singlesnp(dat_gene)
  p_forest <- mr_forest_plot(result_single)[[1]] +
    ggtitle(paste("Forest Plot:", gene_id))
  
  # === 3. 留一图 ===
  result_loo <- mr_leaveoneout(dat_gene)
  p_loo <- mr_leaveoneout_plot(result_loo)[[1]] +
    ggtitle(paste("Leave-One-Out Plot:", gene_id))
  
  # === 4. 漏斗图 ===
  p_funnel <- mr_funnel_plot(result_single)[[1]] +
    ggtitle(paste("Funnel Plot:", gene_id))
  
  # === 保存4种格式的图形 ===
  plot_list <- list(
    scatter = p_scatter,
    forest = p_forest,
    loo = p_loo,
    funnel = p_funnel
  )
  
  for (plot_name in names(plot_list)) {
    file_prefix <- file.path(gene_dir, paste0(gene_id, "_", plot_name))
    
    # JPG格式
    ggsave(paste0(file_prefix, ".jpg"), plot_list[[plot_name]], 
           width = 10, height = 8, dpi = 300)
    
    # PNG格式
    ggsave(paste0(file_prefix, ".png"), plot_list[[plot_name]],
           width = 10, height = 8, dpi = 300)
    
    # PDF格式（矢量图）
    #ggsave(paste0(file_prefix, ".pdf"), plot_list[[plot_name]],
    #       width = 10, height = 8, device = cairo_pdf)
    
    # TIFF格式
    ggsave(paste0(file_prefix, ".tiff"), plot_list[[plot_name]],
           width = 10, height = 8, dpi = 300, compression = "lzw")
  }
  
  # 进度提示
  message(paste("已完成细胞:", gene_id, "| 保存位置:", gene_dir))
}



####----------------------数据及结果保存---------------------------------
# 安装并加载writexl（无Java依赖，不生成临时文件）
if (!require("writexl")) install.packages("writexl")
library(writexl)

# 创建cell_result文件夹（如果不存在）
if (!dir.exists("cell_result")) {
  dir.create("cell_result")
  message("文件夹创建成功: cell_result/")
} else {
  message("文件夹已存在: cell_result/")
}

# 文件1：Egger+PRESSO结果
write_xlsx(
  list(
    Egger_Intercept = egger_intercept_df,
    PRESSO_Global = presso_global_df,
    PRESSO_Corrected = presso_corrected_df
  ),
  path = "cell_result/MR_cell_PRESSO_AND_Egger.xlsx"
)

# 文件2：主分析结果
write_xlsx(
  list(
    `Main Results` = result,
    `Single SNP` = result_single,
    `Leave-One-Out` = result_loo
  ),
  path = "cell_result/MR_cell_Main_Results.xlsx"
)

# 文件3：敏感性分析
write_xlsx(
  list(
    heterogeneity = heterogeneity,
    pleiotropy = pleiotropy
  ),
  path = "cell_result/Sensitivity_Analysis.xlsx"
)

# 确认输出完成
message("Results saved to: ", file.path(getwd(), "MR_cell_Results.xlsx"))
