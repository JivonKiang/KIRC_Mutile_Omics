rm(list = ls())

# 加载必要包
library(tidyverse)
library(rstatix)
library(ggpubr)
library(ggrain)
library(openxlsx)
library(limma)

# 读取原始数据
data <- read.xlsx("qPCR.xlsx")

# 修复批次校正函数
batch_correct_qpcr <- function(data) {
  # 1. 创建长格式数据
  long_data <- data %>%
    pivot_longer(
      cols = c(KIRC, CTL),
      names_to = "Group",
      values_to = "Expression"
    ) %>%
    mutate(
      unique_id = paste(Times, Sepecies, Group, sample_id, sep = "|")
    )
  
  # 2. 按物种分组处理
  species_list <- unique(long_data$Sepecies)
  corrected_list <- list()
  
  for (species in species_list) {
    # 提取当前物种数据
    species_data <- long_data %>% 
      filter(Sepecies == species) %>%
      select(unique_id, Times, Group, Gene, Expression) %>%
      pivot_wider(
        names_from = Gene,
        values_from = Expression
      )
    
    # 准备表达矩阵
    expr_mat <- as.matrix(select(species_data, -c(unique_id, Times, Group)))
    rownames(expr_mat) <- species_data$unique_id
    
    # 设置批次和设计矩阵
    batch <- factor(species_data$Times)
    design <- model.matrix(~ Group, data = species_data)
    
    # 批次校正
    corrected_mat <- removeBatchEffect(
      t(expr_mat),  # 转置为基因×样本
      batch = batch,
      design = design
    )
    
    # 转换回长格式
    corrected_long <- as.data.frame(t(corrected_mat)) %>%
      rownames_to_column("unique_id") %>%
      pivot_longer(
        cols = -unique_id,
        names_to = "Gene",
        values_to = "Expression_corrected"
      )
    
    # 合并元数据
    meta_data <- species_data %>% select(unique_id, Times, Group)
    corrected_list[[species]] <- left_join(corrected_long, meta_data, by = "unique_id") %>%
      mutate(Sepecies = species) %>%
      separate(unique_id, into = c("Times", "Sepecies", "Group", "sample_id"), sep = "\\|")
  }
  
  # 合并所有物种数据
  bind_rows(corrected_list) %>%
    select(sample_id, Sepecies, Times, Group, Gene, Expression = Expression_corrected)
}

# 执行批次校正
data_corrected <- batch_correct_qpcr(data)

# 2. 分组统计检验
# 统计检验（使用校正后数据）
stat_results <- data_corrected %>%
  group_by(Sepecies, Gene) %>%
  wilcox_test(Expression ~ Group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(
    p_label = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# 准备绘图数据
plot_data <- data_corrected %>%
  filter(!is.na(Expression)) %>%  # 移除NA值
  mutate(Group = factor(Group, levels = c("CTL", "KIRC")))

# 创建分面云雨图函数
create_rainplot <- function(species_name, gene_name) {
  # 提取当前统计结果
  current_stat <- stat_results %>%
    filter(Sepecies == species_name, Gene == gene_name)
  
  # 计算y轴位置
  y_max <- max(plot_data$Expression[plot_data$Sepecies == species_name & 
                                      plot_data$Gene == gene_name], na.rm = TRUE)*2
  y_min <- min(plot_data$Expression[plot_data$Sepecies == species_name & 
                                      plot_data$Gene == gene_name], na.rm = TRUE)
  y_range <- y_max - y_min
  
  # 生成云雨图
  p <- plot_data %>%
    filter(Sepecies == species_name, Gene == gene_name) %>%
    ggplot(aes(x = Group, y = Expression*2.1, fill = Group)) +
    ggrain::geom_rain(
      alpha = 0.7, 
      point.args = list(size = 2, alpha = 0.6),
      boxplot.args = list(width = 0.12, outlier.shape = NA),
      violin.args = list(alpha = 0.4)
    ) +
    geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.7) +
    scale_fill_manual(values = c("CTL" = "#4daf4a", "KIRC" = "#e41a1c")) +
    labs(
      title = paste(species_name, "-", gene_name),
      y = "Corrected Expression (log2)"
    ) +
    theme_classic(base_size = 12) +  # 增大基础字号
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 12),
      legend.position = "none",
      panel.border = element_rect(fill = NA, color = "black", linewidth = 1)
    ) +
    # 添加统计标注 - 使用p值而非显著性标记
    geom_signif(
      annotations = paste0("p = ", formatC(current_stat$p.adj, 
                                           #format = "e", 
                                           digits = 2)),
      y_position = y_max*1.8,
      xmin = 1, 
      xmax = 2, 
      textsize = 4.5, 
      tip_length = 0.01,
      vjust = 1.5                # 文本下移避免重叠[7](@ref)
    )
  
  return(p)
}

# 生成四个子图
plots <- list(
  hs_crhbp = create_rainplot("Homo sapiens", "CRHBP"),
  hs_ucn2 = create_rainplot("Homo sapiens", "UCN2"),
  mouse_crhbp = create_rainplot("Mouse", "CRHBP"),
  mouse_ucn2 = create_rainplot("Mouse", "UCN2")
)

# 合并子图并保存
final_plot <- ggarrange(
  plots$hs_crhbp, plots$hs_ucn2, 
  plots$mouse_crhbp, plots$mouse_ucn2,
  ncol = 2, nrow = 2,
  labels = c("A", "B", "C", "D"),
  font.label = list(size = 16, face = "bold")
)

# 多格式保存（调整输出尺寸）
formats <- c("jpg", "png", "tiff", "pdf")
for (fmt in formats) {
  ggsave(
    filename = paste0("RaincloudPlot_", Sys.Date(), ".", fmt),
    plot = final_plot,
    width = 5,
    height = 6,
    dpi = ifelse(fmt %in% c("tiff", "pdf"), 300, 600)
  )
}
