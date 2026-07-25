#!/usr/bin/env Rscript

# =============================================================================
# 脚本功能: 绘制堆叠柱状图（数值版 + 百分比版）
# 运行环境: mamba activate r452
# 输入文件格式: 第一列为分组信息(绘图x轴, 需要列名, 对应group_col), 后面列为特征值数值信息(列名对应value_cols)
# =============================================================================

# 加载所需包
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
source("/data/med-wangcq/01CondaEnv/02Git_repo/00MyGit_wchenqi/Multi_omics/RNA-seq/00Functions/ColorPanel.r")

# ============================================================
# 主函数：绘制堆叠柱状图
# ============================================================
plot_stacked_bar <- function(infile, 
                             group_col = "Group",
                             value_cols = "default",
                             colors = "default",
                             outdir = "./",
                             prefix = "DiffExpress",
                             ymax_counts = NULL,      # 数值图 y 轴最大值
                             yby_counts = NULL) {     # 数值图 y 轴刻度间隔
    
    # ---- 1. 读取数据 ----
    data <- fread(infile) %>% as.data.frame()
    cat("\n=== 原始数据 ===\n")
    print(data)
    
    # 处理 value_cols
    if(value_cols == "default"){
        value_cols <- colnames(data)[2:ncol(data)]
    } else {
        value_cols <- unlist(strsplit(value_cols, ","))
    }
    cat("\n=== 数值列 ===\n")
    print(value_cols)
    
    # ---- 2. 计算总数和百分比 ----
    # 计算每行总数
    data$Total <- rowSums(data[, value_cols, drop = FALSE])
    colnames(data)[1] <- "Group"
    data$Group <- factor(data$Group, levels = unique(data$Group))
    
    # 计算百分比（确保是小数形式）
    data_percent_raw <- data
    for(col in value_cols) {
        data_percent_raw[[col]] <- data_percent_raw[[col]] / data_percent_raw$Total * 100
    }
    
    cat("\n=== 百分比数据 ===\n")
    print(data_percent_raw)
    
    # ---- 3. 颜色设置 ----
    if(colors == "default"){
        if(exists("cols")) {
            colors_use <- cols[1:length(value_cols)]
        } else {
            library(RColorBrewer)
            colors_use <- brewer.pal(max(3, length(value_cols)), "Set2")[1:length(value_cols)]
        }
        names(colors_use) <- value_cols
    } else if(is.character(colors) && length(colors) == 1 && grepl(",", colors)) {
        colors_use <- unlist(strsplit(colors, ","))
        if(length(colors_use) != length(value_cols)) {
            warning("颜色数量不匹配，使用默认颜色")
            if(exists("cols")) {
                colors_use <- cols[1:length(value_cols)]
            } else {
                colors_use <- brewer.pal(max(3, length(value_cols)), "Set2")[1:length(value_cols)]
            }
        }
        names(colors_use) <- value_cols
    } else {
        colors_use <- colors
        if(length(colors_use) != length(value_cols)) {
            stop("颜色数量与数值列数量不一致")
        }
        if(is.null(names(colors_use))) {
            names(colors_use) <- value_cols
        }
    }
    
    cat("\n=== 颜色映射 ===\n")
    print(colors_use)
    
    # ---- 4. 数据转换为长格式 ----
    # 数值数据长格式
    data_long <- data %>%
        pivot_longer(cols = all_of(value_cols), names_to = "Category", values_to = "Count")
    
    # 百分比数据长格式（使用计算好的百分比数据）
    data_percent_long <- data_percent_raw %>%
        pivot_longer(cols = all_of(value_cols), names_to = "Category", values_to = "Percentage")
    
    # 确保 Category 是因子，保持顺序
    data_long$Category <- factor(data_long$Category, levels = value_cols)
    data_percent_long$Category <- factor(data_percent_long$Category, levels = value_cols)
    
    # ---- 5. 计算最大总数用于y轴调整 ----
    max_total <- max(data$Total)
    
    # ---- 6. 构建数值图 y 轴 ----
    if(!is.null(ymax_counts)) {
        y_limits_counts <- c(0, ymax_counts)
        if(!is.null(yby_counts)) {
            y_breaks_counts <- seq(0, ymax_counts, yby_counts)
        } else {
            y_breaks_counts <- waiver()
        }
    } else {
        y_limits_counts <- NULL
        y_breaks_counts <- waiver()
    }
    
    # ---- 7. 数值堆叠柱状图 ----
    p_counts <- ggplot(data_long, aes(x = Group, y = Count, fill = Category)) +
                geom_col(position = "stack", width = 0.8) +
                geom_text(
                    aes(label = ifelse(Count > 0, Count, "")),
                    position = position_stack(vjust = 0.5),
                    color = "white",
                    size = 4,
                    fontface = "bold"
                ) +
                geom_text(
                    data = data,
                    aes(x = Group, y = Total + max_total * 0.03, 
                        label = paste0("Total:", Total)),
                    inherit.aes = FALSE,
                    size = 3.5,
                    color = "black"
                ) +
                scale_fill_manual(values = colors_use) +
                labs(
                    title = paste(prefix, "(Stacked Bar - Counts)"),
                    x = "Group",
                    y = "Count",
                    fill = "Category"
                ) +
                theme_bw() +
                theme(
                    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                    axis.text.y = element_text(size = 10),
                    legend.position = "top",
                    panel.grid.major.x = element_blank()
                )
    
    # ---- 7.1 添加数值图 y 轴控制 ----
    if(!is.null(y_limits_counts)) {
        p_counts <- p_counts + scale_y_continuous(
            limits = y_limits_counts,
            breaks = y_breaks_counts,
            expand = expansion(mult = c(0, 0.05))
        )
    } else {
        p_counts <- p_counts + scale_y_continuous(
            expand = expansion(mult = c(0, 0.1))
        )
    }
    
    # ---- 8. 百分比堆叠柱状图（固定 0-100%） ----
    p_percent <- ggplot(data_percent_long, aes(x = Group, y = Percentage, fill = Category)) +
        geom_col(position = "stack", width = 0.6) +
        geom_text(
            aes(label = ifelse(Percentage > 3, paste0(round(Percentage, 1), "%"), "")),
            position = position_stack(vjust = 0.5),
            color = "white",
            size = 4,
            fontface = "bold"
        ) +
        geom_text(
            data = data,
            aes(x = Group, y = 102, label = paste0("n=", Total)),
            inherit.aes = FALSE,
            size = 3.5,
            color = "black"
        ) +
        scale_fill_manual(values = colors_use) +
        labs(
            title = paste(prefix, "(Stacked Bar - Percentage)"),
            x = "Group",
            y = "Percentage (%)",
            fill = "Category"
        ) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 10),
            legend.position = "top",
            panel.grid.major.x = element_blank()
        ) +
        scale_y_continuous(
            limits = c(0, 105),
            expand = expansion(mult = c(0, 0.02))
        )
    
    # ---- 9. 保存图片 ----
    if(!dir.exists(outdir)) {
        dir.create(outdir, recursive = TRUE)
    }
    
    outfile_counts <- file.path(outdir, paste0(prefix, "_counts.pdf"))
    outfile_percent <- file.path(outdir, paste0(prefix, "_percentage.pdf"))
    
    ggsave(outfile_counts, p_counts, width = 5, height = 6)
    ggsave(outfile_percent, p_percent, width = 5, height = 6)
    
    # ---- 10. 输出信息 ----
    cat("\n===== 完成 =====\n")
    cat("数值堆叠柱状图:", outfile_counts, "\n")
    cat("百分比堆叠柱状图:", outfile_percent, "\n")
    cat("\n=== Y轴设置 ===\n")
    cat("数值图 - 最大值:", ifelse(is.null(ymax_counts), "auto", ymax_counts), "\n")
    cat("数值图 - 间隔:", ifelse(is.null(yby_counts), "auto", yby_counts), "\n")
    cat("百分比图 - 固定: 0-105\n")
}

# ============================================================
# 使用示例
# ============================================================

# 示例1: 默认自动
# plot_stacked_bar(
#     infile = "diff_genes.csv",
#     group_col = "Group",
#     value_cols = "A3SS,A5SS,MXE,RI,SE",
#     outdir = "./results/",
#     prefix = "AS_events"
# )

# 示例2: 固定数值图 y 轴
# plot_stacked_bar(
#     infile = "diff_genes.csv",
#     group_col = "Group",
#     value_cols = "A3SS,A5SS,MXE,RI,SE",
#     outdir = "./results/",
#     prefix = "AS_events",
#     ymax_counts = 300,
#     yby_counts = 50
# )

# 示例3: 只控制最大值，间隔自动
# plot_stacked_bar(
#     infile = "diff_genes.csv",
#     group_col = "Group",
#     value_cols = "A3SS,A5SS,MXE,RI,SE",
#     outdir = "./results/",
#     prefix = "AS_events",
#     ymax_counts = 250
# )