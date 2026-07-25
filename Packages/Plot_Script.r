
##### 画图脚本 =============================================
    ##1) 绘制火山图 [单剪接事件处理]
    #'
    #' @param summary_table 数据框，事件汇总表（需包含EvenType, IncLevelDifference, negLogFDR, geneSymbol, ID列）
    #' @param label 字符，标注方式："topN"（N为数字）标注前N个上调和下调事件，
    #'              或逗号分隔的基因列表标注特定基因
    #' @param evt 字符，事件类型（"SE","A3SS","A5SS","RI","MXE"）
    #' @param outdir 字符，输出目录路径
    #' @param fdr 数值，FDR阈值，默认0.05
    #' @param psi 数值，PSI阈值，默认0.1
    #'
    #' @return 无返回值，火山图保存为PDF文件
    #' @export
    volcano_plot <- function(summary_table,label="top5",evt=c("SE","RI","A3SS","A5SS","MXE"),outdir="./",fdr=0.05,psi=0.1){
            # summary_table: 所有差异剪接分析输出结果
            # label: 指定top N或者eventLabel格式
            ## 建立文件夹
            if(!file.exists(outdir)){
                dir.create(outdir,recursive=TRUE)
            }
            print(outdir)
            data <- summary_table %>% mutate(UpDown = ifelse(FDR <= fdr & abs(IncLevelDifference) >= psi,
                                                      ifelse(IncLevelDifference >= psi, "Up", "Down"), "Stable"),
                                             Rank = IncLevelDifference/(-log10(FDR)))       # 越显著分母越大,Rank顺序排列
            data_sig <- data[which(data$UpDown %in% c("Up","Down")),]

            ## 加标签
            if (grepl("^top", label)) {
                # 如果展示topN
                data_sig <- data_sig[which(data_sig$EvenType %in% evt),]
                data_sig <- data_sig[order(data_sig$Rank, decreasing=FALSE),]       # 显著且高倍
                top_n <- as.numeric(gsub("top", "", label))
                up_lab <- data_sig[data_sig$UpDown == "Up",] %>% head(top_n)
                down_lab <- data_sig[data_sig$UpDown == "Down",] %>% tail(top_n)
                top_label <- rbind(up_lab, down_lab) %>% mutate(lab = paste0(geneSymbol, "_", ID))
                plot_suffix <- label
            } else {
                # 如果展示指定基因列表
                genes_vec <- unlist(strsplit(label, ","))
                # 和显著差异结果取交集
                evt_sig <- data_sig[which(data_sig$geneSymbol %in% genes_vec),]
                data$lab <- ifelse(data$eventLabel %in% evt_sig$eventLabel, data$eventLabel, "")
                top_label <- data[data$lab != "", ]
                plot_suffix <- "custom"
            }
            # write.table(top_label,file.path(outdir, paste0("Volcano_", plot_suffix, ".tsv")),sep="\t",row.names=FALSE,quote=FALSE)
            wb <- createWorkbook()
            addWorksheet(wb, "Volcano_label")
            addWorksheet(wb, "Volcano_all")
            writeData(wb, "Volcano_label", top_label)
            writeData(wb, "Volcano_all", data)
            saveWorkbook(wb, file.path(outdir, paste0("Volcano_", plot_suffix, ".xlsx")), overwrite = TRUE)

            ## 画图
            if (nrow(top_label) > 0) {
                volcano_p <- ggplot(data, aes(x = IncLevelDifference, y = negLogFDR, color = UpDown)) +
                                    geom_point(alpha = 0.5, size = 1.5) +
                                    scale_color_manual(values = c("Stable" = "grey70", "Up" = "#D62728", "Down" = "#272ad6")) +
                                    geom_hline(yintercept = -log10(fdr), linetype = "dashed", color = "grey40") +
                                    geom_vline(xintercept = c(-psi, psi), linetype = "dashed", color = "grey40") +
                                    geom_text_repel(data = top_label,
                                                    aes(x = IncLevelDifference, y = negLogFDR, label = lab),
                                                    color = "black", size = 2,
                                                    inherit.aes = FALSE, max.overlaps = 15,
                                                    box.padding = 0.5, min.segment.length = 0.1) + 
                                    labs(x = paste0(evt,": Delta PSI ",psi), y = paste0("-log10(FDR) ",fdr)) +
                                    theme_bw() + theme(legend.position = "none")

                pdf(file.path(outdir, paste0("Volcano_", plot_suffix, ".pdf")), width = 4, height = 5)
                print(volcano_p)
                dev.off()
                message("    火山图已保存 (", nrow(top_label), " 个标注)")
            }
        }
    
    ##2) PCA和PSI箱线图绘制
    #'
    #' @param maser_sig MASER对象，显著差异剪接事件对象
    #' @param evt 字符，事件类型
    #' @param outdir 字符，输出目录路径
    #'
    #' @return 无返回值，PCA图和箱线图保存为PDF文件
    #' @export
    PCA_PSI_plot <- function(maser_sig,evt,outdir){
            pca_result <- tryCatch(maser::pca(maser_sig, type = evt), error = function(e) NULL)
            pca_dim <- pca_result@data

            p <- ggplot(pca_dim, aes(PC1, PC2, color = Condition, label = Samples)) +
                        geom_point(size = 5) + 
                        geom_text_repel(max.overlaps = 20, size = 3) +  # 自动避免重叠
                        scale_colour_manual(values = c("blue","red")) + 
                        theme_bw() + 
                        theme(axis.text.x = element_text(size = 12),
                            axis.text.y = element_text(size = 12), 
                            axis.title.x = element_text(face = "plain", colour = "black", size = 12), 
                            axis.title.y = element_text(face = "plain", colour = "black", size = 12), 
                            panel.grid.minor = element_blank(), 
                            panel.grid.major = element_blank(), 
                            plot.background = element_blank(),
                            legend.title = element_blank(), 
                            legend.position = "right") +
                        labs(title = paste0("PCA - ", evt)) +
                        xlab(pca_result@labels$x) +
                        ylab(pca_result@labels$y)

            if (!is.null(pca_result)) {
                pdf(file.path(outdir, paste0("PCA_", evt, ".pdf")), width = 6, height = 6)
                print(p)
                dev.off()
            }

            # PSI箱线图
            boxplot_p <- tryCatch(boxplot_PSI_levels(maser_sig, type = evt), error = function(e) NULL)
            if (!is.null(boxplot_p)) {
                pdf(file.path(outdir, paste0("Boxplot_", evt, ".pdf")), width = 6, height = 6)
                print(boxplot_p)
                dev.off()
            }
        }

    ##3) PSI热图绘制 [逐个剪接类别输出结果] ============================================================================
    # 独立运行函数, 未被任何函数调用
    # 测试参数
        # maser_obj <- qread("/scratch/2026-06-08/med-wangcq/Others/AnJQ/AnalysisData/Sam68KOWT_TAC_Day3_4W_RNAseq/06rMATS_filter/Sham_4W/00rds/maser_sig.qs")
        # event_type <- "SE,A5SS,A3SS,MXE,RI"
        # gene_list <- "Fggy,Agtpbp1,Ifi27,Surf1,Myh6,Mybpc3"
        # event_list <- "no"
        # prefix <- "DSG_only_intersect"
        # outdir <- "/scratch/2026-06-08/med-wangcq/Others/AnJQ/AnalysisData/Sam68KOWT_TAC_Day3_4W_RNAseq/06rMATS_filter/Sham_4W/"
    # 使用
        # PSI_heatmap(maser_obj, event_type = event_type, 
        #          gene_list = gene_list, 
        #          event_list = event_list, 
        #          prefix = prefix,
        #          outdir = outdir)
    #'
    #' 为指定基因或事件绘制PSI热图，支持自定义事件顺序。
    #'
    #' @param maser_obj MASER对象（可以是过滤前或过滤后的）
    #' @param event_type 字符，逗号分隔的事件类型，如"SE,A5SS,A3SS"
    #' @param gene_list 字符，逗号分隔的基因列表，默认"no"（不筛选）
    #' @param event_list 字符，逗号分隔的事件ID列表，格式为"Evt:geneSymbol_ID"，默认"no"
    #' @param prefix 字符，输出文件名前缀，默认"no"
    #' @param outdir 字符，输出目录路径
    #' @param cluster_rows 逻辑，是否对行聚类，默认TRUE
    #'
    #' @return 无返回值，热图和PSI数据CSV保存到SPI_plot子目录
    #'
    #' @examples
    #' \dontrun{
    #' PSI_heatmap(maser_obj, 
    #'      event_type = "SE,A5SS,A3SS,MXE,RI",
    #'      gene_list = "Fggy,Agtpbp1,Ifi27",
    #'      outdir = "./results/")
    #' }
    #' @export
    PSI_heatmap_boxplot <- function(summary_table, 
                            gene_list = "no", 
                            event_list = "no", 
                            SigOnly = TRUE,
                            gp_name = "default",
                            outdir = "./",
                            cluster_rows = TRUE,
                            show_significance = FALSE) {
            # 创建输出目录
            outdir_plot <- file.path(outdir,"SPI_plot")
            message("Outdir: ",outdir_plot)
            if (!dir.exists(outdir_plot)) {
                dir.create(outdir_plot, recursive = TRUE)
            }
            if (SigOnly){
                sum_plot <- summary_table[which(summary_table$UpDown != "Stable"),]
            }else{
                sum_plot <- summary_table
            }
            rownames(sum_plot) <- sum_plot$eventLabel

            # 根据基因列表筛选
            if (gene_list != "no") {
                gg <- unlist(strsplit(gene_list, ","))
                gg_use <- base::intersect(gg, sum_plot$geneSymbol)
                # 查看是否有差集基因
                diff_gg <- base::setdiff(gg, sum_plot$geneSymbol)
                if (length(diff_gg) > 0){
                    cat(diff_gg," Not Found !!\n")
                }
                if (length(gg_use) > 0) {
                    print(gg_use)
                    plt_label <- sum_plot$eventLabel[which(sum_plot$geneSymbol %in% gg_use)]
                }
            }
                
            # 根据事件ID列表筛选（仅在未提供 gene_list 或 gene_list 无匹配时）
            if (event_list != "no" && gene_list == "no") {
                evt <- unlist(strsplit(event_list, ","))
                evt_use <- intersect(evt, sum_plot$eventLabel)              # 注意这里格式为Evt:gg_id
                # 查看是否有差集基因
                diff_evt <- base::setdiff(evt, sum_plot$eventLabel)
                if (length(diff_evt) > 0){
                    cat(diff_evt," Not Found !!\n")
                }
                if (length(evt_use) > 0) {
                    plt_label <- evt_use
                }
            }
            # 汇总gene_list和event_list结果后提取矩阵
            psi_lab <- sum_plot[plt_label, c("PSI_1","PSI_2"),drop = FALSE]
            # 拆分矩阵PSI数值
            g1 <- str_split(psi_lab$PSI_1,",",simplify=TRUE)
            g2 <- str_split(psi_lab$PSI_2,",",simplify=TRUE)
            psi_plot <- cbind(g1, g2)
            psi_plot <- matrix(as.numeric(psi_plot), nrow = nrow(psi_plot))
            rownames(psi_plot) <- plt_label
            # 重命名
            if (gp_name != "default"){
                gp <- unlist(strsplit(gp_name,","))
            }else{
                gp <- c("Gp1_rep","Gp2_rep")
            }    
            gp1 <- paste0(rep(gp[1],n=ncol(g1)),"_",seq(ncol(g1)))
            gp2 <- paste0(rep(gp[2],n=ncol(g2)),"_",seq(ncol(g2)))
            colnames(psi_plot) <- c(gp1,gp2)
            # 写出文件
            write.csv(psi_plot, file.path(outdir_plot,ifelse(SigOnly, "PSI_Heatmap_Sig.csv", "PSI_Heatmap.csv")), row.names = TRUE, quote=TRUE)
            ### 开始绘图
            print("Plot Event: ")
            ## 绘制热图
            if(nrow(psi_plot) > 0) {
                # 画热图
                wt <- dim(psi_plot)[2] * 0.6 + 1
                ht <- dim(psi_plot)[1] * 1 + 2
                pdf(file.path(outdir_plot, ifelse(SigOnly,"PSI_Heatmap_Sig.pdf","PSI_Heatmap.pdf")),width=wt,height=ht)
                print(pheatmap::pheatmap(psi_plot,
                                         cellwidth = 20,
                                         cellheight = 15,
                                         scale = "row",
                                         show_rownames = TRUE,
                                         show_colnames = TRUE,
                                         cluster_rows = cluster_rows,
                                         cluster_cols = FALSE,
                                         clustering_method = "complete",
                                         main = "Event PSI Heatmap",
                                         color = colorRampPalette(c("#1b74ce", "white", "#e22238"))(50),
                                         breaks = seq(-1, 1, length.out = 51)))
                dev.off()
                
                # ========== 箱线图部分 ==========
                # 准备长格式数据
                psi_long <- psi_plot %>%
                            as.data.frame() %>%
                            rownames_to_column("eventLabel") %>%
                            pivot_longer(cols = -eventLabel, names_to = "sample", values_to = "PSI")
                            
                # 提取分组信息（假设列名格式为 "Group_数字"，如 "TAC_D3_CTR_1"）
                psi_long <- psi_long %>%
                            mutate(
                                group = case_when(
                                    grepl(gp[1], sample) ~ gp[1],
                                    grepl(gp[2], sample) ~ gp[2],
                                    TRUE ~ "Other"
                                )
                            ) %>%
                            filter(group != "Other")
                
                # 计算每组每个事件的均值和标准差
                psi_summary <- psi_long %>%
                                group_by(eventLabel, group) %>%
                                summarise(
                                    mean_PSI = mean(PSI, na.rm = TRUE),
                                    sd_PSI = sd(PSI, na.rm = TRUE),
                                    n = n(),
                                    .groups = "drop"
                                )
                
                # 计算统计检验（t检验或Wilcoxon检验）
                stat_results <- psi_long %>%
                                group_by(eventLabel) %>%
                                summarise(
                                    p_value = tryCatch({
                                        if (length(unique(group)) == 2) {
                                            g1_data <- PSI[group == gp[1]]
                                            g2_data <- PSI[group == gp[2]]
                                            if (length(g1_data) >= 3 && length(g2_data) >= 3) {
                                                t.test(g1_data, g2_data, paired = FALSE)$p.value
                                            } else {
                                                wilcox.test(g1_data, g2_data, paired = FALSE)$p.value
                                            }
                                        } else {
                                            NA
                                        }
                                    }, error = function(e) NA),
                                    .groups = "drop"
                                ) %>%
                                mutate(
                                    sig_label = case_when(
                                        p_value < 0.001 ~ "***",
                                        p_value < 0.01 ~ "**",
                                        p_value < 0.05 ~ "*",
                                        p_value >= 0.05 ~ "ns",
                                        is.na(p_value) ~ "NA"
                                    )
                                )
                
                # 合并统计数据
                psi_summary <- psi_summary %>% left_join(stat_results, by = "eventLabel")
                
                # ========== 创建分组柱状图 + 误差棒 ==========
                
                # 基础柱状图
                p_bar <- ggplot(psi_summary, aes(x = eventLabel, y = mean_PSI, fill = group)) +
                            geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.7) +
                            geom_errorbar(aes(ymin = mean_PSI - sd_PSI, ymax = mean_PSI + sd_PSI),
                                        position = position_dodge(0.9), width = 0.2) +
                            labs(title = "PSI Values by Group",
                                x = "Event",
                                y = "Percent Spliced In (PSI)",
                                fill = "Group") +
                            theme_bw() +
                            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                                axis.text.y = element_text(size = 10),
                                legend.position = "top",
                                plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
                            scale_fill_manual(values = c("#1b74ce", "#e22238"))
                
                # 根据 show_significance 参数决定是否添加显著性标记
                if (show_significance) {
                    p_bar <- p_bar +
                        geom_text(data = psi_summary %>% 
                                    group_by(eventLabel) %>% 
                                    summarise(max_mean = max(mean_PSI, na.rm = TRUE),
                                            sig_label = first(sig_label),
                                            .groups = "drop"),
                                aes(x = eventLabel, y = max_mean + 0.1, label = sig_label),
                                inherit.aes = FALSE, size = 5, vjust = 0)
                }
                
                # 保存柱状图
                bar_width <- max(8, nrow(psi_summary) * 0.3 + 2)
                ggsave(file.path(outdir_plot, ifelse(SigOnly, "PSI_Boxplot_Sig.pdf", "PSI_Boxplot.pdf")),
                    plot = p_bar,
                    width = bar_width,
                    height = 5,
                    limitsize = FALSE)
                
                # ========== 分面箱线图 ==========
                p_box <- ggplot(psi_long, aes(x = group, y = PSI, fill = group)) +
                            geom_boxplot(width = 0.6, outlier.shape = 21, outlier.size = 2) +
                            facet_wrap(~eventLabel, scales = "free_y", ncol = min(4, nrow(psi_summary))) +
                            labs(title = "PSI Distribution by Event and Group",
                                x = "Group",
                                y = "PSI",
                                fill = "Group") +
                            theme_bw() +
                            theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
                                strip.text = element_text(size = 10, face = "bold"),
                                legend.position = "top") +
                            scale_fill_manual(values = c("#1b74ce", "#e22238"))
                
                # 根据 show_significance 参数决定是否在分面箱线图上添加显著性
                if (show_significance) {
                    p_box <- p_box +
                        stat_compare_means(aes(label = ..p.signif..),
                                        method = ifelse(length(unique(psi_long$group)) >= 3, 
                                                        "kruskal.test", "wilcox.test"),
                                        label.x = 1.5,
                                        label.y = max(psi_long$PSI, na.rm = TRUE) * 0.9,
                                        size = 4)
                }
                
                ggsave(file.path(outdir_plot, ifelse(SigOnly, "PSI_FacetBoxplot_Sig.pdf", "PSI_FacetBoxplot.pdf")),
                    plot = p_box,
                    width = min(9, nrow(psi_summary) * 0.4 + 3),
                    height = 5,
                    limitsize = FALSE)
                
                # 保存统计检验结果（始终保存，无论是否显示）
                write.csv(stat_results, 
                        file.path(outdir_plot, ifelse(SigOnly, "PSI_Stats_Sig.csv", "PSI_Stats.csv")),
                        row.names = FALSE, quote = TRUE)
                
                message("✅ 箱线图已保存到: ", outdir_plot)
        }
    }

    #### ===================================================

    ##4) Venn图 ============================================
    Venn_plot <- function(sets,outdir,prefix,width = 8,height = 8){
            # 计算每个集合的大小
            set_sizes <- sapply(sets, length)
            message("集合大小：")
            for (name in names(set_sizes)) {
                message("  ", name, ": ", set_sizes[name])
            }

            # 绘制Venn图
            message("正在绘制Venn图...")
            venn_plot <- ggvenn(sets, 
                                fill_color = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"),
                                stroke_size = 0.5,
                                set_name_size = 5,
                                text_size = 4)

            # 保存图片
            outfile <- file.path(outdir, paste0(prefix,"_venn_plot.pdf"))
            ggsave(filename = outfile, plot = venn_plot,  width = width, height = height, dpi = 300)
            message("Venn图已保存至：", outfile)

            # 可选：同时输出集合交集的详细信息（文本文件）
            intersection_info <- file.path(outdir, paste0(prefix,"_venn_intersect.txt"))
            cat("Venn图集合交集统计\n", file = intersection_info)
            cat("生成时间：", as.character(Sys.time()), "\n\n", file = intersection_info, append = TRUE)
            cat("集合列表：\n", file = intersection_info, append = TRUE)
            for (name in names(sets)) {
                cat(sprintf("  %s: %d 个元素\n", name, length(sets[[name]])), file = intersection_info, append = TRUE)
            }

            # 如果有2-3个集合，输出两两/三三交集的具体元素
            n_sets <- length(sets)
            if (n_sets <= 3) {
                cat("\n交集详细信息：\n", file = intersection_info, append = TRUE)
                # 生成所有非空子集的组合（排除空集）
                all_combinations <- list()
                for (k in 1:n_sets) {
                    combos <- combn(names(sets), k, simplify = FALSE)
                    for (combo in combos) {
                        intersect_elements <- Reduce(intersect, sets[combo])
                        if (length(intersect_elements) > 0) {
                            combo_str <- paste(combo, collapse = " ∩ ")
                            cat(sprintf("\n%s (%d个元素):\n", 
                                combo_str, 
                                length(intersect_elements)), 
                                file = intersection_info, append = TRUE)
                            cat(paste(intersect_elements, collapse = ","), "\n", file = intersection_info, append = TRUE)
                        }
                    }
                }
            }
            message("交集统计信息已保存至：", intersection_info)

            message("完成！")
    }
    #### ===================================================

    ##5) 剪接类型 x 分组堆叠柱状图 ===========================
    ## 参考maser::splicingDistribution进行调整
    ## 输入: maser_sig
    ## 测试参数: 
    # basic_re <- qread("/scratch/2026-06-22/med-wangcq/Others/AnJQ/AnalysisData/Sam68KOWT_TAC_Day3_4W_RNAseq/06rMATS_filter/Sham_4W/00rds/maserfilter_basic_result.qs")
    # maser_sig <- basic_re[["maser_sig"]]
    # ======================================================
    ## 这个好看点,和maser原生结果风格一致,但是不标注数值 =======
    EventType2Group <- function(maser_obj, fdr = 0.05, deltaPSI = 0.1, outdir = "./", wt = 7,ht = 7){
                if (!is(maser_obj, "Maser")) {
                    stop("Parameter events has to be a maser object.")
                }
                events <- as(maser_obj, "list")
                as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
                nevents_cond1 <- rep(0, length(as_types))
                nevents_cond2 <- rep(0, length(as_types))

                for (i in 1:length(as_types)) {
                    stats <- events[[paste0(as_types[i], "_", "stats")]]
                    cond1 <- dplyr::filter(stats, FDR < fdr, IncLevelDifference > deltaPSI)         # 上调剪接事件
                    cond2 <- dplyr::filter(stats, FDR < fdr, IncLevelDifference < (-1 * deltaPSI))  # 下调剪接事件
                    nevents_cond1[i] <- length(cond1$ID)                                            # 统计上调差异剪接事件数量
                    nevents_cond2[i] <- length(cond2$ID)                                            # 统计下调差异剪接事件数量
                }
                # 计算占比
                nevents_prop1 <- nevents_cond1/sum(nevents_cond1)
                nevents_prop2 <- nevents_cond2/sum(nevents_cond2)

                cond1_condition <- rep(events$conditions[1], length(as_types))
                cond2_condition <- rep(events$conditions[2], length(as_types))
                # 构建矩阵
                df.plot <- data.frame(Condition = c(cond1_condition,cond2_condition), 
                                      Type = c(as_types, as_types), 
                                      Count = c(nevents_cond1, nevents_cond2),
                                      Prop = c(nevents_prop1, nevents_prop2))
                # 设置factor
                df.plot$Condition <- factor(df.plot$Condition, levels = c(events$conditions[1], events$conditions[2]))
                df.plot$Type <- factor(df.plot$Type, levels = as_types)
                # write.csv(df.plot,file.path(outdir,"EventType2Group_count_prop.csv"),row.names=FALSE,quote=FALSE)

                p <- ggplot(df.plot, aes(x = Type, y = Count, colour = Condition,fill = Condition)) + 
                                    geom_bar(stat = "identity", alpha = 0.6) +
                                    theme_bw() + 
                                    theme(axis.text.y = element_text(size = 12,
                                          angle = 0, hjust = 0.5, face = "plain"), axis.text.x = element_text(size = 12,
                                          angle = 0, hjust = 0.5, face = "plain"), axis.title.x = element_text(face = "plain",
                                          colour = "black", size = 12), axis.title.y = element_text(face = "plain",
                                          colour = "black", size = 12), legend.title = element_blank(),
                                          legend.text = element_text(face = "plain", colour = "black", size = 12), panel.grid = element_blank()) + 
                                    ylab("Count of splicing events") + xlab("") + 
                                    scale_fill_brewer(palette = "Set2") + scale_color_brewer(palette = "Set2") +
                                    coord_flip()
                                    
                pdf(file.path(outdir, "AS_Event_CountSig.pdf"), width = wt, height = ht)
                print(p)
                dev.off()
            }


    ## 这个自己写的,不好看,但是标注数值占比 ===================
    EventType2Group_test <- function(maser_obj, fdr = 0.05, deltaPSI = 0.1, outdir = "./", wt = 7,ht = 7){
                source("/data/med-wangcq/01CondaEnv/02Git_repo/00MyGit_wchenqi/Multi_omics/RNA-seq/02bulkRNA/03Usable_Script/02PlotScript/08Count_Percent_Barplot.r")
                if (!is(maser_obj, "Maser")) {
                    stop("Parameter events has to be a maser object.")
                }
                events <- as(maser_obj, "list")
                as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
                nevents_cond1 <- rep(0, length(as_types))
                nevents_cond2 <- rep(0, length(as_types))
                
                for (i in 1:length(as_types)) {
                    stats <- events[[paste0(as_types[i], "_", "stats")]]
                    cond1 <- dplyr::filter(stats, FDR < fdr, IncLevelDifference > deltaPSI)         # 上调剪接事件
                    cond2 <- dplyr::filter(stats, FDR < fdr, IncLevelDifference < (-1 * deltaPSI))  # 下调剪接事件
                    nevents_cond1[i] <- length(cond1$ID)                                            # 统计上调差异剪接事件数量
                    nevents_cond2[i] <- length(cond2$ID)                                            # 统计下调差异剪接事件数量
                }

                # 构建矩阵
                df.plot <- data.frame(Gourp = as_types, cond1 = nevents_cond1, cond2 = nevents_cond2)
                colnames(df.plot) <- c("Gourp",events$conditions[1],events$conditions[2])
                write.csv(df.plot,file.path(outdir,"EventType2Group_count.csv"),row.names=FALSE,quote=FALSE)

                plot_stacked_bar(infile = file.path(outdir,"EventType2Group_count.csv"), 
                                 group_col = "Group",
                                 value_cols = "default",
                                 colors = "default",
                                 outdir = outdir,
                                 prefix = "EventType2Group_count",
                                 ymax_counts = NULL,      # 数值图 y 轴最大值
                                 yby_counts = NULL)
            }

    #### ===================================================

    ##5) 画PSI箱线图(没改) ==================================
    plotPSI_boxPlot <- function (events, type = c("A3SS", "A5SS", "SE", "RI", "MXE"), show_replicates = TRUE){
            if (!is(events, "Maser")) {
                stop("Parameter events has to be a maser object.")
            }
            type <- match.arg(type)
            events <- as(events, "list")
            annot <- events[[paste0(type, "_", "events")]]
            if (length(unique(annot$geneSymbol)) > 1) {
                stop(cat("Multiple genes found. Use geneEvents() to select AS events."))
            }
            PSI <- events[[paste0(type, "_", "PSI")]]
            PSI_long <- reshape2::melt(PSI)
            colnames(PSI_long) <- c("ID", "Sample", "PSI")
            Condition <- rep("NA", nrow(PSI_long))
            idx.cond1 <- grep(paste0("^", events$conditions[1]), x = PSI_long$Sample,
                perl = TRUE)
            idx.cond2 <- grep(paste0("^", events$conditions[2]), x = PSI_long$Sample,
                perl = TRUE)
            Condition[idx.cond1] <- events$conditions[1]
            Condition[idx.cond2] <- events$conditions[2]
            PSI_long <- cbind(PSI_long, Condition)
            if (show_replicates) {
                ggplot(PSI_long, aes(x = Condition, y = PSI, fill = Condition,
                    color = Condition)) + geom_violin(trim = FALSE, alpha = 0.6) +
                    geom_jitter(position = position_jitter(0.05), size = 2) +
                    theme_bw() + theme(axis.text.x = element_text(size = 12,
                    angle = 45, hjust = 1), axis.text.y = element_text(size = 12),
                    axis.title.x = element_text(face = "plain", colour = "black",
                        size = 12), axis.title.y = element_text(face = "plain",
                        colour = "black", size = 12), legend.title = element_blank(),
                    legend.text = element_text(face = "plain", colour = "black",
                        size = 12)) + ylab(paste(type, "PSI")) + scale_y_continuous(limits = c(-0.1,
                    1.05)) + scale_fill_manual(values = c("blue", "red")) +
                    scale_color_manual(values = c("blue", "red")) + facet_grid(. ~
                    ID)
            }
            else {
                ggplot(PSI_long, aes(x = Condition, y = PSI, fill = Condition,
                    color = Condition)) + geom_violin(trim = FALSE) +
                    stat_summary(fun.y = mean, geom = "point", size = 2,
                        color = "black") + theme_bw() + theme(axis.text.x = element_text(size = 12,
                    angle = 45, hjust = 1), axis.text.y = element_text(size = 12),
                    axis.title.x = element_text(face = "plain", colour = "black",
                        size = 12), axis.title.y = element_text(face = "plain",
                        colour = "black", size = 12), legend.title = element_blank(),
                    legend.text = element_text(face = "plain", colour = "black",
                        size = 12)) + ylab(paste(type, "PSI")) + scale_y_continuous(limits = c(-0.1,
                    1.05)) + scale_fill_manual(values = c("blue", "red")) +
                    scale_color_manual(values = c("blue", "red")) + facet_grid(. ~
                    ID)
            }
        }
