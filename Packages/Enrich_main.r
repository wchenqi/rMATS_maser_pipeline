#### conda activate Seurat
## v4版本更新：
# 调整对GO富基结果的similarity可视化部分
# 参考链接: https://zhuanlan.zhihu.com/p/659582431

#keytypes(org.Mm.eg.db)
#[1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
#[6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"   
#[11] "GENETYPE"     "GO"           "GOALL"        "IPI"          "MGI"        
#[16] "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"       
#[21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UNIPROT"

#cols <- c("#E4191C","#66C2A5","#1E78B5","#B2DE69","#006026","#FA9B99","#6A3D9A","#D76100","#80B1D3","#FFD930","#D34062",
#            "#492900","#FAD6BC","#17F227","#C0AC53","#993F5D","#D995CC","#0FE3FA","#D87695","#E8E419","#0D48F5","#D603F5","#5574DE","#5D525D",
#            "#011840","#EA5A5D","#670D5B","#265A7A","#B2B2FC","#F4CDE1","#ED766C","#0273B4","#FF8600","#C60A74","#1E8F87","#81C0B8","#2E9736",
#            "#F4AF5E","#623A8A","#C4B2CA","#A2CB86","#7CB1D8","#9E0242","#FEAF61","#68C3A6","#4E9222","#2266AC","#552789","#E31A1A","#FDB964",
#            "#A7D96B","#4576B5","#820E7D")

### 脚本需要的加载包
# suppressPackageStartupMessages({
#         library(dplyr)
#         library(tidyr)
#         library(proxy)
#         library(igraph)
#         library(pheatmap)
#         library(ggplot2)
# })


# source("/data/med-wangcq/01CondaEnv/02Git_repo/00MyGit_wchenqi/Multi_omics/SplicingAnalysis/00Packages/01bulk/01_LoadPackages_v2.r")
# source("/data/med-wangcq/01CondaEnv/02Git_repo/00MyGit_wchenqi/Multi_omics/SplicingAnalysis/00Packages/ColorPanel.r")

## 整合的富集分析和可视化函数
## 功能：执行富集分析、识别关键通路和核心基因、生成多种可视化图表
analyze_enrichment <- function(species, genes_use, tp, outdir1 = "./", 
                               spEnrch = FALSE, 
                               pval_cutoff = 0.05, 
                               count_cutoff = 3,
                               top_pathways = 20,
                               top_genes = 20,
                               min_genes = 5) {
    
    # 检查输入
    if (length(genes_use) < min_genes) {
        warning(paste("基因数", length(genes_use), "小于", min_genes, "，跳过富集分析"))
        return(NULL)
    }
    
    # 创建输出目录
    outdir_enrich <- file.path(outdir1, "EnrichmentAnalysis")
    if(!dir.exists(outdir_enrich)) dir.create(outdir_enrich, recursive = TRUE)
    
    # ========== 1. 执行富集分析 ==========
    cat("\n========== 1. 执行富集分析 ==========\n")
    
    # 调用已有富集函数，添加错误处理
    Enrich_result <- tryCatch({
        go_kegg(species = species, genes_use = genes_use, tp, outdir1 = outdir1, spEnrch = spEnrch)
    }, error = function(e) {
        warning("富集分析执行失败: ", e$message)
        return(NULL)
    })
    
    # 检查富集结果是否有效
    if (is.null(Enrich_result)) {
        warning("富集分析返回空结果，跳过")
        return(NULL)
    }
    
    # 提取 GO 和 Reactome 结果（根据 go_kegg 返回的格式调整）
    if (is.list(Enrich_result) && length(Enrich_result) >= 2) {
        enrich.go <- Enrich_result[[1]]
        enrich.reactome <- Enrich_result[[2]]
    } else {
        warning("富集结果格式不正确，跳过")
        return(NULL)
    }
    
    # ========== 2. 筛选显著通路 ==========
    cat("\n========== 2. 筛选显著通路 ==========\n")
    
    # GO 筛选
    go_sig <- data.frame()
    if (!is.null(enrich.go) && !is(enrich.go, "try-error")) {
        go_sig <- tryCatch({
            df <- as.data.frame(enrich.go)
            if(nrow(df) > 0) {
                df <- df %>%
                    filter(p.adjust < pval_cutoff, Count >= count_cutoff) %>%
                    arrange(p.adjust)
            }
            df
        }, error = function(e) data.frame())
    }
    
    # Reactome 筛选
    reactome_sig <- data.frame()
    if (!is.null(enrich.reactome) && !is(enrich.reactome, "try-error")) {
        reactome_sig <- tryCatch({
            df <- as.data.frame(enrich.reactome)
            if(nrow(df) > 0) {
                df <- df %>%
                    filter(p.adjust < pval_cutoff, Count >= count_cutoff) %>%
                    arrange(p.adjust)
            }
            df
        }, error = function(e) data.frame())
    }
    
    cat("显著 GO 通路数:", nrow(go_sig), "\n")
    cat("显著 Reactome 通路数:", nrow(reactome_sig), "\n")
    
    # 如果没有显著结果，跳过后续分析
    if (nrow(go_sig) == 0 && nrow(reactome_sig) == 0) {
        warning("没有显著的富集通路，跳过核心通路/基因识别")
        
        # 创建空结果文件说明
        empty_note <- data.frame(
            Message = "No significant pathways found",
            Gene_Count = length(genes_use),
            FDR_Threshold = pval_cutoff,
            Count_Threshold = count_cutoff
        )
        write.csv(empty_note, file.path(outdir_enrich, paste0(tp, "_NO_SIGNIFICANT_RESULTS.csv")), row.names = FALSE)
    }
    
    # ========== 3. 识别关键通路（使用多种方法） ==========
    cat("\n========== 3. 识别关键通路 ==========\n")
    
    core_pathways <- character()
    pathway_sim <- NULL
    g <- NULL
    gene_pathway_mat <- NULL
    
    # 方法1：基于 Reactome 的网络分析（优先）
    if (nrow(reactome_sig) >= 3) {
        cat("  使用方法1: Reactome 网络分析\n")
        tryCatch({
            # 构建基因-通路矩阵
            gene_pathway_mat <- enrich.reactome@result %>%
                filter(ID %in% reactome_sig$ID) %>%
                dplyr::select(ID, geneID) %>%
                tidyr::separate_rows(geneID, sep = "/") %>%
                mutate(value = 1) %>%
                tidyr::pivot_wider(names_from = ID, values_from = value, values_fill = 0)
            
            if (nrow(gene_pathway_mat) > 1 && ncol(gene_pathway_mat) > 1) {
                gene_pathway_mat <- as.data.frame(gene_pathway_mat)
                rownames(gene_pathway_mat) <- gene_pathway_mat$geneID
                gene_pathway_mat$geneID <- NULL
                
                # 计算相似度
                pathway_sim <- as.matrix(proxy::simil(t(gene_pathway_mat), method = "Jaccard"))
                pathway_sim[is.na(pathway_sim)] <- 0
                diag(pathway_sim) <- 0
                
                # 构建网络
                g <- graph_from_adjacency_matrix(pathway_sim, 
                                                 mode = "undirected", 
                                                 weighted = TRUE, 
                                                 diag = FALSE)
                
                # 计算中心性
                V(g)$degree <- igraph::degree(g)
                V(g)$betweenness <- betweenness(g, weights = NA)
                V(g)$closeness <- closeness(g, weights = NA)
                
                # 综合评分（标准化的加权和）
                degree_norm <- scale(V(g)$degree)
                betweenness_norm <- scale(V(g)$betweenness)
                closeness_norm <- scale(V(g)$closeness)
                V(g)$core_score <- degree_norm + betweenness_norm + closeness_norm
                
                # 选择核心通路（综合评分 Top）
                core_pathways <- V(g)$name[order(V(g)$core_score, decreasing = TRUE)]
                core_pathways <- core_pathways[1:min(top_pathways, length(core_pathways))]
                
                cat("  网络分析识别到", length(core_pathways), "个核心通路\n")
            }
        }, error = function(e) {
            cat("  网络分析失败:", e$message, "\n")
        })
    }
    
    # 方法2：基于 p.adjust 排序（备选）
    if (length(core_pathways) == 0 && nrow(reactome_sig) > 0) {
        cat("  使用方法2: 按显著性排序\n")
        core_pathways <- reactome_sig$ID[1:min(top_pathways, nrow(reactome_sig))]
    }
    
    # 方法3：基于 Count 排序（最后的备选）
    if (length(core_pathways) == 0 && nrow(go_sig) > 0) {
        cat("  使用方法3: 按基因数排序\n")
        core_pathways <- go_sig$ID[1:min(top_pathways, nrow(go_sig))]
    }
    
    # ========== 4. 识别核心基因 ==========
    cat("\n========== 4. 识别核心基因 ==========\n")
    
    core_genes <- character()
    
    # 方法1：基于核心通路的基因-通路矩阵
    if (length(core_pathways) > 0 && exists("gene_pathway_mat") && !is.null(gene_pathway_mat) && ncol(gene_pathway_mat) > 0) {
        cat("  使用方法1: 基于核心通路矩阵\n")
        tryCatch({
            core_pathways_exist <- core_pathways[core_pathways %in% colnames(gene_pathway_mat)]
            if (length(core_pathways_exist) > 0) {
                core_gene_pathway_mat <- gene_pathway_mat[, core_pathways_exist, drop = FALSE]
                gene_pathway_count <- rowSums(core_gene_pathway_mat)
                
                # 至少出现在2个核心通路中
                core_genes <- names(gene_pathway_count[gene_pathway_count >= 2])
                core_genes <- sort(gene_pathway_count[core_genes], decreasing = TRUE)
                cat("  识别到", length(core_genes), "个核心基因\n")
            }
        }, error = function(e) {
            cat("  方法1失败:", e$message, "\n")
        })
    }
    
    # 方法2：基于显著通路中出现的频率
    if (length(core_genes) == 0 && nrow(reactome_sig) > 0) {
        cat("  使用方法2: 基于显著通路频率\n")
        tryCatch({
            all_genes <- strsplit(reactome_sig$geneID, "/") %>% unlist()
            gene_freq <- table(all_genes)
            core_genes <- sort(gene_freq[gene_freq >= 2], decreasing = TRUE)
            cat("  识别到", length(core_genes), "个核心基因\n")
        }, error = function(e) {
            cat("  方法2失败:", e$message, "\n")
        })
    }
    
    # 方法3：基于 GO 通路频率
    if (length(core_genes) == 0 && nrow(go_sig) > 0) {
        cat("  使用方法3: 基于 GO 通路频率\n")
        tryCatch({
            all_genes <- strsplit(go_sig$geneID, "/") %>% unlist()
            gene_freq <- table(all_genes)
            core_genes <- sort(gene_freq[gene_freq >= 2], decreasing = TRUE)
            cat("  识别到", length(core_genes), "个核心基因\n")
        }, error = function(e) {
            cat("  方法3失败:", e$message, "\n")
        })
    }
    
    # Top 核心基因
    top_core_genes <- names(head(core_genes, top_genes))
    
    # ========== 5. 保存结果（仅当非空时） ==========
    cat("\n========== 5. 保存结果文件 ==========\n")
    
    # 保存 GO 结果（如果非空）
    if (nrow(go_sig) > 0) {
        write.csv(go_sig, file.path(outdir_enrich, paste0(tp, "_GO_Significant.csv")), row.names = FALSE)
        cat("  已保存 GO 显著通路\n")
    }
    
    # 保存 Reactome 结果（如果非空）
    if (nrow(reactome_sig) > 0) {
        write.csv(reactome_sig, file.path(outdir_enrich, paste0(tp, "_Reactome_Significant.csv")), row.names = FALSE)
        cat("  已保存 Reactome 显著通路\n")
    }
    
    # 保存核心通路（如果非空）
    if (length(core_pathways) > 0) {
        core_pathways_df <- data.frame(
            Pathway = core_pathways, 
            Rank = 1:length(core_pathways),
            Source = ifelse(core_pathways %in% reactome_sig$ID, "Reactome", "GO")
        )
        write.csv(core_pathways_df,
                  file.path(outdir_enrich, paste0(tp, "_Core_Pathways.csv")), 
                  row.names = FALSE)
        cat("  已保存", length(core_pathways), "个核心通路\n")
    }
    
    # 保存核心基因（如果非空）
    if (length(core_genes) > 0) {
        core_genes_df <- data.frame(
            Gene = names(core_genes), 
            Pathway_Count = as.numeric(core_genes),
            Rank = 1:length(core_genes)
        )
        write.csv(core_genes_df,
                  file.path(outdir_enrich, paste0(tp, "_Core_Genes.csv")),
                  row.names = FALSE)
        cat("  已保存", length(core_genes), "个核心基因\n")
    }
    
    # 保存通路网络图（如果存在）
    if (!is.null(g) && igraph::vcount(g) > 1) {
        pdf(file.path(outdir_enrich, paste0(tp, "_Pathway_Network.pdf")), width = 12, height = 10)
        plot(g, 
             vertex.size = pmin(pmax(V(g)$degree * 2, 5), 20),
             vertex.label.cex = 0.7,
             vertex.label.color = "black",
             edge.width = pmax(E(g)$weight * 2, 0.5),
             main = paste(tp, "- Pathway Similarity Network"))
        dev.off()
        cat("  已保存通路网络图\n")
    }
    
    # 保存通路相似度热图
    if (!is.null(pathway_sim) && ncol(pathway_sim) > 1 && ncol(pathway_sim) <= 50) {
        pdf(file.path(outdir_enrich, paste0(tp, "_Pathway_Heatmap.pdf")), width = 10, height = 8)
        pheatmap::pheatmap(pathway_sim,
                           main = paste(tp, "- Pathway Similarity Heatmap"),
                           fontsize_row = 8,
                           fontsize_col = 8)
        dev.off()
        cat("  已保存通路热图\n")
    }
    
    # 保存核心基因条形图
    if (length(core_genes) > 0) {
        top_genes_df <- data.frame(
            Gene = names(head(core_genes, min(top_genes, length(core_genes)))),
            Count = as.numeric(head(core_genes, min(top_genes, length(core_genes))))
        )
        
        p_bar <- ggplot(top_genes_df, aes(x = reorder(Gene, Count), y = Count, fill = Count)) +
            geom_col(show.legend = FALSE) +
            coord_flip() +
            scale_fill_gradient(low = "lightblue", high = "darkred") +
            theme_bw() +
            labs(title = paste(tp, "- Core Genes"),
                 x = "Gene", y = "Number of Associated Pathways")
        
        pdf(file.path(outdir_enrich, paste0(tp, "_CoreGenes_Barplot.pdf")), width = 8, height = max(5, nrow(top_genes_df)/2))
        print(p_bar)
        dev.off()
        cat("  已保存核心基因条形图\n")
    }
    
    # 保存 R 对象
    saveRDS(list(
        go = enrich.go,
        reactome = enrich.reactome,
        go_sig = go_sig,
        reactome_sig = reactome_sig,
        core_pathways = core_pathways,
        core_genes = core_genes,
        pathway_network = g,
        gene_pathway_matrix = gene_pathway_mat
    ), file = file.path(outdir_enrich, paste0(tp, "_Enrichment_Results.rds")))
    cat("  已保存 R 对象\n")
    
    # ========== 6. 返回结果 ==========
    cat("\n========== 分析完成！ ==========\n")
    cat("结果保存目录:", outdir_enrich, "\n")
    cat("核心通路数:", length(core_pathways), "\n")
    cat("核心基因数:", length(core_genes), "\n")
    
    return(list(
        go = enrich.go,
        reactome = enrich.reactome,
        go_sig = go_sig,
        reactome_sig = reactome_sig,
        core_pathways = core_pathways,
        core_genes = core_genes,
        top_core_genes = top_core_genes
    ))
}

# ========== 使用示例 ==========
# 
# result <- analyze_enrichment(
#     species = "mouse",
#     genes_use = c("Cacna1c", "Adgrg6", "Mybpc3"),
#     tp = "up",
#     outdir1 = "./Enrichment",
#     spEnrch = FALSE,
#     pval_cutoff = 0.05,
#     count_cutoff = 3,
#     top_pathways = 20,
#     top_genes = 20,
#     min_genes = 5
# )
# 
# # 查看结果
# print(result$core_pathways)  # 核心通路
# print(result$core_genes)      # 核心基因

# 生成的输出文件说明
# {tp}_GO_Significant.csv        - 显著 GO 通路完整表格
# {tp}_Reactome_Significant.csv  - 显著 Reactome 通路完整表格
# {tp}_Core_Pathways.csv         - 识别出的核心通路列表
# {tp}_Core_Genes.csv            - 识别出的核心基因列表
# {tp}_Pathway_Network.pdf       - 通路相似度网络图
# {tp}_Pathway_Heatmap.pdf       - 通路相似度热图
# {tp}_CoreGenes_Barplot.pdf     - 核心基因条形图
# {tp}_Enrichment_Results.rds    - 完整的 R 对象