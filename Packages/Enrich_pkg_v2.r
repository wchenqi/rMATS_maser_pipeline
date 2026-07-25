#!/usr/bin/env Rscript

##### 运行环境: mamba activate r452
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

### 脚本使用需要的包
# suppressPackageStartupMessages({
#         library(clusterProfiler)
#         library(org.Hs.eg.db)
#         library(org.Mm.eg.db)
#         library(org.Rn.eg.db)
#         library(ReactomePA)
#         library(dplyr)
#         library(ggplot2)
#         library(stringr)
#         library(tidyr)
# })

## 参数说明 ===================
    # species         物种, mouse, human, rat
    # genes_use       基因向量
    # tp              输出文件前缀
    # outdir1         指定输出路径
    # spEnrch         是否简化富集分析结果, 默认 FALSE
## ===========================

## 本函数对每个基因单独功能注释, 在部分调用
# ========= 基于数据库对单基因涉及通路注释 ======================
gene_Anno <- function(genes, species = "mouse", anno_col = "default") {
    cat(" ---- 对基因进行功能注释 ---- \n")
    # 支持单个和多个基因提取
    ## 对基因进行功能注释, 单独存储为list结构,包含: gene_func_anno, GO_summary, KEGG_summary
    # columns(org.Mm.eg.db)     # 查询可用列名
    # [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
    # [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"
    # [11] "GENETYPE"     "GO"           "GOALL"        "IPI"          "MGI"
    # [16] "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"
    # [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UNIPROT"
    # columns(org.Hs.eg.db)
    # [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
    # [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"
    # [11] "GENETYPE"     "GO"           "GOALL"        "IPI"          "MAP"
    # [16] "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"
    # [21] "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"
    # [26] "UNIPROT"

    # 去重，避免重复查询
    genes <- unique(genes)
    
    org <- switch(species,
                  mouse = org.Mm.eg.db,
                  human = org.Hs.eg.db,
                  stop("species must be 'mouse' or 'human'"))
    
    col_all <- columns(org)
    col <- if (length(anno_col) > 1) {
        anno_col
    } else if (anno_col == "default") {
        c("SYMBOL", "ENSEMBL", "ENTREZID", "GENETYPE", "GENENAME", 
          "UNIPROT", "PFAM", "PROSITE", "ENZYME", "GO", "ONTOLOGY", "PATH")
    } else if (anno_col == "all") {
        col_all
    } else {
        stop('anno_col must be "default", "all", or a character vector')
    }
    
    col_invalid <- setdiff(col, col_all)
    if (length(col_invalid) > 0) {
        col <- intersect(col, col_all)
        message("Invalid columns removed: ", paste(col_invalid, collapse = ", "))
    }
    
    # mapIds 查询
    anno_list <- lapply(col, function(c) {
        mapIds(org, keys = genes, column = c, keytype = "SYMBOL", multiVals = "list")
    })
    names(anno_list) <- col
    
    # 按基因重组
    result <- lapply(genes, function(g) {
        gene_info <- lapply(anno_list, function(x) {
            v <- x[[g]]
            if (is.null(v) || all(is.na(v))) return(NA)
            if (length(v) == 1) return(v[[1]])
            return(unlist(v))
        })
        names(gene_info) <- names(anno_list)
        gene_info
    })
    names(result) <- genes
    
    # 过滤 GO/ONTOLOGY/PATH 全为 NA 的基因
    result <- result[sapply(result, function(x) {
        has_valid <- function(col_name) {
            if (!col_name %in% names(x)) return(FALSE)
            v <- x[[col_name]]
            !all(is.na(v))
        }
        has_valid("GO") || has_valid("ONTOLOGY") || has_valid("PATH")
    })]
    
    genes_kept <- names(result)
    message("  有效注释基因: ", length(genes_kept), "/", length(genes))
    
    # ========== GO 汇总（用 Map 保持配对，避免 stack 问题） ==========
    go_summary <- if ("GO" %in% col && length(result) > 0) {
        # 用 Map 并行迭代，保持 gene-GO-ONTOLOGY 的原始配对
        go_df <- Map(function(g, go, ont) {
            if (all(is.na(go))) return(NULL)
            # 确保 GO 和 ONTOLOGY 长度一致
            n_go <- length(go)
            n_ont <- if (is.null(ont) || all(is.na(ont))) 0 else length(ont)
            
            # 如果 ONTOLOGY 缺失或长度不匹配，用 NA 填充
            if (n_ont == 0) {
                ont_vec <- rep(NA, n_go)
            } else if (n_ont != n_go) {
                # 长度不匹配时，取第一个或 NA
                ont_vec <- rep(NA, n_go)
                ont_vec[1:min(n_go, n_ont)] <- ont[1:min(n_go, n_ont)]
            } else {
                ont_vec <- ont
            }
            
            data.frame(
                gene = g,
                GO = go,
                ONTOLOGY = ont_vec,
                stringsAsFactors = FALSE
            )
        }, 
        names(result),
        lapply(result, function(x) x$GO),
        if ("ONTOLOGY" %in% col) lapply(result, function(x) x$ONTOLOGY) else lapply(result, function(x) NA)
        )
        
        go_long <- do.call(rbind, go_df)
        rownames(go_long) <- NULL
        
        # 按 GO 汇总
        go_summary <- go_long %>%
            dplyr::filter(!is.na(GO)) %>%
            group_by(GO, ONTOLOGY) %>%
            summarise(
                count = n(),
                genes = paste(unique(gene), collapse = ";"),
                .groups = "drop"
            ) %>%
            arrange(desc(count))
        
        go_summary
    } else NULL
    
    # ========== KEGG 汇总（同样用 Map） ==========
    kegg_summary <- if ("PATH" %in% col && length(result) > 0) {
        path_df <- Map(function(g, path) {
            if (all(is.na(path))) return(NULL)
            data.frame(
                gene = g,
                PATH = path,
                stringsAsFactors = FALSE
            )
        }, 
        names(result),
        lapply(result, function(x) x$PATH)
        )
        
        path_long <- do.call(rbind, path_df)
        rownames(path_long) <- NULL
        
        path_long %>%
            dplyr::filter(!is.na(PATH)) %>%
            group_by(PATH) %>%
            summarise(
                count = n(),
                genes = paste(unique(gene), collapse = ";"),
                .groups = "drop"
            ) %>%
            arrange(desc(count))
    } else NULL
    
    attr(result, "GO_summary") <- go_summary
    attr(result, "KEGG_summary") <- kegg_summary
    attr(result, "genes") <- genes_kept
    attr(result, "columns") <- col
    attr(result, "filtered") <- length(genes) - length(genes_kept)
    class(result) <- c("geneAnnoList", class(result))
    return(result)
}

## 辅助函数,提取gene_Anno输出结果转换矩阵写出文件 =================
gene_Anno_to_df <- function(fun_anno,outdir = "./",prefix){
        data_flat <- lapply(fun_anno,function(x) { 
                        sapply(x, function(y) {
                            if (all(is.na(y))) return(NA_character_)
                                paste0(y,collapse=",")
                            })
                        })
        data_df <- do.call(rbind, data_flat) |> as.data.frame(stringsAsFactors = FALSE)
        data_tb <- data.frame(gene = rownames(data_df), data_df)
        print(head(data_tb))
        # write.table(data_tb, file.path(outdir,paste0(prefix,"_GeneAnno.xls")), row.names=FALSE, sep="\t", quote=FALSE)
        # 写真正的 Excel 文件
        wb <- createWorkbook()
        addWorksheet(wb, "GeneAnno")
        writeData(wb, "GeneAnno", data_tb)
        
        # 自动调整列宽
        setColWidths(wb, "GeneAnno", cols = 1:ncol(data_tb), widths = "auto")
        
        saveWorkbook(wb, 
                    file.path(outdir, paste0(prefix, "_GeneAnno.xlsx")), 
                    overwrite = TRUE)
        ## 如果要在已有xlsx文件中新增sheet
        # wb <- loadWorkbook("已有文件.xlsx")
        # addWorksheet(wb, "新Sheet名")
        return(data_tb)
    }
# =======================================

# ========== GOKEGG富集执行函数 ======================
go_kegg <- function(species, genes_use, tp, outdir1 = "./", spEnrch = FALSE) {
    print(species)
    
    # 物种设置
    if (species == "human") {
        suppressMessages(library(org.Hs.eg.db))
        orgdb <- 'org.Hs.eg.db'
        org <- 'hsa'
        ref <- 'human'
        ENTREZID <- mapIds(org.Hs.eg.db, keys = genes_use, keytype = "SYMBOL", column = "ENTREZID")
    } else if (species == "mouse") {
        suppressMessages(library(org.Mm.eg.db))
        orgdb <- 'org.Mm.eg.db'
        org <- 'mmu'
        ref <- 'mouse'
        ENTREZID <- mapIds(org.Mm.eg.db, keys = genes_use, keytype = "SYMBOL", column = "ENTREZID")
    } else if (species == 'rat') {
        suppressMessages(library(org.Rn.eg.db))
        orgdb <- 'org.Rn.eg.db'
        org <- 'rno'
        ref <- 'rat'
        ENTREZID <- mapIds(org.Rn.eg.db, keys = genes_use, keytype = "SYMBOL", column = "ENTREZID")
    }

    ENTREZID <- ENTREZID[!is.na(ENTREZID)]
    print(orgdb)
    print(org)
    print(length(ENTREZID))

    # ========== GO 富集 ==========
    print("RUN GO!!!")
        outdir_go <- file.path(outdir1, "GO")
        if (!file.exists(outdir_go)) {
            dir.create(outdir_go, recursive = TRUE)
            print(outdir_go)
        }
        
        enrich.go <- enrichGO(gene = genes_use,
                              OrgDb = orgdb,
                              keyType = 'SYMBOL',
                              ont = 'ALL',
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.1,
                              readable = FALSE)

        GO_all <- as.data.frame(enrich.go)
        print(head(GO_all))
        
        if (nrow(GO_all) >= 1) {
            # 计算 RichFactor
            GeneRatio <- as.numeric(str_split(GO_all$GeneRatio, "/", simplify = TRUE)[, 1])
            BgRatio <- as.numeric(str_split(GO_all$BgRatio, "/", simplify = TRUE)[, 1])
            GO_all$term <- paste(GO_all$ID, GO_all$Description, sep = ': ')
            GO_all$RichFactor <- GeneRatio / BgRatio
            
            # 排序并处理 geneID
            GO_all <- GO_all[order(GO_all$ONTOLOGY, GO_all$RichFactor, decreasing = TRUE), ]
            GO_all$geneID <- gsub("/", ",", GO_all$geneID)
            
            # 按 ontology 拆分输出
            for (ont in c("BP", "CC", "MF")) {
                sub_go <- GO_all[GO_all$ONTOLOGY == ont, ]
                if (nrow(sub_go) > 0) {
                    write.table(sub_go, 
                            file.path(outdir_go, paste0(tp, '_go.', ont, '.xls')), 
                            sep = "\t", row.names = FALSE, quote = FALSE)
                    # 拆分输出通路-基因网络图
                    pdf(file.path(,"PathGene_NetPlot.pdf"),width=7,height=7)
                    print(cnetplot(enrich.go, showCategory=10))
                    dev.off()
                    pdf(file.path(,"PathGene_ClusterNet.pdf"),width=7,height=7)
                    print(emapplot(pairwise_termsim(enrich.go), showCategory=30))
                    dev.off()
                }
            }

            # 输出完整 GO 结果
            write.table(GO_all, 
                    file.path(outdir_go, paste0(tp, '_go.ALL.xls')), 
                    sep = "\t", row.names = FALSE, quote = FALSE)

            # GO similarity 分析（可选）
            if (spEnrch) {
                outdir_go_sim <- file.path(outdir_go, "GO_similarity")
                if (!file.exists(outdir_go_sim)) {
                    dir.create(outdir_go_sim, recursive = TRUE)
                }
                print(outdir_go_sim)
                
                for (k in unique(GO_all$ONTOLOGY)) {
                    print(k)
                    GO_sub <- GO_all[which(GO_all$ONTOLOGY == k), ]
                    go_id <- GO_sub$ID
                    print(length(go_id))
                    print(head(go_id))
                    
                    mat <- GO_similarity(go_id, ont = k, db = orgdb, measure = "Rel")
                    write.csv(mat, file.path(outdir_go_sim, paste0(tp, "_similarityGO_", k, "_mat.csv")), row.names = TRUE)
                    print(head(mat))
                    print(paste0("ncol(mat): ", ncol(mat)))
                    
                    if (ncol(mat) > 3) {
                        print(dim(mat))
                        pdf(file.path(outdir_go_sim, paste0(tp, "_simplifyGO_", k, ".pdf")), width = 10, height = 10)
                        print(simplifyGO(mat, plot = TRUE))
                        dev.off()
                        
                        test <- simplifyGO(mat, plot = FALSE)
                        write.csv(test, file.path(outdir_go_sim, paste0(tp, "_simplifyGO_", k, ".csv")))
                        print(head(test))
                    }
                }
            }
        }

    # ========== KEGG 富集（带容错） ==========
        kegg <- NULL
        kegg_df <- data.frame()  # 空数据框占位
        
        if (use_kegg) {
            print("RUN KEGG!!!")
            outdir_kegg <- file.path(outdir1, "KEGG")
            if (!file.exists(outdir_kegg)) {
                dir.create(outdir_kegg, recursive = TRUE)
            }
            
            # 设置超时
            old_timeout <- getOption("timeout")
            options(timeout = kegg_timeout)
            on.exit(options(timeout = old_timeout), add = TRUE)
            
            kegg <- tryCatch({
                enrichKEGG(
                    gene = ENTREZID,
                    organism = org,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.1
                )
            }, error = function(e) {
                message("  KEGG 查询失败: ", e$message)
                NULL
            })
            
            if (!is.null(kegg)) {
                kegg_df <- as.data.frame(kegg)
                print(head(kegg_df))
                
                if (nrow(kegg_df) >= 1) {
                    # 计算 RichFactor
                    kegg_df$RichFactor <- sapply(strsplit(kegg_df$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
                    kegg_df <- kegg_df[order(kegg_df$RichFactor, decreasing = TRUE), ]
                    
                    # 替换 GeneID
                    kegg_df$geneID <- sapply(strsplit(kegg_df$geneID, "/"), function(genes) {
                        symbols <- mapIds(get(orgdb), keys = genes, keytype = "ENTREZID", column = "SYMBOL")
                        paste(na.omit(symbols), collapse = ",")
                    })
                    
                    # 清理 Description
                    kegg_df$Description <- gsub(" - .*? \\(.*\\)", "", kegg_df$Description)
                    
                    # 输出
                    write.table(kegg_df, file.path(outdir_kegg, paste0(tp, '_kegg.xls')), 
                            sep = "\t", quote = FALSE, row.names = FALSE)
                    
                    # 可视化
                    topn <- min(10, nrow(kegg_df))
                    if (topn >= 3) {
                        plot_df <- kegg_df[1:topn, ]
                        plot_df$term <- paste(plot_df$ID, plot_df$Description, sep = ": ")
                        
                        p_kegg <- ggplot(plot_df, aes(x = RichFactor, y = reorder(term, RichFactor))) +
                            geom_point(aes(size = Count, color = p.adjust)) +
                            scale_colour_gradient(low = "red", high = "blue") +
                            theme_bw() +
                            labs(title = paste('Top', topn, 'KEGG Pathways'),
                                x = 'Rich Factor', y = '') +
                            theme(plot.title = element_text(hjust = 0.5))
                        
                        pdf(file.path(outdir_kegg, paste0(tp, '_KEGG_DotPlot.pdf')), width = 8, height = 6)
                        print(p_kegg)
                        dev.off()
                    }
                }
            } else {
                message("  KEGG 富集失败，跳过 KEGG 分析")
            }
        }
        
    # ========== REACTOME 富集 ==========
    print("RUN REACTOME!!!")
    outdir_reactome <- file.path(outdir1, "REACTOME")
    if (!file.exists(outdir_reactome)) {
        dir.create(outdir_reactome, recursive = TRUE)
        print(outdir_reactome)
    }
    
    Reactome <- enrichPathway(
        gene = ENTREZID,
        pvalueCutoff = 1,
        readable = TRUE,
        organism = ref
    ) %>% as.data.frame()
    
    # 输出完整 REACTOME 结果
    write.table(Reactome, 
               file.path(outdir_reactome, paste0(tp, '_Reactome.xls')), 
               sep = "\t", quote = FALSE, row.names = FALSE)

    # REACTOME 可视化
    sigdf <- subset(Reactome, p.adjust < 0.1)
    topn <- min(10, nrow(sigdf))
    
    if (topn >= 3) {
        plotdf <- separate(data = sigdf, col = BgRatio, into = c("BgCount", "BgAll"), sep = "/")
        plotdf$BgCount <- as.integer(plotdf$BgCount)
        plotdf$RichFactor <- plotdf$Count / plotdf$BgCount
        plotdf <- plotdf[order(plotdf$p.adjust), ]
        plotdf <- head(plotdf, topn)
        
        p <- ggplot(data = plotdf, mapping = aes(x = RichFactor, y = reorder(Description, RichFactor))) +
            geom_point(aes(color = p.adjust, size = Count)) +
            scale_colour_gradient(low = "red", high = "blue") +
            theme_bw() +
            labs(title = paste('Top', nrow(plotdf), 'Reactome Pathways'),
                 x = 'Rich factor',
                 y = 'Pathway') +
            theme(plot.title = element_text(hjust = 0.5),
                  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
            scale_y_discrete(labels = function(x) str_wrap(x, width = 50))

        pdf(file.path(outdir_reactome, paste0(tp, '_REACTOME_DotPlot.pdf')), width = 8, height = 6)
        print(p)
        dev.off()
    }

    # ========== GO 可视化（简化） ==========
    if (nrow(GO_all) > 3) {
        print("Plot GO")
        
        # 按 p.adjust 取 top5
        go_padj <- GO_all %>%
            group_by(ONTOLOGY) %>%
            slice_min(order_by = p.adjust, n = 5) %>%
            ungroup() %>%
            arrange(ONTOLOGY, p.adjust)
        go_padj$term <- factor(go_padj$term, levels = go_padj$term)
        
        p_go_padj <- ggplot(go_padj, aes(term, -log10(p.adjust))) +
            geom_col(aes(fill = ONTOLOGY), width = 0.5, show.legend = FALSE) +
            scale_fill_manual(values = c('BP' = '#D06660', 'CC' = '#5AAD36', 'MF' = '#6C85F5')) +
            facet_grid(ONTOLOGY ~ ., scale = 'free_y', space = 'free_y') +
            theme_bw() +
            theme(panel.grid = element_blank(),
                  axis.text = element_text(size = 10),
                  strip.background = element_rect(colour = "white", fill = "white")) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
            scale_x_discrete(labels = function(x) str_wrap(x, width = 70)) +
            coord_flip() +
            labs(x = '', y = '-Log10 P-Value\n')
        
        pdf(file.path(outdir_go, paste0(tp, '_GOplot_Padj.pdf')), width = 10, height = 8)
        print(p_go_padj)
        dev.off()
        
        # 按 RichFactor 取 top5
        go_rich <- GO_all %>%
            group_by(ONTOLOGY) %>%
            slice_max(order_by = RichFactor, n = 5) %>%
            ungroup() %>%
            arrange(ONTOLOGY, desc(RichFactor))
        go_rich$term <- factor(go_rich$term, levels = go_rich$term)
        
        p_go_rich <- ggplot(go_rich, aes(term, RichFactor)) +
            geom_col(aes(fill = ONTOLOGY), width = 0.5, show.legend = FALSE) +
            scale_fill_manual(values = c('BP' = '#D06660', 'CC' = '#5AAD36', 'MF' = '#6C85F5')) +
            facet_grid(ONTOLOGY ~ ., scale = 'free_y', space = 'free_y') +
            theme_bw() +
            theme(panel.grid = element_blank(),
                  axis.text = element_text(size = 10),
                  strip.background = element_rect(colour = "white", fill = "white")) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
            scale_x_discrete(labels = function(x) str_wrap(x, width = 70)) +
            coord_flip() +
            labs(x = '', y = 'Rich Factor\n')
        
        pdf(file.path(outdir_go, paste0(tp, '_GOplot_Genenb.pdf')), width = 10, height = 8)
        print(p_go_rich)
        dev.off()
    }

    return(list(GO = enrich.go, KEGG = kegg, Reactome = Reactome))
}

# ========================================

# ========== 新函数: 实现富集分析 ==========
run_enrich_v2 <- function(genes, bg_matrix = "no", species = "mouse", 
                            label = "", outdir = "./") {
        if (length(genes) < 5) {
            message("    跳过: ", label, " (基因数 ", length(genes), " < 5)")
            return(NULL)
        }
        if(!file.exists(outdir)){
            dir.create(outdir, recursive = TRUE)
        }

        # 物种设置
        if (species == "human") {
            org <- org.Hs.eg.db
            organism <- "hsa"
        } else if (species == "mouse") {
            org <- org.Mm.eg.db
            organism <- "mmu"
        }

        # ID转换
        gene_conv <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org)
        gene_entrez <- na.omit(unique(gene_conv$ENTREZID))
        
        # 背景基因
        if (bg_matrix != "no" && file.exists(bg_matrix)) {
            mtr <- read.csv(bg_matrix, row.names = 1)
            bg_conv <- bitr(rownames(mtr), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org)
            bg_entrez <- na.omit(unique(bg_conv$ENTREZID))
        } else {
            bg_entrez <- NULL
        }

        cat("    ", label, ": ", length(genes), " genes -> ", length(gene_entrez), " ENTREZID\n")

        # GO富集
        go_results <- list()
        for (ont in c("BP", "MF", "CC")) {
            ego <- enrichGO(
                gene = gene_entrez,
                OrgDb = org,
                keyType = "ENTREZID",
                ont = ont,
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2,
                universe = bg_entrez,
                readable = TRUE
            )
            if (!is.null(ego) && nrow(ego@result) > 0) {
                go_results[[ont]] <- ego@result
                write.csv(ego@result, file.path(outdir, paste0("GO_", ont, "_", label, ".csv")), row.names = FALSE)
            }
        }

        # KEGG富集
        ekegg <- enrichKEGG(
            gene = gene_entrez,
            organism = organism,
            pvalueCutoff = 0.05,
            universe = bg_entrez
        )
        ekegg <- setReadable(ekegg, OrgDb = org, keyType = "ENTREZID")
        
        kegg_result <- NULL
        if (!is.null(ekegg) && nrow(ekegg@result) > 0) {
            kegg_result <- ekegg@result
            write.csv(ekegg@result, file.path(outdir, paste0("KEGG_", label, ".csv")), row.names = FALSE)
        }

        # Reactome富集
        ereactome <- enrichPathway(
            gene = gene_entrez,
            organism = ifelse(species == "mouse", "mouse", "human"),
            pvalueCutoff = 0.05,
            universe = bg_entrez
        )
        ereactome <- setReadable(ereactome, OrgDb = org, keyType = "ENTREZID")
        
        reactome_result <- NULL
        if (!is.null(ereactome) && nrow(ereactome@result) > 0) {
            reactome_result <- ereactome@result
            write.csv(ereactome@result, file.path(outdir, paste0("Reactome_", label, ".csv")), row.names = FALSE)
        }

        return(list(
            label = label,
            genes = genes,
            gene_count = length(genes),
            GO = go_results,
            KEGG = kegg_result,
            Reactome = reactome_result
        ))
    }
# ========================================

# ========== 新函数: 构建通路网络 ==========
build_pathway_network <- function(all_enrich_results, outdir) {
        # 提取所有显著通路ID
        pathway_list <- list()
        db_types <- c("GO", "KEGG", "Reactome")
        for (result_name in names(all_enrich_results)) {
            result <- all_enrich_results[[result_name]]
            if (is.null(result)) next
            
            for (db in db_types) {
                if (db == "GO") {
                    for (ont in names(result[[db]])) {
                        df <- result[[db]][[ont]]
                        if (!is.null(df) && nrow(df) > 0) {
                            sig_df <- df[df$p.adjust < 0.05, ]
                            if (nrow(sig_df) > 0) {
                                key <- paste0(result_name, "_", db, "_", ont)
                                pathway_list[[key]] <- data.frame(
                                    result_name = result_name,
                                    db = db,
                                    sub_db = ont,
                                    ID = sig_df$ID,
                                    Description = sig_df$Description,
                                    p.adjust = sig_df$p.adjust,
                                    gene_ratio = sig_df$GeneRatio,
                                    stringsAsFactors = FALSE
                                )
                            }
                        }
                    }
                } else {
                    df <- result[[db]]
                    if (!is.null(df) && nrow(df) > 0) {
                        sig_df <- df[df$p.adjust < 0.05, ]
                        if (nrow(sig_df) > 0) {
                            key <- paste0(result_name, "_", db)
                            pathway_list[[key]] <- data.frame(
                                result_name = result_name,
                                db = db,
                                sub_db = db,
                                ID = sig_df$ID,
                                Description = sig_df$Description,
                                p.adjust = sig_df$p.adjust,
                                gene_ratio = sig_df$GeneRatio,
                                stringsAsFactors = FALSE
                            )
                        }
                    }
                }
            }
        }

        if (length(pathway_list) == 0) {
            message("  无显著通路，跳过网络分析")
            return(NULL)
        }

        # 合并所有通路
        all_pathways <- dplyr::bind_rows(pathway_list)
        
        # 统计高频通路（被多个分析命中的通路）
        pathway_freq <- all_pathways %>%
            group_by(ID, Description, db, sub_db) %>%
            summarise(
                n_results = n_distinct(result_name),
                result_names = paste(unique(result_name), collapse = ";"),
                min_padj = min(p.adjust),
                avg_gene_ratio = mean(as.numeric(sapply(strsplit(gene_ratio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])))),
                .groups = "drop"
            ) %>%
            arrange(desc(n_results), min_padj)

        # 保存高频通路
        write.csv(pathway_freq, file.path(outdir, "Pathway_Frequency.csv"), row.names = FALSE)

        # 构建通路-分析二分网络
        network_edges <- all_pathways %>%
            select(result_name, ID, p.adjust) %>%
            mutate(weight = -log10(p.adjust))

        write.csv(network_edges, file.path(outdir, "Pathway_Network_Edges.csv"), row.names = FALSE)

        # 提取核心高频通路（被≥3个分析命中）
        core_pathways <- pathway_freq %>% filter(n_results >= 3)
        if (nrow(core_pathways) > 0) {
            write.csv(core_pathways, file.path(outdir, "Core_HighFreq_Pathways.csv"), row.names = FALSE)
            message("  高频通路数 (≥3): ", nrow(core_pathways))
        }

        # 提取跨数据库高频通路
        cross_db <- all_pathways %>%
            group_by(ID, Description) %>%
            summarise(
                dbs = paste(unique(db), collapse = ";"),
                n_db = n_distinct(db),
                n_results = n_distinct(result_name),
                .groups = "drop"
            ) %>%
            filter(n_db >= 2) %>%
            arrange(desc(n_results))

        if (nrow(cross_db) > 0) {
            write.csv(cross_db, file.path(outdir, "CrossDB_HighFreq_Pathways.csv"), row.names = FALSE)
            message("  跨数据库高频通路数: ", nrow(cross_db))
        }

        return(list(
            all_pathways = all_pathways,
            pathway_freq = pathway_freq,
            network_edges = network_edges,
            core_pathways = core_pathways,
            cross_db_pathways = cross_db
        ))
    }

# ========================================

# ========== 新函数: 保存富集汇总 ==========
save_enrich_summary <- function(all_enrich_results, pathway_network, outdir) {
        summary_list <- list()
        for (name in names(all_enrich_results)) {
            result <- all_enrich_results[[name]]
            if (is.null(result)) next
            
            n_go <- sum(sapply(result$GO, function(x) ifelse(is.null(x), 0, nrow(x))))
            n_kegg <- ifelse(is.null(result$KEGG), 0, nrow(result$KEGG))
            n_reactome <- ifelse(is.null(result$Reactome), 0, nrow(result$Reactome))
            
            summary_list[[name]] <- data.frame(
                Analysis = name,
                Gene_Count = result$gene_count,
                GO_Pathways = n_go,
                KEGG_Pathways = n_kegg,
                Reactome_Pathways = n_reactome,
                Total_Pathways = n_go + n_kegg + n_reactome,
                stringsAsFactors = FALSE
            )
        }
        
        summary_df <- dplyr::bind_rows(summary_list)
        write.csv(summary_df, file.path(outdir, "Enrichment_Summary.csv"), row.names = FALSE)
        
        message("  富集汇总已保存")
    }
