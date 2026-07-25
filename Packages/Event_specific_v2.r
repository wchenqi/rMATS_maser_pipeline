# source("/data/med-wangcq/01CondaEnv/02Git_repo/00MyGit_wchenqi/Multi_omics/SplicingAnalysis/00Packages/01bulk/03_Enrich_main.r")
# source("/data/med-wangcq/01CondaEnv/02Git_repo/00MyGit_wchenqi/Multi_omics/SplicingAnalysis/00Packages/01bulk/04_rMATS_filter_pkg/Plot_Script.r")

##### event_func: 实现单个剪接类型的分析 [参数简化: 从cache读取summary_table]
    event_func <- function(cache, evt, label = "top5", outdir, do_enrich= TRUE,
                           fdr = 0.05, psi = 0.1, species = "mouse") {
        if (!file.exists(outdir)) {
            dir.create(outdir, recursive = TRUE)
        }
        cat("\nNew Dir: ", outdir, "\n")
        ## 分事件提取summary_table
        summary_table <- cache$by_event_type[[evt]]
        if (is.null(summary_table) || nrow(summary_table) == 0) {
            cat("No event found/n")
            next
        }
        ## 火山图
        print("Volcano Plot")
        volcano_plot(summary_table, label = label, evt = evt, outdir = outdir, fdr = fdr, psi = psi)
        
        ## 判断是否进行富集分析
        if(do_enrich){
            # 提取上下调基因
            evt_data <- summary_table[summary_table$FDR <= fdr & abs(summary_table$IncLevelDifference) >= psi, ]
            gg_up <- evt_data$geneSymbol[evt_data$IncLevelDifference > 0]
            gg_down <- evt_data$geneSymbol[evt_data$IncLevelDifference < 0]
            gg_all <- unique(c(gg_up, gg_down))
            cat("    上调基因数: ", length(gg_up), "\n")
            cat("    下调基因数: ", length(gg_down), "\n")
            # 建立文件夹
            outdir_enrich <- file.path(outdir, "Enrichment")
            print(outdir_enrich)
            if (!file.exists(outdir_enrich)) {
                dir.create(outdir_enrich, recursive = TRUE)
            }
            # 使用新版统一富集接口
            analyze_enrichment(species = species, genes_use = gg_up, tp = "Up", outdir1 = file.path(outdir_enrich, "Up"), 
                            spEnrch = TRUE, pval_cutoff = 0.05, count_cutoff = 3, top_pathways = 20,
                            top_genes = 20, min_genes = 5)
            analyze_enrichment(species = species, genes_use = gg_down, tp = "Down", outdir1 = file.path(outdir_enrich, "Down"), 
                            spEnrch = TRUE, pval_cutoff = 0.05, count_cutoff = 3, top_pathways = 20,
                            top_genes = 20, min_genes = 5)
            analyze_enrichment(species = species, genes_use = gg_all, tp = "All", outdir1 = file.path(outdir_enrich, "All"), 
                            spEnrch = TRUE, pval_cutoff = 0.05, count_cutoff = 3, top_pathways = 20,
                            top_genes = 20, min_genes = 5)
        }
    }
##### =======================================================