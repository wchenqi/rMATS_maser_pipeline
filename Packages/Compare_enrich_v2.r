#!/usr/bin/env Rscript

#### 脚本说明：
#1) 运行环境: conda activate /data/med-hancs/apps/anaconda3/2022.10/envs/BasicR
#2) 应用: 基于cache的独立富集分析接口
#3) 更新: v2版本，统一调用run_enrich_v2，支持从cache直接读取

### 加载包
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(enrichplot)

### 从cache直接进行富集分析（最小参数）
run_enrich_from_cache <- function(cache,
                                  deseq2_path = NULL,
                                  bg_matrix = "no",
                                  species = "mouse",
                                  outdir = "./") {
    
    # 从cache读取DSG基因
    if (!"summary_table" %in% names(cache$by_gene)) {
        stop("cache中未找到summary_table")
    }
    summary_table <- cache$by_gene$summary_table
    
    dsg_up <- summary_table$geneSymbol[summary_table$UpDown == "Up"] %>% unique()
    dsg_down <- summary_table$geneSymbol[summary_table$UpDown == "Down"] %>% unique()
    dsg_all <- unique(c(dsg_up, dsg_down))
    
    # 调用DEG_DSG（使用新逻辑）
    DEG_DSG(cache = cache,
            deseq2_path = deseq2_path,
            fdr = 0.05,
            logfc = 1,
            outdir = outdir,
            species = species,
            bg_matrix = bg_matrix)
}

### 旧接口保留（兼容性）
run_enrich <- function(genes, bg_matrix = "no", species = "mouse") {
    return(run_enrich_v2(genes = genes, bg_matrix = bg_matrix, species = species))
}