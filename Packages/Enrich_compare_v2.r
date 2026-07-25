#!/usr/bin/env Rscript

#### 脚本说明：
#1) 运行环境: mamba activate r452
#2) 应用: 包含两个部分
#         基于两个enrichGO和两个enrichKEGG结果对富集到通路进行筛选
#         提供基因列表,进行单基因的功能注释，输出注释
#3) 更新记录:
#         添加pairedGSEA对差异表达转录本和差异剪接基因进行功能差异分析, 在DEG_DSG中调用

### 加载包
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(org.Mm.eg.db)
# library(enrichplot)
library(qs)
library(writexl,lib.loc = "/data/med-wangcq/01CondaEnv/00DataBase/00Tools/seeksoultools.1.2.0/lib/R/library")
library(pairedGSEA)
library(dplyr)

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

### 使用pairedGSEA进行比较富集分析
# 输入要求:
# DESeq2输出DESeqDataSet对象必须是转录本水平结果
# 输出结果在04Enrich文件夹中建立pairedGSEA
# 测试参数
# indir <- "/scratch/2026-07-13/med-wangcq/Others/AnJQ_source/AnalysisData/Sam68KOWT_TAC_Day3_4W_RNAseq/08DESeq2_diffGene/pairedGSEA/"
# outdir <- "/scratch/2026-07-13/med-wangcq/Others/AnJQ_source/AnalysisData/Sam68KOWT_TAC_Day3_4W_RNAseq/06rMATS_filter/02_New/"
# groups <- "Sham_4W_KOvsSham_4W_CTR,TAC_4W_CTRvsTAC_D3_CTR,TAC_4W_KOvsTAC_4W_CTR,TAC_D3_CTRvsSham_4W_CTR,TAC_D3_KOvsTAC_D3_CTR"
# species <- "Mus musculus"             # "Homo sapiens"
# collection <- "M5"
# metadata <- NULL
# covariates <- NULL

# # 基本处理
#1) 适用转录本层面的差异剪接和表达富集结果对比
paired_GSEA <- function(deseq2_path, conditions, outdir, species){
        # 这里的outdir对应main_func里面定义的outdir_enrich
        cat("Compare Groups: ",conditions,"\n")
        object <- qread(deseq2_path)
        gs <- unlist(strsplit(conditions,"vs"))
        # 建立文件夹
        outdir_enrich <- file.path(outdir, "04Enrich", "pairedGSEA")    # 存储输出富集结果文件
        outdir_rds <- file.path(outdir, "00rds")                        # 用来存储输出rds文件
        if(!file.exists(outdir_enrich)){
            dir.create(outdir_enrich, recursive = TRUE)
        }
        cat("Out Path: ",outdir_enrich, "\n")
        diff_results <- paired_diff(object,
                                    group_col = "group",
                                    sample_col = "sample",
                                    baseline = gs[1],
                                    case = gs[2],
                                    metadata = NULL,
                                    covariates = NULL,
                                    experiment_title = conditions,
                                    store_results = FALSE,
                                    use_limma = FALSE,
                                    prefilter = 10,
                                    quiet = FALSE,
                                    parallel = FALSE,
                                    BPPARAM = Bioconditionsarallel::bpparam(),
                                    expression_only = FALSE,
                                    custom_design = FALSE)
        ## 输出结果为dataframe
        qsave(diff_results, file.path(outdir_rds, "Paired_diff.qs"))            # 存储路径再定

        ## Define gene sets in your preferred way
        # 建立存储对象
        ora_list <- list()
        if (species == "mouse"){
            db_species <- "MM"
            collection_all <- c("MH","M1","M2","M3","M5","M7","M8")
        }else if(species == "human"){
            db_species <- "HS"
            collection_all <- c("H","C1","C2","C3","C4","C5","C6","C7","C8","C9")
        }
        for (collection in collection_all){
            print(collection)
            ## 获取通路集合
            gene_sets <- pairedGSEA::prepare_msigdb(
                                            species = species,
                                            db_species = db_species,
                                            collection = collection,
                                            gene_id_type = "ensembl_gene")
            ## Running over-representation analyses, Joining result
            ora <- paired_ora(paired_diff_result = diff_results,
                            gene_sets = gene_sets,
                            cutoff = 0.05,
                            min_size = 25,
                            experiment_title = conditions,
                            expression_only = FALSE,
                            quiet = FALSE)
            ## 存到ora_list里面
            ora_list[[collection]] <- ora

            ## 可视化结果
            # 首先需要判断一下是否可绘图
            ora_sig <- ora[which(ora$padj_expression <= 0.05 | ora$padj_splicing <= 0.05),] %>% as.data.frame()
            if(nrow(ora_sig) > 0){
                writexl::write_xlsx(ora_sig, path = file.path(outdir_enrich, paste0("pairedGSEA_ora_", db_species, "_", collection, "_sig.xlsx")))
                p <- plot_ora(ora,
                            pattern = NULL,
                            paired = FALSE,
                            plotly = FALSE,
                            cutoff = 0.05,
                            lines = TRUE,
                            colors = c("gray40", "#D62728", "#1F77B4", "#FF7F0E")
                            ) + ggplot2::theme_classic()

                pdf(file.path(outdir_enrich,paste0("pairedGSEA_ora_",db_species,"_",collection,".pdf")),width=6,height=5)
                print(p)
                dev.off()
            }
            qsave(ora_list, file.path(outdir_rds, "Paired_oralist.qs"))
        }
}
#2) 可以对接多组富集分析上游输出结果, 实现多组间对比
# EnrichTable为富集输出表格文件
paired_GO <- function(EnrichTable, outdir="./", species="mouse", prefix="no"){
    # 1. 检查并创建输出目录
    if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
    
    # 2. 将你的大表格拆分为按通路 ID 的基因列表
    # 假设你的表格包含列：ID (通路ID), geneID (逗号分隔的基因名)
    # 将原始 dataframe 按 ID 列拆分成一个 list (每行是一个通路)
    # 注意：strsplit 的第二个参数必须是分隔符","
    genelist <- split(EnrichTable$geneID, EnrichTable$ID)
    
    # 3. 清理格式：把逗号分隔的字符串变为纯向量，并去重
    genelist <- lapply(genelist, function(x) {
        # x 是字符串 "GeneA,GeneB,GeneC"
        unique(unlist(strsplit(as.character(x), ",")))
    })
    
    # 4. 根据物种选择数据库
    if(species == "mouse"){
        org <- org.Mm.eg.db
    } else if(species == "human"){
        org <- org.Hs.eg.db
    } else {
        stop("Species must be 'mouse' or 'human'")
    }
    
    # 5. 使用 compareCluster 批量跑 GO 富集 (注意 fun 是单个字符串)
    cat("正在进行多组 GO 富集比较...\n")
    xx <- compareCluster(geneClusters = genelist,
                         fun = "enrichGO",
                         OrgDb = org,
                         keyType = "SYMBOL",
                         ont = "ALL",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2)
    
    # 6. 计算所有富集结果之间的通路相似度（画 emapplot 的前提）
    xx2 <- pairwise_termsim(xx)
    
    # 7. 画图并保存
    # 节点大小按连线数量计算，多组时节点会自动变成对比饼图
    p <- emapplot(xx2, showCategory = 15) 
    pdf(file = file.path(outdir, ifelse(prefix=="no", "PathGene_emapplot.pdf", paste0(prefix,"_PathGene_emapplot.pdf"))), width = 8, height = 8)
    print(p)
    dev.off()

    p <- cnetplot(genelist,
                  layout = "nicely",
                  showCategory = 5,
                  color_category = "#E5C494",
                  size_category = 1,
                  color_item = "#B3B3B3",
                  size_item = 1,
                  color_edge = "grey",
                  size_edge = 0.5,
                  node_label = "all",
                  foldChange = NULL,
                  fc_threshold = NULL,
                  hilight = "none",
                  hilight_alpha = 0.3)
    pdf(file = file.path(outdir, ifelse(prefix=="no", "PathGene_cnetplot.pdf", paste0(prefix,"_PathGene_cnetplot.pdf"))), width = 10, height = 10)
    print(p)
    dev.off()
    cat("绘图完成！图片已保存至:", outdir, "\n")
}