#!/usr/bin/env Rscript

### 脚本说明
#1) mamba activate r452
#2) 应用: 对接rMATS输出结果, 进行结果过滤, 相关可视化
#3) 更新说明 (v10):
    #           main_func的DSG_DEG部分read_deseq2_results函数:
    #                   调整差异表达基因读入方式, 添加DESeq2输出的qs文件读入方式, 下游接入pairedGSEA
#4) 待更新说明:
    #           富集逻辑调整: 对DSGs上下all差异基因富集，内置在DSGs_DEGs,筛选交集/差集/指定基因通路SBNGview
    #           新增通路网络分析: 基于org.Hs.eg.db注释通路，标注高频富集通路ID
    #           WGCNA

###   模块	                    依赖数据	                     是否可独立运行	    重新运行条件
    # basic_func	            rMATS原始输出	                 ❌ 独立	         更换数据源
    # build_gene_events_cache	summary_table, GR_list_all	    ✅ 独立	        更换阈值或数据
    # event_func (火山图)	    cache	                         ✅ 独立	         更换阈值或标签
    # PCA_PSI_plot	            maser_sig	                    ✅ 独立	        更换数据
    # DEG_DSG	                cache, DESeq2结果(可选)	        ✅ 独立	         更换DEG或阈值
    # batch_gene_specific	    cache, gene_list	            ✅ 独立	        更换基因列表
    # generate_gene_specific	cache, 单基因	                ✅ 独立	         更换单基因

### ========================================================

### 基本分析框架 ============================================
    # main_func
    #     │
    #     基本准备过程 
    #     │── basic_func(infile) 对接rMATS输出结果进行显著性筛选
    #     │       │
    #     │       ├── maser_obj 基于覆盖度筛选后的所有结果, maser格式
    #     │       ├── maser_sig 基于fdr和psi筛选后的结果, maser格式
    #     │       ├── summary_table 汇总成表格格式的所有剪接事件信息
    #     │       └── GR_list_all 显著剪接事件的坐标位点信息, 可对接Gviz画图
    #     │
    #     ├── build_gene_events_cache(basic_result)
    #     │       │
    #     │       └── gene_events_cache 
    #     │                 │
    #     │                 ├── genes 显著差异剪接的基因
    #     │                 ├── event_types 剪接类型
    #     │                 ├── eventLabel 剪接类型:基因_ID
    #     │                 ├── group: 分组设定
    #     │                 ├── by_gene: 依据基因拆分的剪接分析结果
    #     │                 │       ├── gene1: 基因1的事件信息
    #     │                 │       ├── gene2: 基因2的事件信息
    #     │                 │       └── summary_table: 所有事件汇总表 ← NEW
    #     │                 ├── by_event_type: 依据剪接类型拆分信息
    #     │                 ├── fingerprint: 按事件类型的指纹映射 ← UNCHANGED
    #     │                 └── build_time: 构建时间
    #     │
    #     ├── Event_specific.r (cache) ← 参数简化，从cache读取
    #     │       │
    #     │       ├── volcano_plot
    #     │       └── 每个剪接类型的差异上下调基因富集分析 
    #     │
    #     ├── PCA_PSI_plot (maser_sig)
    #     ├── DEG_DSG (cache, deseq2_path[可选]) ← 参数简化
    #     │       │
    #     │       ├── DSG_up/down 分别富集
    #     │       ├── DEG_up/down 分别富集 (如果提供deseq2)
    #     │       ├── DSG+DEG合并富集 (如果提供deseq2)
    #     │       └── 通路网络分析: 高频通路提取 ← NEW
    #     │
    #     └── Gene_specific.r (cache, gene_list, maser_obj[可选])
    #             │
    #             ├── generate_gene_specific
    #             └── generate_sashimi_input
    #
### =======================================================

##### 调用依赖包 ===========================================
    # source("/data/med-wangcq/01CondaEnv/02Git_repo/00MyGit_wchenqi/Multi_omics/SplicingAnalysis/00Packages/01bulk/04_rMATS_filter_pkg/00_Loadpackages.r")
    # source("/data/med-wangcq/01CondaEnv/02Git_repo/00MyGit_wchenqi/Multi_omics/SplicingAnalysis/00Packages/01bulk/03_Enrich_main.r")
    # source("/data/med-wangcq/01CondaEnv/02Git_repo/00MyGit_wchenqi/Multi_omics/SplicingAnalysis/00Packages/01bulk/04_rMATS_filter_pkg/Plot_Script.r")
    # source("/data/med-wangcq/01CondaEnv/02Git_repo/00MyGit_wchenqi/Multi_omics/SplicingAnalysis/00Packages/01bulk/04_rMATS_filter_pkg/Integration.r")
    # source("/data/med-wangcq/01CondaEnv/02Git_repo/00MyGit_wchenqi/Multi_omics/SplicingAnalysis/00Packages/01bulk/04_rMATS_filter_pkg/Event_specific.r")
    # source("/data/med-wangcq/01CondaEnv/02Git_repo/00MyGit_wchenqi/Multi_omics/RNA-seq/02bulkRNA/03Usable_Script/02PlotScript/08Count_Percent_Barplot.r")
    # source("/data/med-wangcq/01CondaEnv/02Git_repo/00MyGit_wchenqi/Multi_omics/SplicingAnalysis/00Packages/01bulk/脚本调整/Enrich_compare_v2.r")
    
    # event_colors <- c(
    #         "SE" = list(c(Inclusion = "#E69F00", Skipping = "#CC79A7")),
    #         "RI" = list(c(Retention = "#E69F00", Non_Retention = "#CC79A7")),
    #         "A3SS" = list(c(A3SS_Short = "#E69F00", A3SS_Long = "#CC79A7")),
    #         "A5SS" = list(c(A5SS_Short = "#E69F00", A5SS_Long = "#CC79A7")),
    #         "MXE" = list(c(MXE_Exon1 = "#E69F00", MXE_Exon2 = "#CC79A7"))
    #     )
##### =====================================================

##### main_func: 分析主框架 [包含所有参数]
    main_func <- function(infile, outdir, gtf_path, conditions,
                          deseq2_path = NULL,
                          enrich_DSG = TRUE,
                          enrich_DEG = FALSE,
                          fdr = 0.05, psi = 0.1, logfc = 1, avg_reads = 5,
                          types = c("SE","A3SS","A5SS","RI","MXE"),
                          ftype = "JC",
                          species = "mouse",
                          gene_list = NULL,
                          cache_file = NULL,
                          rebuild_cache = FALSE,
                          skip_event_analysis = FALSE,
                          verbose = TRUE) {

        cat(" ... 建立文件夹 ... ")
        outdir_rds <- file.path(outdir,"00rds")
        outdir_event <- file.path(outdir,"01EventSpecific")
        outdir_gene <- file.path(outdir,"02GeneSpecific")
        outdir_category <- file.path(outdir, "03CategoryAnalysis")
        outdir_enrich <- file.path(outdir, "04Enrich")

        dir.create(outdir_rds,recursive=TRUE)
        dir.create(outdir_event,recursive=TRUE)
        dir.create(outdir_gene,recursive=TRUE)
        dir.create(outdir_category,recursive=TRUE)
        dir.create(outdir_enrich,recursive=TRUE)
        
        cat(" --- 1. basic_func (数据提取) --- \n")
        basic_result <- basic_func(infile = infile,
                                   outdir = outdir,
                                   conditions = conditions,
                                   fdr = fdr,
                                   psi = psi,
                                   avg_reads = avg_reads,
                                   types = types,
                                   ftype = ftype,
                                   species = species)

        qsave(basic_result, file.path(outdir_rds, "maserfilter_basic_result.qs"))
        
        maser_sig <- basic_result[["maser_sig"]]
        summary_table <- basic_result[["summary_table"]]
        GR_list_all <- basic_result[["GR_list_all"]]

        print(dim(summary_table))
        print(table(summary_table$EvenType,summary_table$UpDown))

        cat(" ---- 2. 构建或加载缓存 ---- \n")
        if (is.null(cache_file)) {
            cache_file <- file.path(outdir_rds, "gene_events_cache.qs")
        }
        if (rebuild_cache && file.exists(cache_file)) {
            file.remove(cache_file)
        }

        # 传入basic_result，内部提取summary_table存入cache
        gene_events_cache <- build_gene_events_cache(
                                    basic_result = basic_result,
                                    cache_file = cache_file,
                                    verbose = verbose,
                                    do_anno = TRUE,
                                    species = species)
        
        # 提取上下调基因
        gg_up <- gene_events_cache$genes$Up
        gg_down <- gene_events_cache$genes$Down
        gg_all <- unique(c(gg_up, gg_down))
        if(enrich_DSG){
        # 对差异剪接分子进行富集分析
            ## 单独对所有Up, Down, All基因进行富集分析
            cat("\n ---- 对所有上下调剪接基因进行富集分析 ---- \n")
            analyze_enrichment(species = species, genes_use = gg_up, tp = "Up", outdir1 = file.path(outdir_enrich, "DSG_sig_Up"), 
                                spEnrch = TRUE, pval_cutoff = 0.05, count_cutoff = 3, top_pathways = 20,
                                top_genes = 20, min_genes = 5)
            analyze_enrichment(species = species, genes_use = gg_down, tp = "Down", outdir1 = file.path(outdir_enrich, "DSG_sig_Down"), 
                                spEnrch = TRUE, pval_cutoff = 0.05, count_cutoff = 3, top_pathways = 20,
                                top_genes = 20, min_genes = 5)
            analyze_enrichment(species = species, genes_use = gg_all, tp = "All", outdir1 = file.path(outdir_enrich, "DSG_sig_All"), 
                                spEnrch = TRUE, pval_cutoff = 0.05, count_cutoff = 3, top_pathways = 20,
                                top_genes = 20, min_genes = 5)
        }
        # 对差异表达分子进行富集分析
        if (enrich_DEG){
            ## 单独对所有Up, Down, All基因进行富集分析
            cat("\n ---- 对所有上下调表达基因进行富集分析 ---- \n")
            deseq2_result <- read_deseq2_results(deseq2_path, fdr = fdr, logfc = logfc)
            deg_up <- deseq2_result[["up_genes"]]
            deg_down <- deseq2_result[["down_genes"]]
            deg_all <- deseq2_result[["all_genes"]]
            cat("  DEG总基因数: ", length(deg_all), "\n")
            cat("  上调DEG基因数: ", length(deg_up), "\n")
            cat("  下调DEG基因数: ", length(deg_down), "\n")
            analyze_enrichment(species = species, genes_use = deg_up, tp = "Up", outdir1 = file.path(outdir_enrich, "DEG_sig_Up"), 
                                spEnrch = TRUE, pval_cutoff = 0.05, count_cutoff = 3, top_pathways = 20,
                                top_genes = 20, min_genes = 5)
            analyze_enrichment(species = species, genes_use = deg_down, tp = "Down", outdir1 = file.path(outdir_enrich, "DEG_sig_Down"), 
                                spEnrch = TRUE, pval_cutoff = 0.05, count_cutoff = 3, top_pathways = 20,
                                top_genes = 20, min_genes = 5)
            analyze_enrichment(species = species, genes_use = deg_all, tp = "All", outdir1 = file.path(outdir_enrich, "DEG_sig_All"), 
                                spEnrch = TRUE, pval_cutoff = 0.05, count_cutoff = 3, top_pathways = 20,
                                top_genes = 20, min_genes = 5)
            ## 这里添加对差异表达和差异剪接分子合并富集
            All_gg <- c(gg_all,deg_all) %>% unique()
            analyze_enrichment(species = species, genes_use = All_gg, tp = "All", outdir1 = file.path(outdir_enrich, "DSG_DEG_sig_All"), 
                               spEnrch = TRUE, pval_cutoff = 0.05, count_cutoff = 3, top_pathways = 20,
                               top_genes = 20, min_genes = 5)
        }
        # 如果deseq2_path是qs文件就启用pairedGSEA分析
        if(tolower(tools::file_ext(deseq2_path)) == "qs"){
            paired_GSEA(deseq2_path = deseq2_path,
                        conditions = conditions,
                        outdir = outdir,
                        species = species)
        }

        cat(" ---- 3. event_func (可跳过) ---- \n")
        if (!skip_event_analysis) {
            # 参数简化：只传cache，函数内部读取summary_table
            for (evt in types) {
                if (verbose) cat("Event specific Analysis: ",evt,"\n")
                event_func(cache = gene_events_cache, 
                           evt = evt, 
                           label = "top5", 
                           outdir = file.path(outdir_event, evt), 
                           fdr = fdr, 
                           psi = psi, 
                           do_enrich = FALSE,
                           species = species)
                PCA_PSI_plot(maser_sig, evt, outdir=file.path(outdir_event,evt))
            }
        } else {
            if (verbose) message("  跳过事件级分析 (火山图/PCA/PSI)")
        }
        if (!is.null(deseq2_path)){
            cat(" ---- DEG_DSG ---- \n")
            # 参数简化：只传cache和deseq2_path
            DEG_DSG(cache = gene_events_cache, 
                    summary_table = summary_table,
                    deseq2_path = deseq2_path, 
                    fdr = fdr, 
                    logfc = logfc, 
                    outdir = outdir_category, 
                    species = species,
                    bg_matrix = "no")
        }
        
        if (verbose) message("\n===== 所有分析完成 =====")
    }

##### =====================================================

##### basic_func: 创建maser对象,基本过滤 [保持不变]
    basic_func <- function(infile, outdir, conditions, fdr=0.05, psi=0.1, 
                           avg_reads=5, types=c("SE","A3SS","A5SS","RI","MXE"), 
                           ftype="JC", species="mouse"){
            condi <- unlist(strsplit(conditions,"vs"))
            print(condi)
            
            if(!file.exists(outdir)){
                dir.create(outdir,recursive=TRUE)
            }
            
            maser_obj <- maser(path=infile, cond_labels=condi, ftype = ftype)
            maser_filter <- filterByCoverage(maser_obj, avg_reads = avg_reads)
            maser_sig <- topEvents(maser_filter,fdr = fdr, deltaPSI = psi)
            
            p <- splicingDistribution(maser_sig)
            pdf(file.path(outdir, "AS_Event_PropSig.pdf"), width = 7, height = 7)
            print(p)
            dev.off()
            
            # Plot_Script.r 脚本调用
            EventType2Group(maser_obj, fdr = 0.05, deltaPSI = 0.1, outdir = outdir, wt = 7,ht = 7)
            EventType2Group_test(maser_obj, fdr = 0.05, deltaPSI = 0.1, outdir = outdir, wt = 7,ht = 7)     

            cat("\n[1/x] 剪接事件输出分析结果合并\n")
            evt_sum <- list()
            GR_list <- list()
            for (evt in types){
                cat("\n 处理剪接类型: ",evt,"\n")
                data <- summary(maser_filter, type = evt)
                if(nrow(data) == 0){
                    cat("\n 没有显著差异事件，跳过\n")
                    next
                }
                data <- data %>% mutate(EvenType = evt,
                                        eventLabel = paste0(EvenType,":",geneSymbol,"_",ID),
                                        negLogFDR = -log10(FDR),
                                        UpDown = "Stable")
                evt_sum[[evt]] <- data

                GR_list[[evt]] <- granges(maser_sig, type = evt) %>% unlist()
                GR_list[[evt]]$eventType <- evt
                GR_list[[evt]]$AS_anno <- names(GR_list[[evt]])
                GR_list[[evt]]$eventLabel <- paste0(GR_list[[evt]]$eventType,":",GR_list[[evt]]$geneSymbol,"_",GR_list[[evt]]$ID)
                print(head(GR_list[[evt]]))
            }
            if(length(evt_sum) == 0) {
                stop("No events found for any type. Check input data.")
            }
            summary_integ <- dplyr::bind_rows(evt_sum)
            GR_list_all <- do.call(c,unname(GR_list))
            sig_labels <- unique(GR_list_all$eventLabel)

            # 标注上下调基因
            if(length(sig_labels) > 0) {
                print("label UpDown!")
                summary_integ$UpDown <- ifelse(summary_integ$eventLabel %in% sig_labels,
                                               ifelse(summary_integ$IncLevelDifference > 0, "Up", "Down"),
                                               "Stable")
            }
            print(head(summary_integ))
            print(table(summary_integ$EvenType, summary_integ$UpDown))
            write.csv(table(summary_integ$EvenType, summary_integ$UpDown),
                      file.path(file.path(outdir,"Event_UpDown_count.csv")),
                      row.names=TRUE,quote=FALSE)

            return(list(maser_obj = maser_filter,
                        maser_sig = maser_sig,
                        summary_table = summary_integ,
                        GR_list_all = GR_list_all
                        ))
    }

##### =============================================================

##### build_gene_events_cache: 构建缓存 [核心改造: summary_table存入by_gene]
    make_fingerprint <- function(evt_type, df) {
            # df 是summary table格式
            # 返回结果为eventLabel为name的向量
            if(evt_type == "SE") {
                fp <- paste0(df$Chr,":", 
                             df$Strand,":",
                             df$exon_target,":",
                             df$exon_upstream,":", 
                             df$exon_downstream)
            } else if(evt_type %in% c("A3SS", "A5SS")) {
                fp <- paste0(df$Chr,":",
                             df$Strand,":",
                             df$exon_long,":",
                             df$exon_short,":",
                             df$exon_flanking)
            } else if(evt_type == "RI") {
                fp <- paste0(df$Chr,":",
                             df$Strand,":",
                             df$exon_ir,":",
                             df$exon_upstream,":",
                             df$exon_downstream)
            } else if(evt_type == "MXE") {
                    fp <- paste0(df$Chr,":",
                             df$Strand,":",
                             df$exon_1,":",
                             df$exon_2,":",
                             df$exon_upstream,":",
                             df$exon_downstream)
            }
            names(fp) <- df$eventLabel
            return(fp)
    }

    build_gene_events_cache <- function(basic_result,
                                        cache_file = NULL,
                                        verbose = TRUE,
                                        do_anno = FALSE,
                                        species = species) {
            # 更新: 在summary table中基于 org.Hs.eg.db/org.Mm.eg.db 添加基因相关通路注释
            if (verbose) message("\n===== 构建基因事件缓存 (v10) =====")
            summary_table <- basic_result[["summary_table"]]
            GR_list_all <- basic_result[["GR_list_all"]]
            maser_sig <- basic_result[["maser_sig"]]
            condition <- factor(
                    c(rep(maser_sig@conditions[1], maser_sig@n_cond1),
                    rep(maser_sig@conditions[2], maser_sig@n_cond2)), 
                    levels = maser_sig@conditions)

            summary_sig <- summary_table[summary_table$UpDown %in% c("Up", "Down"), ]
            sig_ids <- unique(summary_sig$eventLabel)
            sig_ids_up <- unique(summary_sig$eventLabel[which(summary_sig$UpDown == "Up")])
            sig_ids_down <- unique(summary_sig$eventLabel[which(summary_sig$UpDown == "Down")])
            GR_list_all <- GR_list_all[GR_list_all$eventLabel %in% sig_ids, ]
            gr_ids <- unique(GR_list_all$eventLabel)

            if(length(sig_ids) != length(gr_ids)) {
                stop("GR_list_all 和 summary_sig 显著剪接事件数量不一致")
            }

            all_genes <- list()
            all_genes[["Up"]] <- unique(summary_sig$geneSymbol[which(summary_sig$UpDown == "Up")])
            all_genes[["Down"]] <- unique(summary_sig$geneSymbol[which(summary_sig$UpDown == "Down")])
            sig_genes <- unique(summary_sig$geneSymbol)
            event_types <- unique(summary_sig$EvenType)
            
            if (verbose) message("  总基因数: ", length(all_genes))
            cat("Sig Event types: ", event_types, "\n")

            # ---- 初始化缓存结构 ----
            cache <- list(
                genes = all_genes,
                event_types = event_types,
                eventLabel = list(Up = sig_ids_up, Down = sig_ids_down),
                group = condition,
                fingerprint = list(),
                by_gene = list(),
                by_event_type = list(),
                build_time = Sys.time()
            )

            # ---- 1. 构建by_gene ----
            if (verbose) message("  按基因索引事件...")
            cache$by_gene[["gr"]] <- split(GR_list_all, GR_list_all$geneSymbol)
            cache$by_gene[["summary"]] <- split(summary_sig, summary_sig$geneSymbol)         
            ## 这里按照基因拆分,得到每个基因对应的eventlabel,同时添加简介类型标注
            cache$by_gene[["eventLabel"]] <- split(summary_sig$eventLabel, summary_sig$geneSymbol) 
            cache$by_gene[["eventLabel"]] <- by(summary_sig[, c("eventLabel", "EvenType")],
                                                summary_sig$geneSymbol,
                                                FUN = function(x) setNames(x$eventLabel, x$EvenType))
            # 逐个基因添加基因注释信息,仅在有显著差异剪接基因的时候构建 cache$by_gene[["fun_anno"]]
            ## 原来的存储结构是单独构建list存储分Up, Down, All的注释结果
            # genelist <- list(Up = unique(all_genes$Up),
            #                  Down = unique(all_genes$Down),
            #                  All = c(all_genes$Up,all_genes$Down) %>% unique())
            # cache$fun_anno <- setNames(lapply(genelist, function(gene) {
            #                             if (length(gene) == 0) next
            #                             if (verbose) message("  逐个基因获取功能注释 (", length(gene), " 个基因)...")
            #                             gene_Anno(gene, species = species, anno_col = c("ENSEMBL", "ENTREZID", "GO", "ONTOLOGY", "PATH"))
            #                         }), names(genelist))
            
            ## 调整后结构: 存到cache$by_gene[["fun_anno"]]里面
            gene <- c(all_genes$Up,all_genes$Down) %>% unique()
            cache$by_gene[["fun_anno"]] <- gene_Anno(gene, species = species, anno_col = c("ENSEMBL", "ENTREZID", "GO", "ONTOLOGY", "PATH"))

            # ---- 2. 构建 by_event_type ----
            if (verbose) message("  构建事件类型索引...")
            # 后面event specific分析只需要调用summary table
            cache$by_event_type <- split(summary_table, summary_table$EvenType)     # 注意这里为了对接火山图, 需要保留所有结果
            for (evt in event_types) {
                cat("Event type: ",evt,"\n")
                summary_evt <- cache$by_event_type[[evt]]
                if (is.null(summary_evt) || nrow(summary_evt) == 0) next
                # ---- 3. 构建fingerprint ----
                cache$fingerprint[[evt]] <- make_fingerprint(evt, summary_evt)
                # 如果要提取特定基因的fp,可以从cache$by_gene[["eventLabel"]][[gene]]中提取 eventLabel
                # 然后和cache$fingerprint[[evt]]的names取交集, 得到坐标
            }

        if (verbose) {
            message("  上下调差异剪接基因数:")
            print(sapply(cache$genes,length))
            message("  所有剪接事件类型分布:")
            print(sapply(cache$fingerprint,length))
        }

        # ---- 5. 保存到磁盘 ----
        if (!is.null(cache_file)) {
            dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
            qsave(cache, cache_file)
            if (verbose) message("  缓存已保存至: ", cache_file)
        }
        
        return(cache)
    }

##### =====================================================

##### DEG/DSG分类分析
    # 输出结果: 
    #         1) 交集并集韦恩图
    #         2) 交集并集通路整合
    #         3) WGCNA分析结果
    read_deseq2_results <- function(deseq2_path,
                                    fdr = 0.05,
                                    logfc = 1,
                                    species = species) {
        if (is.null(deseq2_path) || deseq2_path == "" || !file.exists(deseq2_path)) {
            message("  DESeq2结果文件未提供或不存在，跳过DEG分析")
            return(NULL)
        }

        message("  读取DESeq2结果: ", deseq2_path)

        ext <- tolower(tools::file_ext(deseq2_path))
        if (ext %in% c("tsv", "txt")) {
            deseq_df <- read.delim(deseq2_path, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
            if (!"geneSymbol" %in% colnames(deseq_df)) {
                deseq_df$geneSymbol <- rownames(deseq_df)
            }
        } else if(ext == "qs"){
            diff <- qread(deseq2_path)
            deseq_df <- results(diff, alpha=0.1)
            gg <- str_split(rownames(deseq_df),":",simplify=TRUE)[,1]
            # 暂时现在只支持mouse和human
            if (species == "mouse"){
                deseq_df$geneSymbol <- convertid::todisp2(gg, biom.data.set = "mmusculus_gene_ensembl")
            }else if(species == "human"){
                deseq_df$geneSymbol <- convertid::todisp2(gg)
            }
        } else{
            deseq_df <- read.csv(deseq2_path, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
            if (!"geneSymbol" %in% colnames(deseq_df)) {
                deseq_df$geneSymbol <- rownames(deseq_df)
            }
        }

        message("  DESeq2总行数: ", nrow(deseq_df))

        col_map <- list(
            padj = c("padj","FDR", "fdr", "adj_pvalue","p_val_adj"),
            log2FoldChange = c("log2FoldChange", "log2fc", "log2FC", "log2_fold_change", "avg_log2FC")
        )

        for (std_name in names(col_map)) {
            found <- base::intersect(colnames(deseq_df), col_map[[std_name]])
            if (length(found) > 0) {
                colnames(deseq_df)[colnames(deseq_df) == found[1]] <- std_name
            }
        }

        deseq_df <- deseq_df[!is.na(deseq_df$padj) & !is.na(deseq_df$log2FoldChange), ]

        deg_df <- deseq_df[deseq_df$padj <= fdr & abs(deseq_df$log2FoldChange) >= logfc, ]

        message("  DEG筛选条件: padj < ", fdr, ", |log2FC| >= ", logfc)
        message("  差异表达基因数: ", nrow(deg_df))

        up_genes <- unique(deg_df$geneSymbol[deg_df$log2FoldChange > 0])
        down_genes <- unique(deg_df$geneSymbol[deg_df$log2FoldChange < 0])
        all_genes <- unique(c(up_genes, down_genes))

        result <- list(
            raw_df = deseq_df,
            deg_df = deg_df,
            up_genes = up_genes,
            down_genes = down_genes,
            all_genes = all_genes,
            fdr = fdr,
            logfc = logfc
        )

        return(result)
    }

    categorize_deg_dsg <- function(dsg_genes, deg_genes, outdir, prefix) {
        # 建立输出文件夹
        outdir1 <- file.path(outdir, prefix)
        if(!file.exists(outdir1)){
            dir.create(outdir1, recursive=TRUE)
        }

        dsg_genes <- unique(dsg_genes)
        deg_genes <- unique(deg_genes)

        both_genes <- base::intersect(dsg_genes, deg_genes)
        only_deg <- base::setdiff(deg_genes, dsg_genes)
        only_dsg <- base::setdiff(dsg_genes, deg_genes)

        # 创建list用于后续绘图
        result <- list(
            Both = both_genes,
            Only_DEG = only_deg,
            Only_DSG = only_dsg
            )
        
        # 整理结果输出文本
        summary_df <- data.frame(
            Category = c("Both", "Only_DEG", "Only_DSG"),
            Genes = c(
                paste0(result$Both, collapse = ","),
                paste0(result$Only_DEG, collapse = ","),
                paste0(result$Only_DSG, collapse = ",")
            ),
            Gene_Count = c(
                length(both_genes),
                length(only_deg),
                length(only_dsg)
            ),
            Description = c(
                "有差异表达和差异剪接",
                "有差异表达无差异剪接",
                "无差异表达有差异剪接"
            )
        )
        write.csv(summary_df, file.path(outdir1, paste0(prefix, "_Category_Summary.csv")), row.names = FALSE)

        return(result)
    }

    # 主脚本:     
    DEG_DSG <- function(cache, summary_table, deseq2_path = NULL,
                        fdr, logfc, outdir, species, bg_matrix = "no") {
        
        cat("\n ---- DEG/DSG分类与富集分析 ---- \n")

        cat("\n ---- 针对基因的处理 ---- \n")
        ## 提取DSG基因（按上下调拆分）
        dsg_up <- summary_table$geneSymbol[summary_table$UpDown == "Up"] %>% unique()
        dsg_down <- summary_table$geneSymbol[summary_table$UpDown == "Down"] %>% unique()
        dsg_all <- c(dsg_up,dsg_down) %>% unique()
        cat("  DSG总基因数: ", length(dsg_all), "\n")
        cat("  上调DSG基因数: ", length(dsg_up), "\n")
        cat("  下调DSG基因数: ", length(dsg_down), "\n")

        # 写出显著差异剪接矩阵表格
        outdir_all <- file.path(outdir, "DSG_DEG_All")
        if(!file.exists(outdir_all)){
            dir.create(outdir_all, recursive = TRUE)
        }
        summary_sig <- summary_table[summary_table$UpDown %in% c("Up","Down"),]
        if (!is.null(summary_sig) && nrow(summary_sig) != 0) {
            write.table(summary_sig, file.path(outdir_all, "DSG_sig_result.xls"), sep = "\t", quote = FALSE, row.names = FALSE)
        }

        ## 读取差异基因文件, 提取上下调差异基因
        deseq2_result <- read_deseq2_results(deseq2_path, fdr = fdr, logfc = logfc)
        deg_up <- deseq2_result[["up_genes"]]
        deg_down <- deseq2_result[["down_genes"]]
        deg_all <- deseq2_result[["all_genes"]]
        cat("  DEG总基因数: ", length(deg_all), "\n")
        cat("  上调DEG基因数: ", length(deg_up), "\n")
        cat("  下调DEG基因数: ", length(deg_down), "\n")
        
        # 写出显著差异表达矩阵表格
        if (!is.null(deseq2_result$deg_df)) {
            write.csv(deseq2_result$deg_df, file.path(outdir_all, "DEG_sig_result.csv"), quote=FALSE, row.names = TRUE)
        }
        ## 对差异剪接和差异表达基因取交集画Venn图
        DSG_DEG_Up <- categorize_deg_dsg(dsg_genes = dsg_up, deg_genes = deg_up, outdir, prefix = "DSG_DEG_Up")
        DSG_DEG_Down <- categorize_deg_dsg(dsg_genes = dsg_down, deg_genes = deg_down, outdir, prefix = "DSG_DEG_Down")
        DSG_DEG_All <- categorize_deg_dsg(dsg_genes = dsg_all, deg_genes = deg_all, outdir, prefix = "DSG_DEG_All")
        # 韦恩图
        Venn_plot(sets = list(DSG_up = dsg_up, DEG_up = deg_up), outdir = file.path(outdir,"DSG_DEG_Up"),prefix = "DSG_DEG_Up",width = 8,height = 8)
        Venn_plot(sets = list(DSG_down = dsg_down, DEG_down = deg_down), outdir = file.path(outdir,"DSG_DEG_Down"),prefix = "DSG_DEG_Down",width = 8,height = 8)
        Venn_plot(sets = list(DSG_all = dsg_all, DEG_all = deg_all), outdir = file.path(outdir,"DSG_DEG_All"),prefix = "DSG_DEG_All",width = 8,height = 8)

        cat("\n ---- 针对功能注释结果的处理 ---- \n")
        cat("\n ---- 注释 DEG_sig 的功能, 写出文件 ---- \n")
        deg_FuncAnno <- gene_Anno(genes = deg_all, species = "mouse", anno_col = c("ENSEMBL", "ENTREZID", "GO", "ONTOLOGY", "PATH"))
        deg_tb <- gene_Anno_to_df(fun_anno = deg_FuncAnno , outdir = file.path(outdir,"DSG_DEG_All"), prefix = "DEG_all")
        
        cat("\n ---- 提取 DSG_sig 注释功能, 写出文件 ---- \n")
        dsg_FuncAnno <- cache$by_gene$fun_anno       # dsg的注释结果可以直接提取
        dsg_tb <- gene_Anno_to_df(fun_anno = dsg_FuncAnno , outdir = file.path(outdir,"DSG_DEG_All"), prefix = "DSG_all")

        ## 对差异剪接和差异表达基因功能注释结果取交集画Venn图
        # ---- 定义分析组合 ----
        gene_sets <- list(Up = list(dsg = dsg_up, deg = deg_up),
                          Down = list(dsg = dsg_down, deg = deg_down),
                          All = list(dsg = dsg_all, deg = deg_all)
                        )
        anno_types <- list(GO = "GO", KEGG = "PATH")
        # ---- 嵌套循环执行 ----
        for (anno_name in names(anno_types)) {
            type_enrich <- anno_types[[anno_name]]
            cat("\n---- Venn of ",anno_name," - ",type_enrich," ----\n")
            for (set_name in names(gene_sets)) {
                print(set_name)
                dsg_genes <- gene_sets[[set_name]]$dsg
                deg_genes <- gene_sets[[set_name]]$deg
                
                # 提取注释
                dsg_anno <- dsg_tb[which(dsg_tb$gene %in% dsg_genes),type_enrich] %>% 
                                na.omit() %>% unique()
                dsg_anno <- unlist(strsplit(dsg_anno,","))
                deg_anno <- deg_tb[which(deg_tb$gene %in% deg_genes),type_enrich] %>% 
                                na.omit() %>% unique()
                deg_anno <- unlist(strsplit(deg_anno,","))
                # 建立list
                plot_list <- setNames(list(dsg_anno, deg_anno),c(paste0("DSG_",set_name,"_",anno_name),paste0("DEG_",set_name,"_",anno_name)))
                # 绘图（直接用 paste0 生成文件名）
                Venn_plot(sets = plot_list, 
                          outdir = file.path(outdir, paste0("DSG_DEG_", set_name)), 
                          prefix = paste0("DSG_DEG_", set_name, "_", anno_name), 
                          width = 8, 
                          height = 8)
            }
        }
        if (verbose) message("\n===== DEG/DSG分析完成 =====")
    }

##### =======================================================

##### 推荐分析流 =========================================================
    # # 1. 首次完整分析
    
    # 
    # # 2. 仅重新运行DEG_DSG分析（使用已有cache）
    
    # 
    # # 3. 仅重新运行event_func（使用已有cache）
    
##### ====================================================================
