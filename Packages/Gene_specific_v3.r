##### 脚本实现基因特异分析 ===========================================

##1) 主函数 [兼容改造: 排除summary_table]
    generate_gene_specific <- function(cache = NULL,
                                       maser_obj = NULL,
                                       gene_name,
                                       ens_gtf = NULL,
                                       outdir,
                                       event_types = c("SE","A3SS","A5SS","RI","MXE"),
                                       verbose = TRUE) {
                if (verbose) message("\n===== 处理基因: ", gene_name, " =====")
                valid_genes <- c(cache$genes$Up,cache$genes$Down)
                
                if (!gene_name %in% valid_genes) {
                    stop("  无有效基因可分析")
                }
                if (is.null(cache)) {
                    stop("必须提供 cache 对象")
                }
                
                # 提取基因对应剪接事件 eventLabel
                gene_events <- cache$by_gene$eventLabel[[gene_name]]
                if (is.null(gene_events) || length(gene_events) == 0) {
                    stop("  基因 ", gene_name, " 无剪接事件")
                }
                
                if (length(gene_events) > 0) {
                    gene_events <- gene_events[names(gene_events) %in% event_types]
                    if (length(gene_events) == 0) {
                        if (verbose) message("  基因 ", gene_name, " 无指定剪接类型的事件")
                        next   # ✅ 跳过当前循环，继续下一个基因
                    }
                }
                
                Event_count <- table(names(gene_events))
                cat("  发现剪接类型数量分布：\n")
                print(Event_count)
                
                ## 建立输出文件夹
                outdir_gene <- file.path(outdir, gene_name)
                if (!dir.exists(outdir_gene)) {
                    dir.create(outdir_gene, recursive = TRUE)
                }
                # 根据输入文件判断是否基于maser对象绘图
                use_maser <- !is.null(maser_obj) && !is.null(ens_gtf)
                if (use_maser) {
                    # 使用geneEvents(maser内置函数调用指定基因的maser obj)
                    gene_events_maser <- geneEvents(maser_obj, geneS = gene_name)
                    if (verbose) message("  使用 maser_obj 生成转录本结构图")
                    for (evt in names(gene_events)) {
                        # 这里提取基因已有显著剪接类型
                        evt_data <- gene_events[[evt]]
                        # 如果该类型没有剪接事件，跳过
                        if (length(evt_data) == 0) next
                        # 对类型内所有剪接事件循环
                        for (event_id in evt_data) {
                            event_id <- str_split(event_id,"_",simplify=TRUE)[,2]
                            print(event_id)
                            tryCatch({
                                    pdf(file.path(outdir_gene, paste0("Transcripts_", gsub(":","_",event_id), ".pdf")), width = 6, height = 5)
                                    print(plotTranscripts_with_exon(
                                            events = gene_events_maser,
                                            type = evt,
                                            event_id = as.numeric(event_id),
                                            gtf = ens_gtf,
                                            zoom = TRUE,
                                            show_PSI = FALSE,
                                            show_exon_number = TRUE
                                    ))
                                    dev.off()
                                }, error = function(e) {
                                    if (verbose) message("      绘图失败 ", evt, "_", event_id, ": ", e$message)
                                    while (dev.cur() > 1) dev.off()
                            })
                        }
                    }
                } else {
                        if (verbose) message("    ", evt, ": 跳过 track 图绘制")
                }
                
                return(list(
                    gene = gene_name,
                    output_dir = outdir_gene,
                    event_types = names(gene_events)
                ))
            }

##2) 辅助函数 [保持不变]
        # 2-1) 绘制带外显子编号的转录本结构图
    #'
    #' 扩展plotTranscripts函数，增加外显子编号显示功能。
    #'
    #' @param events MASER对象
    #' @param type 字符，事件类型
    #' @param event_id 字符，事件ID
    #' @param gtf GRanges对象，GTF注释
    #' @param zoom 逻辑，是否缩放，默认FALSE
    #' @param show_PSI 逻辑，是否显示PSI，默认FALSE
    #' @param show_exon_number 逻辑，是否显示外显子编号，默认TRUE
    #' @param exon_number_style 字符，编号样式："numeric","Exon_","E"
    #'
    #' @return 无返回值（绘图到设备）或GRanges（如果return_gr=TRUE）
    #' @export
    ## 测试参数:
    # gg <- "Snhg17"
    # event_id <- "3088"
    # type <- "RI"
    # exon_number_style <- "Exon_"
    # gtf_path <- "/data/med-wangcq/01CondaEnv/02Git_repo/01DataBase/Genome_Annotation_Reference/00Download/Ensembl/Mus_musculus.GRCm38.102.gtf"
    # gtf <- import.gff(gtf_path)
    # basic_re <- qread("/scratch/2026-06-15/med-wangcq/SelfUse/Heart/01_GEO/GSE153801/GSE153800_RNAseq/05rMATS_filter/WT_D0vsSam68KO_D0/00rds/maserfilter_basic_result.qs")
    # maser_sig <- basic_re[["maser_sig"]]
    # events <- geneEvents(maser_sig, geneS = gg)
    plotTranscripts_with_exon <- function(events,
                                          type = c("A3SS", "A5SS", "SE", "RI", "MXE"),
                                          event_id,
                                          gtf,
                                          zoom = FALSE,
                                          show_PSI = FALSE,
                                          show_exon_number = TRUE) {
        
        # ---- 1. 参数检查（保持不变） ----
        is_strict <- TRUE
        exon_number_style <- "Ensembl"
        if (!is(events, "Maser")) {
            stop("Parameter events has to be a maser object.")
        }
        if (!class(gtf) == "GRanges") {
            stop(cat("\"gtf\" should be a GRanges class."))
        }
        if (any(!grepl("chr", Seqinfo::seqlevels(gtf)))) {
            Seqinfo::seqlevels(gtf) <- paste0("chr", Seqinfo::seqlevels(gtf))
        }
        std_chr <- c(paste0("chr", seq(1:22)), "chrX", "chrY")      # 标准染色体命名
        if (any(!Seqinfo::seqlevels(gtf) %in% std_chr)) {
            Seqinfo::seqlevels(gtf, pruning.mode = "coarse") <- intersect(seqlevels(gtf), std_chr)      # 粗放修剪, 只保留标准命名染色体
        }
        type <- match.arg(type)
        
        annot <- slot(events, paste0(type, "_events"))
        if (length(unique(annot$geneSymbol)) > 1) {
            stop(cat("Multiple genes found. Use geneEvents() to select \n             gene-specific AS events."))
        }
        
        # ---- 2. 提取基因名 ----
        gene_name <- unique(annot$geneSymbol)[1]
        
        # ---- 3. 提取特定类型和事件ID的GRanges ----
        grl <- slot(events, paste0(type, "_gr"))
        eventGr_list <- lapply(grl, function(exon) {
            return(exon[exon$ID == event_id, ])
        })
        # 构建 GRangesList 并保留名称
        eventGr <- GRangesList(eventGr_list)
        names(eventGr) <- names(eventGr_list)

        # ---- 4. 创建事件轨道 ----
        eventTrack <- maser:::createAnnotationTrack_event(eventGr, type)
        
        # ---- 5. 获取外显子编号（用于标注） ----
        gtf_exons <- gtf[gtf$type == "exon", ]
        
        # ---- 6. 提取 Inclusion 和 Skipping 的片段 ----
        if (type == "SE") {
            inclusion_gr <- eventGr[[grep("exon_target", names(eventGr))]]
            skipping_gr <- unlist(eventGr[grep("exon_upstream|exon_downstream", names(eventGr))])
        } else if (type == "RI") {
            inclusion_gr <- eventGr[[grep("exon_ir", names(eventGr))]]
            skipping_gr <- unlist(eventGr[grep("exon_upstream|exon_downstream", names(eventGr))])
        } else if (type %in% c("A3SS", "A5SS")) {
            inclusion_gr <- eventGr[[grep("exon_short", names(eventGr))]]
            skipping_gr <- eventGr[[grep("exon_long", names(eventGr))]]
        } else if (type == "MXE") {
            inclusion_gr <- eventGr[[grep("exon_1", names(eventGr))]]
            skipping_gr <- eventGr[[grep("exon_2", names(eventGr))]]
        }
        print("inclusion_gr:")
        print(inclusion_gr)
        print("skipping_gr:")
        print(skipping_gr)
        
        # ---- 7. 为片段添加外显子编号（从GTF中获取） ----
        if (show_exon_number && length(inclusion_gr) > 0) {
            # 获取基因的所有外显子用于编号
            gene_exons <- gtf_exons[gtf_exons$gene_name == gene_name]
            # 得到注释轨道基本信息
            strand <- unique(gene_exons$strand)
            Chr_id <- unique(gene_exons$seqnames)
            feature <- gene_name
            gene_exons <- sort(gene_exons)
            gene_exons$exon_label <- gene_exons$exon_id         # 指定外显子标记列
            # ---- 7.1 为 inclusion 片段添加编号 ----
            # 通过坐标匹配找到 inclusion 片段对应的外显子编号
            inclusion_hits <- findOverlaps(inclusion_gr, gene_exons)
            if (length(inclusion_hits) == 0) next
            # nearest_idx <- nearest(inclusion_gr, gene_exons)          # 找到相邻外显子
            inclusion_exon_labels <- gene_exons$exon_label[subjectHits(inclusion_hits)]
            # 如果多个片段重叠同一个外显子，取第一个（通常不会）
            inclusion_labels <- sapply(split(inclusion_exon_labels, queryHits(inclusion_hits)), `[`, 1)
            
            # ---- 7.2 为 skipping 片段添加编号 ----
            skipping_hits <- findOverlaps(skipping_gr, gene_exons)
            if (length(skipping_hits) == 0) next
            skipping_exon_labels <- gene_exons$exon_label[subjectHits(skipping_hits)]
            skipping_labels <- sapply(split(skipping_exon_labels, queryHits(skipping_hits)), `[`, 1)
            
            # ---- 7.3 创建带编号的 AnnotationTrack ----
            AnnotationTrack(
       range = NULL,
       start = NULL,
       end = NULL,
       width = NULL,
       feature,
       group,
       id,
       strand,
       chromosome,
       genome,
       stacking = "squish",
       name = "AnnotationTrack",
       fun,
       selectFun,
       importFunction,
       stream = FALSE,
       ...
     )

            inclusion_track <- AnnotationTrack(
                                    range = inclusion_gr,
                                    name = "Inclusion",
                                    id = inclusion_labels,
                                    showFeatureId = TRUE,
                                    cex = 0.6,
                                    fill = "#E69F00",
                                    col = "black",
                                    stacking = "dense")
                    
            skipping_track <- AnnotationTrack(
                                    range = skipping_gr,
                                    name = "Skipping",
                                    id = skipping_labels,
                                    showFeatureId = TRUE,
                                    fill = "#CC79A7",
                                    stacking = "dense")
            
            txnTracks <- list(
                    inclusionTrack = inclusion_track,
                    skippingTrack = skipping_track)
        } else {
            # 不显示外显子编号，使用原始方法
            txnTracks <- maser:::createAnnotationTrack_transcripts(eventGr, gtf_exons, type, is_strict)
        }
        
        # ---- 8. PSI轨道 ----
        if (show_PSI) {
            PSI <- slot(events, paste0(type, "_PSI"))
            PSI_event <- PSI[idx.event, , drop = FALSE]
            groups <- factor(c(rep(events@conditions[1], events@n_cond1),
                            rep(events@conditions[2], events@n_cond2)), 
                            levels = events@conditions)
            psiTrack <- maser:::createPSITrack_event(eventGr, PSI_event, groups, type, zoom)
            trackList <- list(psiTrack, eventTrack, txnTracks$inclusionTrack, txnTracks$skippingTrack)
        } else {
            trackList <- list(eventTrack, txnTracks$inclusionTrack, txnTracks$skippingTrack)
        }
        
        # ---- 9. 绘图 ----
        if (zoom) {
            Gviz::plotTracks(trackList, col.line = NULL, col = NULL,
                            Inclusion = "orange", Skipping = "purple", 
                            Retention = "orange", Non_Retention = "purple",
                            MXE_Exon1 = "orange", MXE_Exon2 = "purple",
                            A5SS_Short = "orange", A5SS_Long = "purple", 
                            A3SS_Short = "orange", A3SS_Long = "purple",
                            from = start(range(unlist(eventGr))) - 500, 
                            to = end(range(unlist(eventGr))) + 500)
        } else {
            Gviz::plotTracks(trackList, col.line = NULL, col = NULL,
                            Inclusion = "orange", Skipping = "purple", 
                            Retention = "orange", Non_Retention = "purple",
                            MXE_Exon1 = "orange", MXE_Exon2 = "purple",
                            A5SS_Short = "orange", A5SS_Long = "purple", 
                            A3SS_Short = "orange", A3SS_Long = "purple")
        }
    }

##3) 生成sashimi输入 [兼容改造: 从summary_table重建]
    generate_sashimi_input <- function(cache, gene_name, outdir, verbose = TRUE) {
        if (verbose) message("  生成 Sashimi 输入: ", gene_name)
        
        valid_genes <- names(cache$by_gene[["summary"]])
        if (!gene_name %in% valid_genes) {
            stop("基因 ", gene_name, " 不在缓存中")
        }
        gene_summary <- cache$by_gene[["summary"]][[gene_name]]

        if (nrow(gene_summary) == 0) {
            stop("基因 ", gene_name, " 无剪接事件")
        }
        
        event_types <- unique(gene_summary$EvenType)
        
        output_dir <- file.path(outdir, gene_name, "Sashimi_input")
        if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }
        
        event_summary <- data.frame()
        
        for (evt in event_types) {
            evt_data <- gene_summary[gene_summary$EvenType == evt, ]
            if (nrow(evt_data) == 0) next
            
            mats_file <- file.path(output_dir, paste0(evt, ".MATS.txt"))
            write.table(evt_data, file = mats_file, sep = "\t", quote = FALSE, row.names = FALSE)
            
            event_summary <- base::rbind(event_summary, data.frame(
                event_type = evt,
                event_count = nrow(evt_data),
                mats_file = mats_file,
                stringsAsFactors = FALSE
            ))
        }
        
        config_file <- file.path(output_dir, "sashimi_config.txt")
        config_lines <- c(
            paste0("# Sashimi 配置文件"),
            paste0("# 基因: ", gene_name),
            paste0("# 生成时间: ", Sys.time()),
            paste0("# 事件总数: ", nrow(gene_summary)),
            "",
            paste0("GENE = ", gene_name),
            paste0("EVENT_TYPES = ", paste(event_types, collapse = ","))
        )
        writeLines(config_lines, config_file)
        
        return(list(
            gene = gene_name,
            output_dir = output_dir,
            event_types = event_types,
            event_summary = event_summary,
            config_file = config_file,
            status = "success"
        ))
    }
##4) GRanges坐标转换bed12 [单基因处理]
    GR_to_bed12 <- function(gg,gr,outdir){
                        ## 建立文件夹
                        outdir1 <- file.path(outdir,gg)
                        if(!file.exists(outdir1)){
                            dir.create(outdir1, recursive = TRUE)
                        }
                        cat("文件存储路径: ", outdir1,"\n")
                        # 占位添加thick位点
                        mcols(gr)$thick <- IRanges(start=1, end=1)
                        mcols(gr)$score <- 0            # 保证输出为bed12格式
                        mcols(gr)$name <- names(gr)
                        # 根据基因提取显著剪接事件
                        gr_gg <- gr[which(gr$geneSymbol == gg),]
                        # 依据剪接类型拆分
                        gr_gg_list <- split(gr_gg,gr_gg$eventLabel)
                        # 直接调用export.bed12()转换格式
                        export(gr_gg_list, 
                               con = file.path(outdir1, "All_Splice_Events_IGV.bed"), 
                               format = "BED")         # 告诉它厚块在哪里

                        cat("✅ 成功导出 BED 文件，请用 IGV 加载\n")}

##5) 批量基因特异性分析 [兼容改造]
    batch_gene_specific <- function(cache,
                                    maser_obj = NULL,
                                    gene_list,
                                    outdir_gene,
                                    gtf_path = NULL,
                                    event_types = c("SE","A3SS","A5SS","RI","MXE"),
                                    verbose = TRUE,
                                    generate_sashimi = FALSE) {
                genes <- unlist(strsplit(gene_list, ","))
                genes <- trimws(genes)
                
                # 检查基因是否为显著差异剪接分子
                valid_genes <- c(cache$genes$Up,cache$genes$Down)
                if (length(setdiff(genes, valid_genes)) > 0) cat(setdiff(genes, valid_genes)," Not Found!\n")
                genes <- intersect(genes, valid_genes)

                if (length(genes) == 0) {
                    stop("  无有效基因可分析")
                }
                cat("\n ---- Gene Specific Processing ",genes," ---- \n")
                
                if (!is.null(gtf_path)) {
                    if (verbose) message("  读取GTF文件...")
                    ens_gtf <- rtracklayer::import.gff(gtf_path)
                } else {
                    ens_gtf <- NULL
                }
                
                if (verbose) message("\n===== 开始批量基因特异性分析 (", length(genes), " 个基因) =====")
                
                results <- list()
                for (gene in genes) {
                    result <- generate_gene_specific(
                                        cache = cache,
                                        maser_obj = maser_obj,
                                        gene_name = gene,
                                        ens_gtf = ens_gtf,
                                        outdir = outdir_gene,
                                        event_types = event_types,
                                        verbose = verbose)
                    
                    if (generate_sashimi) {
                        sashimi_result <- generate_sashimi_input(
                                                            cache = cache,
                                                            gene_name = gene,
                                                            outdir = outdir_gene,
                                                            verbose = verbose)
                        result$sashimi <- sashimi_result
                    }
                    results[[gene]] <- result
                }
                
                if (verbose) {
                    message("\n===== 批量分析完成 =====")
                    message("  基因数: ", length(genes))
                }
                
                return(results)
            }

##### =======================================================