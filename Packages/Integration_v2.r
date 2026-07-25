#source("/data/med-wangcq/01CondaEnv/02Git_repo/00MyGit_wchenqi/Multi_omics/SplicingAnalysis/00Packages/01bulk/脚本调整/Gene_specific_v2.r")

##### 多组cache交集分析模块 [兼容改造: 从cache读取summary_table]
    compare_caches <- function(cache_files,
                                cache_labels = NULL,
                                gene_list = NULL,
                                gtf_path,
                                outdir = "./",
                                event_types = c("SE", "A5SS", "A3SS", "MXE", "RI"),
                                min_shared = 2,
                                verbose = TRUE) {
        
        if (length(cache_files) < 2) {
            stop("至少需要提供2个cache文件")
        }
        if (is.null(cache_labels)) {
            cache_labels <- paste0("Condition", seq_along(cache_files))
        }
        if (length(cache_labels) != length(cache_files)) {
            stop("cache_labels 长度必须与 cache_files 一致")
        }
        
        if (verbose) message("\n===== 加载cache文件 =====")
        caches <- list()
        for (i in seq_along(cache_files)) {
            if (verbose) message("  加载: ", cache_labels[i], " (", cache_files[i], ")")
            caches[[cache_labels[i]]] <- qread(cache_files[i])
        }
        
        if (verbose) message("\n===== 加载GTF文件 =====")
        ens_gtf <- rtracklayer::import.gff(gtf_path)
        
        outdir_intersect <- file.path(outdir, "Compare_ASevent")
        if (!dir.exists(outdir_intersect)) {
            dir.create(outdir_intersect, recursive = TRUE)
        }
        
        if (verbose) message("\n===== 分析交集事件 =====")
        
        all_results <- list()
        all_intersect_genes <- c()
        all_intersect_events <- list()
        
        for (evt in event_types) {
            if (verbose) message("\n  处理事件类型: ", evt)
            
            fp_maps <- list()
            for (label in cache_labels) {
                if (!evt %in% names(caches[[label]]$by_event_type)) {
                    if (verbose) message("    ", label, ": 无 ", evt, " 事件")
                    fp_maps[[label]] <- c()
                } else {
                    # 兼容改造: fingerprint结构不变
                    fp_maps[[label]] <- caches[[label]]$fingerprint[[evt]]
                    if (verbose) message("    ", label, ": ", length(fp_maps[[label]]), " 个事件")
                }
            }
            
            if (all(sapply(fp_maps, length) == 0)) {
                if (verbose) message("    所有条件均无 ", evt, " 事件，跳过")
                next
            }
            
            non_empty_maps <- fp_maps[sapply(fp_maps, length) > 0]
            if (length(non_empty_maps) < min_shared) {
                if (verbose) message("    有效映射数不足，跳过")
                next
            }
            
            common_fingerprints <- Reduce(intersect, non_empty_maps)
            if (verbose) message("    共有指纹数: ", length(common_fingerprints))
            
            if (length(common_fingerprints) == 0) next
            
            intersect_events <- list()
            for (label in names(non_empty_maps)) {
                map <- non_empty_maps[[label]]
                intersect_events[[label]] <- names(map[map %in% common_fingerprints])
            }
            
            intersect_genes <- c()
            for (label in names(intersect_events)) {
                genes <- unique(gsub("^[^:]+:([^_]+)_.*$", "\\1", intersect_events[[label]]))
                intersect_genes <- union(intersect_genes, genes)
            }
            
            if (!is.null(gene_list)) {
                genes_filter <- unique(trimws(unlist(strsplit(gene_list, ","))))
                if (verbose) message("    应用基因列表筛选 (", length(genes_filter), " 个基因)")

                filtered_events <- list()
                for (label in names(intersect_events)) {
                    evt_vec <- intersect_events[[label]]
                    keep <- sapply(evt_vec, function(x) {
                        gene <- strsplit(strsplit(x, ":")[[1]][2], "_")[[1]][1]
                        gene %in% genes_filter
                    })
                    filtered_events[[label]] <- evt_vec[keep]
                }
                intersect_events <- filtered_events

                intersect_genes <- c()
                for (label in names(intersect_events)) {
                    genes <- unique(gsub("^[^:]+:([^_]+)_.*$", "\\1", intersect_events[[label]]))
                    intersect_genes <- union(intersect_genes, genes)
                }

                if (length(unlist(intersect_events)) == 0) {
                    if (verbose) message("    基因列表筛选后无事件，跳过")
                    next
                }
            }

            event_df <- data.frame(
                eventLabel = unlist(intersect_events),
                Condition = rep(names(intersect_events), times = sapply(intersect_events, length)),
                stringsAsFactors = FALSE
            )
            write.csv(event_df, 
                    file.path(outdir_intersect, paste0("Intersect_", evt, "_events.csv")),
                    row.names = FALSE)
            
            gene_df <- data.frame(geneSymbol = intersect_genes)
            write.csv(gene_df,
                    file.path(outdir_intersect, paste0("Intersect_", evt, "_genes.csv")),
                    row.names = FALSE)
            
            all_results[[evt]] <- list(
                common_fingerprints = common_fingerprints,
                intersect_events = intersect_events,
                intersect_genes = intersect_genes,
                n_intersect_events = length(common_fingerprints),
                n_intersect_genes = length(intersect_genes)
            )
            
            all_intersect_genes <- union(all_intersect_genes, intersect_genes)
            all_intersect_events[[evt]] <- intersect_events
        }
        
        if (verbose) message("\n===== 交集分析汇总 =====")
        summary_df <- data.frame(
            EventType = names(all_results),
            IntersectEvents = sapply(all_results, function(x) x$n_intersect_events),
            IntersectGenes = sapply(all_results, function(x) x$n_intersect_genes)
        )
        print(summary_df)
        write.csv(summary_df,
                file.path(outdir_intersect, "Intersection_Summary.csv"),
                row.names = FALSE)
        
        all_genes_df <- data.frame(geneSymbol = all_intersect_genes)
        write.csv(all_genes_df,
                file.path(outdir_intersect, "All_Intersect_Genes.csv"),
                row.names = FALSE)
        
        if (verbose) message("\n===== 完成 =====")
        
        return(list(
            results_by_type = all_results,
            all_intersect_genes = all_intersect_genes,
            all_intersect_events = all_intersect_events,
            summary = summary_df,
            output_dir = outdir_intersect
        ))
    }