#!/usr/bin/env Rscript

# ============================================================
# 运行环境: mamba activate r452
# 剪接事件坐标 → BED12 格式转换器
# 支持格式: chr:start-end:start-end:... 或 chr:strand:start-end:...
# ============================================================

library(stringr)

#' 将坐标字符串转换为BED12格式
#' 
#' @param coord_str 坐标字符串，如 "chrX:-:50573008-50573144:50569277-50570291:50573224-50573351"
#' @param event_name 事件名称（可选）
#' @return BED12格式的数据框行
#'

  # 测试参数 ===============================
  # AStype <- "MXE"         # 结构 exon1:exon2:flanking1:flanking2
  # event_name <- "MXE:Firre_14051"
  # coord_str <- "chrX:-:50606076-50606162:50608471-50608814:50600744-50600898:50615361-50615515"

  # AStype <- "SE"          # 结构 exon:flanking1:flanking2
  # event_name <- "SE:Firre_72368"
  # coord_str <- "chrX:-:50573008-50573144:50569277-50570291:50573224-50573351"

  # AStype <- "A5SS"        # 结构 long_exon:short_exon:flanking
  # event_name <- "A5SS:Mybpc3_30099"
  # coord_str <- "chr2:+:91121856-91122342:91121856-91121955:91122636-91122647"
  
  # AStype <- "A3SS"        # 结构 long_exon:short_exon:flanking
  # event_name <- "A3SS:Mybpc3_40811"
  # coord_str <- "chr2:+:91122197-91122646:91122639-91122646:91121902-91121955"

  # AStype <- "RI"          # 结构 target:flanking1:flanking2
  # event_name <- "RI:Myh6_6606"
  # coord_str <- "chr14:-:54946313-54946682:54946313-54946516:54946635-54946682"
### =========================================================

# 扩展原有的coord_to_bed12函数，支持事件类型参数
coord_to_bed12 <- function(coord_str, event_name = "event", event_type = "SE") {
    # 解析字符串
    parts <- strsplit(coord_str, ":")[[1]]
    
    # 提取染色体和链
    chrom <- parts[1]
    if (parts[2] %in% c("+", "-")) {
      strand <- parts[2]
      region_parts <- parts[3:length(parts)]
    } else {
      strand <- "+"
      region_parts <- parts[2:length(parts)]
    }
    
    # 根据事件类型解析区域
    if (event_type == "SE") {
      # SE: 3个区域 (upstream, target, downstream)
      regions <- parse_regions(region_parts[1:3])
      # 按基因组位置排序（确保正确的顺序）
      regions <- regions[order(sapply(regions, function(r) r$bed_start))]
      
    } else if (event_type %in% c("A3SS", "A5SS")) {
        # A3SS/A5SS: long和short是同一外显子的不同版本
        # 需要生成两个独立的BED12记录
        # 这里我们为单个事件生成两个转录本
        # 注意：long和short在fingerprint中顺序可能是 long:short:flanking
        flanking <- parse_regions(region_parts[3])[[1]]
        exon_long <- parse_regions(region_parts[1])[[1]]
        exon_short <- parse_regions(region_parts[2])[[1]]
        
        # 生成两个转录本
        regions_long <- list(flanking, exon_long)
        regions_long <- regions_long[order(sapply(regions_long, function(r) r$bed_start))]
        
        regions_short <- list(flanking, exon_short)
        regions_short <- regions_short[order(sapply(regions_short, function(r) r$bed_start))]
        
        # 为长转录本生成BED12
        bed_long <- create_bed12_from_regions(chrom, strand, regions_long, paste0(event_name, "_long"))
        # 为短转录本生成BED12
        bed_short <- create_bed12_from_regions(chrom, strand, regions_short, paste0(event_name, "_short"))
        
        return(list(long = bed_long, short = bed_short))

    } else if (event_type == "RI") {
      # RI: 3个区域 (upstream, intron_retention, downstream)
      regions <- parse_regions(region_parts[1:3])
      regions <- regions[order(sapply(regions, function(r) r$bed_start))]
      
    } else if (event_type == "MXE") {
      # MXE: 4个区域 (upstream, exon1, exon2, downstream)
      # 生成两个独立的BED12记录
      upstream <- parse_regions(region_parts[3])[[1]]
      exon1 <- parse_regions(region_parts[1])[[1]]
      exon2 <- parse_regions(region_parts[2])[[1]]
      downstream <- parse_regions(region_parts[4])[[1]]
      
      # 转录本1: upstream + exon1 + downstream
      regions1 <- list(upstream, exon1, downstream)
      regions1 <- regions1[order(sapply(regions1, function(r) r$bed_start))]
      
      # 转录本2: upstream + exon2 + downstream
      regions2 <- list(upstream, exon2, downstream)
      regions2 <- regions2[order(sapply(regions2, function(r) r$bed_start))]
      
      bed1 <- create_bed12_from_regions(chrom, strand, regions1, 
                                        paste0(event_name, "_exon1"))
      bed2 <- create_bed12_from_regions(chrom, strand, regions2, 
                                        paste0(event_name, "_exon2"))
      
      return(list(tx1 = bed1, tx2 = bed2))
    }
    
    # 对于SE和RI，直接生成BED12
    return(create_bed12_from_regions(chrom, strand, regions, event_name))
}

# 辅助函数：解析区域字符串
parse_regions <- function(region_strings) {
    lapply(region_strings, function(p) {
      coords <- strsplit(p, "-")[[1]]
      if (length(coords) != 2) {
        stop(paste("区域格式错误:", p))
      }
      start <- as.numeric(coords[1])
      end <- as.numeric(coords[2])
      list(start = start, 
          end = end,
          bed_start = start - 1,
          bed_end = end,
          length = end - start + 1)
    })
}

# 辅助函数：从regions列表创建BED12
create_bed12_from_regions <- function(chrom, strand, regions, name) {
    chrom_start <- min(sapply(regions, function(r) r$bed_start))
    chrom_end <- max(sapply(regions, function(r) r$bed_end))
    
    block_count <- length(regions)
    block_sizes <- paste(sapply(regions, function(r) r$length), collapse = ",")
    block_starts <- paste(sapply(regions, function(r) r$bed_start - chrom_start), collapse = ",")
    
    data.frame(
      chrom = chrom,
      chromStart = chrom_start,
      chromEnd = chrom_end,
      name = name,
      score = 0,
      strand = strand,
      thickStart = chrom_start,
      thickEnd = chrom_end,
      itemRgb = 0,
      blockCount = block_count,
      blockSizes = paste0(block_sizes, ","),
      blockStarts = paste0(block_starts, ","),
      stringsAsFactors = FALSE
    )
}

# 批量处理函数
fingerprint_to_bed12 <- function(evt_type, cache, gene, event_id) {
    # 生成fingerprint
    fps <- cache$
    
    # 根据事件类型处理
    if (evt_type == "SE" || evt_type == "RI") {
      # 单个转录本，直接合并
      bed_list <- lapply(names(fps), function(name) {
        coord_to_bed12(fps[name], name, evt_type)
      })
      return(do.call(rbind, bed_list))
      
    } else if (evt_type %in% c("A3SS", "A5SS")) {
      # 两个转录本，分别收集
      long_list <- list()
      short_list <- list()
      
      for (name in names(fps)) {
        result <- coord_to_bed12(fps[name], name, evt_type)
        long_list[[name]] <- result$long
        short_list[[name]] <- result$short
      }
      
      # 合并所有长转录本和短转录本
      bed_long <- do.call(rbind, long_list)
      bed_short <- do.call(rbind, short_list)
      
      # 添加事件类型标记
      bed_long$event_type <- evt_type
      bed_short$event_type <- evt_type
      bed_long$isoform <- "long"
      bed_short$isoform <- "short"
      
      return(rbind(bed_long, bed_short))
      
    } else if (evt_type == "MXE") {
      # 两个互斥外显子转录本
      tx1_list <- list()
      tx2_list <- list()
      
      for (name in names(fps)) {
        result <- coord_to_bed12(fps[name], name, evt_type)
        tx1_list[[name]] <- result$tx1
        tx2_list[[name]] <- result$tx2
      }
      
      bed_tx1 <- do.call(rbind, tx1_list)
      bed_tx2 <- do.call(rbind, tx2_list)
      
      bed_tx1$event_type <- evt_type
      bed_tx2$event_type <- evt_type
      bed_tx1$isoform <- "exon1"
      bed_tx2$isoform <- "exon2"
      
      return(rbind(bed_tx1, bed_tx2))
    }
}