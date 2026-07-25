## 或者使用pacman::p_load()

## processx

## 获得所有安装包的版本
#dd <- installed.packages() %>% as.data.frame() %>% select(c("Package","Version"))

## 本脚本涉及环境
# BasicR = 4.4.3      ## 不建议在这个环境运行脚本
# irGSEA = 4.5.3
# monocle3_R4.5.2 = 4.5.2

# 第一步：先加载Bioconductor核心包（单独加载，不嵌套suppressMessages）
# 先指定库路径，避免多路径冲突
irGSEA_lib <- "/data/med-hancs/apps/anaconda3/2022.10/envs/irGSEA/lib/R/library"
BasicR_lib <- "/data/med-hancs/apps/anaconda3/2022.10/envs/BasicR/lib/R/library"
Monocle_lib <- "/data/med-hancs/apps/anaconda3/2022.10/envs/monocle3_R4.5.2/lib/R/library"
other_lib <- "/data/med-wangcq/01CondaEnv/00DataBase/00Tools/seeksoultools.1.2.0/lib/R/library"
.libPaths(c(.libPaths(), BasicR_lib, Monocle_lib, irGSEA_lib, other_lib)) # 统一库路径优先级

### 运行脚本前必须先加载系统模块
# module load hdf5/1.10.4-gcc-4.8.5       ## 主要是Azimuth需要

# 这里添加检测运行环境
env_nm <- Sys.getenv("CONDA_DEFAULT_ENV")
cat(paste0("当前运行环境: ",env_nm, "\n"))

# 这里添加检测R版本
#R_ver <- paste(R.version$major,R.version$minor,sep=".")
R_ver <- getRversion()
cat(paste0("当前R版本: ",R_ver, "\n"))

# 获取当前运行环境的安装包需要最低R版本 [这里可以再改进一点，获取所有环境的R包版本，分sheet记录]
pkgs <- installed.packages()
df <- data.frame(包名 = rownames(pkgs),
                 版本 = pkgs[, "Version"],
                 依赖R版本 = pkgs[, "Depends"],
                 建议R版本 = pkgs[, "Suggests"],
                 内置路径 = pkgs[, "LibPath"])
# 导出表格
write.csv(df, paste0("环境",env_nm,"_R包版本依赖清单.csv"), row.names = F)

pkgs <- c(
  # 基础核心
  "methods","utils","conflicted","data.table","openxlsx","writexl","tibble","readxl",
  # 数据处理
  "dplyr","tidyr","stringr","proxy","ggpubr",
  # 剪接处理
  "maser","rtracklayer","qs","GenomicFeatures","txdbmaker",
  # 绘图可视化
  "ggplot2","ggvenn","ggrepel","pheatmap","igraph","karyoploteR","Gviz","GenomicRanges","enrichplot","ComplexHeatmap","ggtangle",
  # 统计分析
  "stats","WGCNA","DESeq2",
  # 注释&富集
  "clusterProfiler","pairedGSEA","ReactomePA","org.Mm.eg.db","org.Hs.eg.db","simplifyEnrichment",
  "ReactomePA"# "org.Rn.eg.db"  # 大鼠（如需要取消注释）
)

# 静默批量加载
invisible(suppressMessages(lapply(pkgs, library, character.only=TRUE, quietly=TRUE)))

# ===================== 4. 环境专用包（最后加载） =====================
if(env_nm == "monocle3_R4.5.2") library(monocle3, quietly=TRUE)
if(env_nm == "irGSEA") library(irGSEA, quietly=TRUE)

message("✅ 所有包加载完成，无任何依赖冲突！")

## 指定函数调用包
conflict_prefer("select", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("intersect", "base")
conflicts_prefer(WGCNA::cor)
conflicts_prefer(stats::dist)
conflicts_prefer(stats::as.dist)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::first)
conflicts_prefer(igraph::union)
conflicts_prefer(GenomicRanges::union)
conflicts_prefer(generics::union)
conflicts_prefer(base::union)
conflicts_prefer(base::setdiff)