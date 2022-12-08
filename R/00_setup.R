# Setup doc for lists and important vectors

library(tidyverse)
library(BiocManager)
library(biomaRt)
library(hypeR)
library(NatParksPalettes)
library(DESeq2)
library(purrr)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(ggpubr)
library(ggpubr)
library(ggsignif)
library(ggplot2)
library(ComplexHeatmap)
library(tidyverse)
library(foreach)
library(doParallel)
library(TCGAbiolinks)
library(DEXSeq)
library(beepr)
library(ggbiplot)
library(M3C)
library(circlize)
library(plotrix)
library(Rtsne)
library(beepr)
cl <- makePSOCKcluster(detectCores() - 2)
registerDoParallel(cl)

files <- list.files("/projectnb2/evolution/zwakefield/tcga/runs/")

sample_sheet <- read.delim("/projectnb2/evolution/zwakefield/tcga/metadata/sample_sheet_full_metadata.txt", sep = "\t") %>% filter(File.ID %in% files)

sample_sheet <- sample_sheet[!duplicated(sample_sheet$File.ID),]

sample_sheet <- transform(sample_sheet,id=as.numeric(factor(Project.ID)))

sample_sheet$id[sample_sheet$Sample.Type == "Solid Tissue Normal"] <- 0

sample_sheet$bin_id <- 1

sample_sheet$bin_id[sample_sheet$Sample.Type == "Solid Tissue Normal"] <- 0

matchedSample_sheet <- sample_sheet[sample_sheet$Case.ID %in% names(table(sample_sheet$Case.ID)[table(sample_sheet$Case.ID) >= 2]),] %>% arrange(Case.ID)

matchedSample_sheet <- matchedSample_sheet[matchedSample_sheet$Case.ID %in% unique(matchedSample_sheet$Case.ID)[unlist(lapply(unique(matchedSample_sheet$Case.ID), 
                           function(x) sum(c("Primary Tumor", "Solid Tissue Normal") %in% matchedSample_sheet$Sample.Type[matchedSample_sheet$Case.ID == x]) == 2))],]

kallisto <-paste("/projectnb2/evolution/zwakefield/tcga/runs/", files, "/kallisto/abundance.tsv", sep ="")

se <- paste("/projectnb2/evolution/zwakefield/tcga/runs/", files, "/rmats/SE.MATS.JC.txt", sep = "")

ri <- paste("/projectnb2/evolution/zwakefield/tcga/runs/", files, "/rmats/RI.MATS.JC.txt", sep = "")

mxe <- paste("/projectnb2/evolution/zwakefield/tcga/runs/", files, "/rmats/MXE.MATS.JC.txt", sep = "")

a5ss <- paste("/projectnb2/evolution/zwakefield/tcga/runs/", files, "/rmats/A5SS.MATS.JC.txt", sep = "")

a3ss <- paste("/projectnb2/evolution/zwakefield/tcga/runs/", files, "/rmats/A3SS.MATS.JC.txt", sep = "")

afe <- paste("/projectnb2/evolution/zwakefield/tcga/runs/", files, "/hit/", paste(files, ".AFEPSI", sep = ""), sep = "")

ale <- paste("/projectnb2/evolution/zwakefield/tcga/runs/", files, "/hit/", paste(files, ".ALEPSI", sep = ""), sep = "")

hit <- paste("/projectnb2/evolution/zwakefield/tcga/runs/", files, "/hit/", paste(files, ".exon", sep = ""), sep = "")

# usable genes derived from biolinks

genes <- foreach (i = 1:length(list.files("/projectnb/evolution/zwakefield/tcga/procedural/all_samples/dfs")[grep("biolinks", list.files("/projectnb/evolution/zwakefield/tcga/procedural/all_samples/dfs"))])) %dopar% {
  
  df <- read.delim(paste("/projectnb/evolution/zwakefield/tcga/procedural/all_samples/dfs/", list.files("/projectnb/evolution/zwakefield/tcga/procedural/all_samples/dfs")[grep("biolinks", list.files("/projectnb/evolution/zwakefield/tcga/procedural/all_samples/dfs"))][i], sep = ""), sep = " ")
  if (!is.na(df)) {
    df$gene[as.numeric(apply(df[,!(colnames(df) == "gene")], 1, function(x) sum(x == 0)/length(x))) < .90]
  }
  
}
substanceGenes <- Reduce(intersect, genes[unlist(lapply(genes, length)) != 0])


# Functions

extract_Biolinks <- function(type, fsample_sheet = sample_sheet) {
  q <- GDCquery(
    project = type,
    data.category = "Transcriptome Profiling", 
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  GDCdownload(q)
  q_d <- GDCprepare(q)
  data2 <- q_d
  q_d2 <- data.frame(assay(q_d))
  geneWiseCounts <- apply(q_d2, 1, sum)
  CountGreater0 <- geneWiseCounts > 0
  q_d2 <- q_d2[ CountGreater0, ]
  q_d2 <- log((q_d2) + 1, 2)
  rownames(q_d2) <- unlist(lapply(strsplit(rownames(q_d2), split = "[.]"), "[[", 1))
  q_d2$gene <- rownames(q_d2)
  q_d2 <- q_d2 %>% dplyr::relocate(gene)
  
  qc_p <- q_d2[,c(TRUE, q_d$sample %in% fsample_sheet$Sample.ID[fsample_sheet$File.ID %in% unlist(lapply(strsplit(kallisto[!is.na(kallisto)], split = '[/]'), "[[", 7))])]
  currType <- fsample_sheet[fsample_sheet$Project.ID %in% type,] %>% arrange(Case.ID)
  return(list(gex = qc_p,
              currType = currType))
  
}

extract_Kallisto <- function(type, fKal.list = kallisto, fsample_sheet = sample_sheet) {
  currType <- fsample_sheet[fsample_sheet$Project.ID %in% type,] %>% arrange(Case.ID)
  kal.vals <- foreach(i = 1:length(currType$File.ID)) %dopar% {
    
    library(tidyverse)
    library(ComplexHeatmap)
    library(reshape2)
    library(purrr)
    
    temp <- read.delim(fKal.list[grep(currType$File.ID[i], fKal.list)[1]], sep = "\t")
    temp <- temp[temp$tpm > 0,]
    temp$gene <- unlist(lapply(strsplit(unlist(lapply(strsplit(temp$target_id, split = "[|]"), "[[", 2)), split = "[.]"), "[[", 1))
    tempm <- temp %>% dplyr::select(gene, tpm)
    tempa <- aggregate(tpm~gene, tempm, sum)
    colnames(tempa)[2] <- paste(colnames(tempa)[2], ".", currType$File.ID[i], sep="")
    tempa
  }
  
  kal.t <- suppressWarnings(lapply(kal.vals, function(x) {if (!is.na(x)) {
    x
  }
  }) %>% purrr::reduce(full_join, by=c('gene')))
  
  kal.t[is.na(kal.t)] <- 0
  kal.t[,2:length(colnames(kal.t))] <- log((kal.t[,2:length(colnames(kal.t))])+1, 2)
  return(list(kal.t = kal.t,
              currType = currType))
}

extract_rMATS <- function(type, ex_t, frM.list, fsample_sheet = sample_sheet) {
  currType <- fsample_sheet[fsample_sheet$Project.ID %in% type,] %>% arrange(Case.ID)
  all.names <- c()
  rM.vals <- foreach(i = 1:length(currType$File.ID)) %dopar% {
    library(tidyverse)
    library(ComplexHeatmap)
    library(reshape2)
    library(purrr)
    
    temp <- read.delim(frM.list[grep(currType$File.ID[i], frM.list)[1]])
    
    if (ex_t == "SE") {
      temp <- temp %>% dplyr::select('GeneID', "chr", "strand", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE", "IncLevel1")
      temp$id <- paste(temp$GeneID, temp$chr, temp$strand, temp$exonStart_0base, temp$exonEnd, temp$upstreamES,
                       temp$upstreamEE, temp$downstreamES, temp$downstreamEE,sep = ";")
      
    } else if (ex_t == "A3SS" | ex_t == "A5SS") {
      temp <- temp %>% dplyr::select('GeneID', "chr", "strand", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE", "IncLevel1", "IncLevel2")
      temp$id <- paste(temp$GeneID, temp$chr, temp$strand, temp$longExonStart_0base, temp$longExonEnd, temp$shortES,
                       temp$shortEE, temp$flankingES, temp$flankingEE,sep = ";")
      
    } else if (ex_t == "MXE") {
      temp <- temp %>% dplyr::select('GeneID', "chr", "strand", "X1stExonStart_0base", "X1stExonEnd", "X2ndExonStart_0base", "X2ndExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE", "IncLevel1", "IncLevel2")
      temp$id <- paste(temp$GeneID, temp$chr, temp$strand, temp$X1stExonStart_0base, temp$X1stExonEnd, temp$X2ndExonStart_0base, temp$X2ndExonEnd,temp$upstreamES,
                       temp$upstreamEE, temp$downstreamES, temp$downstreamEE,sep = ";")
      
    } else if (ex_t == "RI") {
      temp <- temp %>% dplyr::select('GeneID', "chr", "strand", "riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE", "IncLevel1", "IncLevel2")
      temp$id <- paste(temp$GeneID, temp$chr, temp$strand, temp$riExonStart_0base, temp$riExonEnd, temp$upstreamES,
                       temp$upstreamEE, temp$downstreamES, temp$downstreamEE,sep = ";")
    }
    
    if (!(length(temp$IncLevel1) == sum(is.na(temp$IncLevel1)))) {
      temp <- temp %>% dplyr::select(id, IncLevel1)
    } else {
      temp <- temp %>% dplyr::select(id, IncLevel2)
    }
    colnames(temp)[2] <- "IncLevel" #paste("IncLevel", ".", currType$File.ID[i], sep="")
    temp <- temp[!(is.na(temp$IncLevel)),]
    temp <- temp %>% arrange(id)
    temp <- temp[!duplicated(temp$id),]
    temp
    
  }
  
  
  
  
  ids <- unique(unlist(lapply(rM.vals, function(x) x$id)))
  inco <- data.frame(id = ids, IncLevel = NA)
  try2 <- foreach(i =1:length(rM.vals)) %dopar% {
    hm <- rM.vals[[i]]
    try <-rbind(hm, inco)
    try <- try[!(duplicated(try$id)),] %>% arrange(id)
    # progress(table(is.na(try2[[i]]$IncLevel))[[1]], as.numeric(table(is.na(hm$IncLevel))))
    colnames(try)[2] <- paste("IncLevel", ".", currType$File.ID[i], sep="")
    try
  }
  table(unlist(lapply(try2, function(x) setequal(try2[[1]]$id, x$id))))
  comb.df <- do.call(cbind, try2)
  rM.comb <- comb.df[c(1, seq(2, length(colnames(comb.df)), by = 2))]
  lapply(try2, function(x) dim(x))
  rM.comb[is.na(rM.comb)] <- 0
  
  return(list(rM.t = rM.comb,
              currType = currType))
}

extract_afe <- function(type, fAFEPSI.list = afe, fsample_sheet = sample_sheet) {
  if (length(type) > 1 ) {
    title_type <- "TCGA-PanCancer"
  } else {
    title_type <- type
  }
  currType <- fsample_sheet[fsample_sheet$Project.ID %in% type,] %>% arrange(Case.ID)
  library(reshape2)
  library(purrr)
  type.df <- foreach (i=1:length(currType$File.ID)) %dopar% {
    library(tidyverse)
    t.df <- read.delim(fAFEPSI.list[grep(currType$File.ID[i], fAFEPSI.list)[1]]) %>% dplyr::select("gene", "exon", "strand", "AFEPSI")
    colnames(t.df)[4] <- paste(colnames(t.df)[4], ".", currType$File.ID[i], sep="")
    t.df
  }
  type.cAFEPSI <- suppressWarnings(lapply(type.df, function(x) {if (!is.na(x)) {
    x
  }
  }) %>% purrr::reduce(full_join, by=c('gene', 'exon', 'strand')))
  type.cAFEPSI[is.na(type.cAFEPSI)] <- 0
  return(list(AFE = type.cAFEPSI,
              curr = currType))
}

extract_ale <- function(type, fALEPSI.list = ale, fsample_sheet = sample_sheet) {
  if (length(type) > 1 ) {
    title_type <- "TCGA-PanCancer"
  } else {
    title_type <- type
  }
  currType <- fsample_sheet[fsample_sheet$Project.ID %in% type,] %>% arrange(Case.ID)
  library(reshape2)
  library(purrr)
  type.df <- foreach (i=1:length(currType$File.ID)) %dopar% {
    library(tidyverse)
    library(ComplexHeatmap)
    t.df <- read.delim(fALEPSI.list[grep(currType$File.ID[i], fALEPSI.list)[1]]) %>% dplyr::select("gene", "exon", "strand", "ALEPSI")
    colnames(t.df)[4] <- paste(colnames(t.df)[4], ".", currType$File.ID[i], sep="")
    t.df
  }
  type.cALEPSI <- suppressWarnings(lapply(type.df, function(x) {if (!is.na(x)) {
    x
  }
  }) %>% purrr::reduce(full_join, by=c('gene', 'exon', 'strand')))
  type.cALEPSI[is.na(type.cALEPSI)] <- 0
  return(list(ale = type.cALEPSI,
              curr = currType))
}

extract_hit <- function(type, fexon.list = hit, fsample_sheet = sample_sheet) {
  if (length(type) > 1 ) {
    title_type <- "TCGA-PanCancer"
  } else {
    title_type <- type
  }
  currType <- fsample_sheet[fsample_sheet$Project.ID %in% type,] %>% arrange(Case.ID)
  library(reshape2)
  library(purrr)
  type.df <- foreach (i=1:length(currType$File.ID)) %dopar% {
    library(tidyverse)
    library(ComplexHeatmap)
    t.df <- read.delim(fexon.list[grep(currType$File.ID[i], fexon.list)[1]]) %>% dplyr::select("gene", "exon", "strand", "HITindex")
    colnames(t.df)[4] <- paste(colnames(t.df)[4], ".", currType$File.ID[i], sep="")
    t.df
  }
  type.cExon <- suppressWarnings(lapply(type.df, function(x) {if (!is.na(x)) {
    x
  }
  }) %>% purrr::reduce(full_join, by=c('gene', 'exon', 'strand')))
  type.cExon[is.na(type.cExon)] <- 0
  return(list(hit = type.cExon,
              curr = currType))
}


extract_id <- function(type, fexon.list = hit, fsample_sheet = sample_sheet) {
  if (length(type) > 1 ) {
    title_type <- "TCGA-PanCancer"
  } else {
    title_type <- type
  }
  currType <- fsample_sheet[fsample_sheet$Project.ID %in% type,] %>% arrange(Case.ID)
  
  type.df <- foreach (i=1:length(currType$File.ID)) %dopar% {
    library(reshape2)
    library(purrr)
    library(dplyr)
    library(tidyverse)
    library(ComplexHeatmap)
    t.df <- read.delim(fexon.list[grep(currType$File.ID[i], fexon.list)[1]]) %>% dplyr::select("gene", "exon", "strand", "ID")
    colnames(t.df)[4] <- paste(colnames(t.df)[4], ".", currType$File.ID[i], sep="")
    t.df
  }
  type.cID <- suppressWarnings(lapply(type.df, function(x) {if (!is.na(x)) {
    x
  }
  }) %>% purrr::reduce(full_join, by=c('gene', 'exon', 'strand')))
  return(list(id = type.cID,
              curr = currType))
}



## Colorwheel

n <- 33
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

n <- 33
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector2 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))