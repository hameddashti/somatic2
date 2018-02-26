{
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(SomaticSignatures)
  library(GenomicRanges)
  library(readr)
  library(ggplot2)
  library(dplyr)
  library(VariantAnnotation)
  library(tidyr)
  library(csv)
  library(dtplyr)
  library(readxl)
}

{
  b <- 'Bladder'
  
  inputfile <-
    paste0(
      '/home/san/halinejad/Desktop/Masroor Bayati/DeepCancer project/Data/ICGC/',
      b,
      '/simple_somatic_mutation.open.tsv'
    )
  
  Cancer <- data.frame(read_tsv(inputfile, col_names = TRUE))
  gene_list <-
    read_excel(
      '/home/san/halinejad/Desktop/Masroor Bayati/DeepCancer project/Data/annotation/FANTOM5_gene_list.xlsx'
    )
  annot <- as.data.frame(gene_list[, c(1, 3, 4, 5, 6, 8)])
  colnames(annot) <-
    c('id', 'chromosome', 'start', 'end', 'strand', 'geneClass')
  gene <- subset(annot, annot$geneClass == 'coding_mRNA')
  lncRNA <- subset(annot, annot$geneClass == 'lncRNA')
  pseudogene <- subset(annot, annot$geneClass == 'pseudogene')
  sense_overlap_RNA <-
    subset(annot, annot$geneClass == 'sense_overlap_RNA')
  short_ncRNA <- subset(annot, annot$geneClass == 'short_ncRNA')
  small_RNA <- subset(annot, annot$geneClass == 'small_RNA')
  structural_RNA <- subset(annot, annot$geneClass == 'structural_RNA')
  uncertain_coding <-
    subset(annot, annot$geneClass == 'uncertain_coding')
  tmp <- as.data.frame(Cancer[, c(1, 5, 9, 10, 11, 12, 16, 17)])
  tmp1 <- data.frame(subset(tmp, tmp$icgc_mutation_id != "NA"))
  tmp3 <-
    data.frame(
      subset(
        tmp1,
        nchar(tmp1$mutated_from_allele) == 1 &
          nchar(tmp1$mutated_to_allele) == 1
        &
          tmp1$mutated_from_allele != '-' & tmp1$mutated_to_allele != "-"
        &
          tmp1$mutated_from_allele != '_' & tmp1$mutated_to_allele != "_"
        &
          tmp1$chromosome != "MT" & tmp1$chromosome != "M"
      )
    )
  
  tmp2 <- tmp3[!duplicated(tmp3), ]
  tmp2 <-
    subset(tmp2, tmp2$mutated_from_allele != tmp2$mutated_to_allele)
  
  tmp2$chromosome_strand[which(tmp2$chromosome_strand == '1')] <-
    '+'
  tmp2$chromosome_strand[which(tmp2$chromosome_strand == '2')] <-
    '-'
  
  tmp2$chromosome <- paste0("chr", tmp2$chromosome)
  
  gr <- makeGRangesFromDataFrame(
    tmp2,
    keep.extra.columns = T,
    seqnames.field = 'chromosome',
    start.field = 'chromosome_start',
    end.field = 'chromosome_end',
    strand.field = 'chromosome_strand'
  )
  
  
  idx <-
    match(c(
      'mutated_from_allele',
      'mutated_to_allele',
      'icgc_sample_id'
    ),
    names(mcols(gr)))
  mcols(gr) <- cbind(mcols(gr)[idx], mcols(gr)[-idx])
  
  vr <- makeVRangesFromGRanges(
    gr,
    ref.field = 'mutated_from_allele',
    alt.field = 'mutated_to_allele',
    sampleNames.field = 'icgc_sample_id',
    keep.extra.columns = T
  )
  
  vr <- mutationContext(vr, Hsapiens)
  
  
  variations <-
    data.frame(
      icgc_samlpe_id = mcols(gr)$icgc_sample_id,
      chromosome = as.character(seqnames(vr)),
      position = start(vr),
      strand = as.character(strand(gr)),
      motif = paste0(
        as.character(mcols(vr)$alteration),
        '-',
        as.character(mcols(vr)$context)
      )
    )
}
x <- as.data.frame(variations[, c(1, 5, 2, 3)], header = T)
rm(tmp,
   Cancer,
   tmp1,
   tmp3,
   vr,
   variations,
   clear_var,
   gr,
   annot,
   tmp2,
   gene_list)
gc()

{
  df1 <- structure(x)
  df2 <- structure(as.data.frame(gene[, c(1, 2, 3, 4)]))
  q <- df1 %>% inner_join(df2, "chromosome") %>%
    mutate(
      geneID_motif = paste(id, motif, sep = ","),
      n = if_else(position >= start &
                    position <= end, 1, 0)
    ) %>%
    select(icgc_samlpe_id, geneID_motif, n) %>%
    group_by(icgc_samlpe_id, geneID_motif) %>%
    summarise(n = sum(n)) %>%
    spread(key = geneID_motif,
           value = n,
           fill = 0)
  q0 = q[, !grepl("N", colnames(q))]
  row_name <- q0$icgc_samlpe_id
  q0 <- q0[, -1]
  rownames(q0) <- row_name
  col_name <- colnames(q0)
  q0 <- colSums(q0)
  q0 <- as.data.frame(q0)
  k <- vector(max(q0), mode = 'list')
  ttt <- subset(q0, q0 >= 1)
  rn<-rownames(ttt)
  temporary<- q[,c(as.list(rn))]
  print(i)
  w <- paste0('/home/san/halinejad/Desktop/Masroor Bayati/DeepCancer project/Data/ICGC/',b, "_gene_bigger_than_", as.character(i),".csv")
  write.csv(temporary, file = w, row.names = TRUE)
}
