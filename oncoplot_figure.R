# oncoplot
# whole exome seq analysis 

library(data.table)
library(maftools)

# -------- READ MUTATION TABLE --------

# Tab-delimited table of annotated somatic coding variants identified from whole-exome sequencing samples. 
# Each row represents a single variant and includes sample ID, gene symbol, functional consequence annotation, genomic coordinates, reference/alternate alleles, and variant type. 
# This table served as the input for mutation frequency analyses and oncoplot generation.

mut <- fread("oncoplot_mutations_all_raw.tsv")

# -------- BUILD CLINICAL DATA --------

clinical <- data.frame(
  Tumor_Sample_Barcode = unname(rename_map),

  Self_Renewal = c(
    "LowSR", "LowSR", "LowSR", "LowSR", "LowSR", "LowSR",
    "HighSR", "HighSR", "HighSR", "HighSR"
  ),
  
  Smoking_Status = c(
    "Former", "Former", "Current", "Former", "Current", "Former",
    "Current", "Former", "Former", "Former"
  ),
  
  stringsAsFactors = FALSE
)

write.table(
  clinical,
  file = "clinical_data.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

clinical

# -------- CREATE MAF DATAFRAME --------

maf_df <- data.frame(
  Hugo_Symbol = mut$gene,
  Tumor_Sample_Barcode = mut$sample,
  Variant_Classification = mut$effect,
  Variant_Type = mut$variant_type,
  Chromosome = mut$chrom,
  Start_Position = mut$pos,
  End_Position = mut$pos,
  Reference_Allele = mut$ref,
  Tumor_Seq_Allele2 = mut$alt,
  stringsAsFactors = FALSE
)

# -------- CLEAN MUTATION CLASS NAMES --------

maf_df$Variant_Classification[
  grepl("missense", maf_df$Variant_Classification)
] <- "Missense_Mutation"

maf_df$Variant_Classification[
  grepl("stop_gained", maf_df$Variant_Classification)
] <- "Nonsense_Mutation"

maf_df$Variant_Classification[
  grepl("splice", maf_df$Variant_Classification)
] <- "Splice_Site"

maf_df$Variant_Classification[
  grepl("frameshift", maf_df$Variant_Classification) &
    maf_df$Variant_Type == "DEL"
] <- "Frame_Shift_Del"

maf_df$Variant_Classification[
  grepl("frameshift", maf_df$Variant_Classification) &
    maf_df$Variant_Type == "INS"
] <- "Frame_Shift_Ins"

maf_df$Variant_Classification[
  grepl("inframe_deletion", maf_df$Variant_Classification)
] <- "In_Frame_Del"

maf_df$Variant_Classification[
  grepl("inframe_insertion", maf_df$Variant_Classification)
] <- "In_Frame_Ins"

maf_df$Variant_Classification[
  grepl("start_lost|initiator_codon_variant", maf_df$Variant_Classification)
] <- "Translation_Start_Site"

maf_df$Variant_Classification[
  grepl("stop_lost", maf_df$Variant_Classification)
] <- "Nonstop_Mutation"

maf_df <- maf_df[
  !grepl("stop_retained_variant|start_retained_variant",
         maf_df$Variant_Classification),
]

# -------- FILTER SHARED VARIANTS --------

maf_df$variant_id <- paste(
  maf_df$Hugo_Symbol,
  maf_df$Chromosome,
  maf_df$Start_Position,
  maf_df$Reference_Allele,
  maf_df$Tumor_Seq_Allele2,
  sep = "_"
)

variant_sample_counts <- aggregate(
  Tumor_Sample_Barcode ~ variant_id,
  data = maf_df,
  FUN = function(x) length(unique(x))
)

shared_common <- variant_sample_counts$variant_id[
  variant_sample_counts$Tumor_Sample_Barcode >= 8
]

maf_df_filtered <- maf_df[
  !(maf_df$variant_id %in% shared_common),
]

maf_df_filtered$variant_id <- NULL

# -------- FILTER TO LUNG SQUAMOUS / LUNG CANCER GENES --------
# From CBioPortal

lung_genes <- c(
  
  # classic LUSC tumor suppressors / oncogenes
  "TP53",
  "PIK3CA",
  "CDKN2A",
  "NFE2L2",
  "KEAP1",
  "SOX2",
  "FAT1",
  "KMT2D",
  "NOTCH1",
  "FBXW7",
  "PTEN",
  "ARID1A",
  "RB1",
  "EP300",
  "CREBBP",
  "NF1",
  
  # RTK/RAS signaling
  "KRAS",
  "EGFR",
  "ERBB2",
  "BRAF",
  "MET",
  "FGFR1",
  "FGFR2",
  "FGFR3",
  "MAP2K1",
  "NRAS",
  "HRAS",
  
  # chromatin / epigenetic
  "KDM6A",
  "KDM5C",
  "SETD2",
  "SMARCA4",
  "SMARCB1",
  "PBRM1",
  "BAP1",
  
  # DNA damage / repair
  "ATM",
  "ATR",
  "BRCA1",
  "BRCA2",
  "CHEK2",
  
  # squamous differentiation
  "TP63",
  "ZNF750",
  "IRF6",
  "AJUBA",
  
  # oxidative stress / signaling
  "CUL3",
  "STK11",
  "RICTOR",
  "AKT1",
  
  # cell cycle
  "CCND1",
  "CDK4",
  "CDK6",
  
  # additional recurrent lung genes
  "NAV3",
  "DDR2",
  "MYC",
  "MCL1",
  "TERT"
)

maf_df_lung <- maf_df_filtered[
  maf_df_filtered$Hugo_Symbol %in% lung_genes,
]

# check remaining genes
sort(unique(maf_df_lung$Hugo_Symbol))

# -------- WRITE FILTERED MAF --------

write.table(
  maf_df_lung,
  file = "lung_driver_filtered.maf",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -------- READ FILTERED MAF WITH CLINICAL DATA --------

maf <- read.maf(
  maf = "lung_driver_filtered.maf",
  clinicalData = "clinical_data.tsv"
)

getClinicalData(maf)
getGeneSummary(maf)

# -------- ANNOTATION COLORS --------

ann_colors <- list(
  Smoking_Status = c(
    "Current" = "red",
    "Former" = "orange"
  ),
  
  Self_Renewal = c(
    "HighSR" = "blue",
    "LowSR" = "pink"
  )
)

# -------- CREATE LUNG DRIVER ONCOPLOT --------

while (!is.null(dev.list())) dev.off()

pdf(
  "lung_driver_oncoplot.pdf",
  width = 10,
  height = 7
)

oncoplot(
  maf = maf,
  top = 25,
  titleText = "Lung squamous cancer gene mutations",
  clinicalFeatures = c("Self_Renewal", "Smoking_Status"),
  annotationColor = ann_colors,
  sortByAnnotation = TRUE,
  annotationOrder = c("HighSR", "LowSR"),
  groupAnnotationBySize = FALSE,
  removeNonMutated = TRUE
)

dev.off()

getGeneSummary(maf)[1:30, ]

table(maf_df$Variant_Classification)


