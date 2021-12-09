library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
sd =  999
set.seed(sd)

young <- read.table(
  file = "data/youngPeaks.txt",
  col.names = c("chr", "start", "end")
)
aged <- read.table(
  file = "data/agedPeaks.txt",
  col.names = c("chr", "start", "end")
)
young <- makeGRangesFromDataFrame(young)
aged <- makeGRangesFromDataFrame(aged)
combined.peaks <- reduce(x = c(young, aged))

md.young <- read.table(
  file = "GSM5723631_Young_HSC_singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row
md.aged <- read.table(
  file = "GSM5723632_Aged_HSC_singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

md.young <- md.young[md.young$passed_filters > 500, ]
md.aged <- md.aged[md.aged$passed_filters > 500, ]

frags.young <- CreateFragmentObject(
  path = "GSM5723631_Young_HSC_fragments.tsv.gz",
  cells = rownames(md.young)
)
frags.aged <- CreateFragmentObject(
  path = "GSM5723632_Aged_HSC_fragments.tsv.gz",
  cells = rownames(md.aged)
)

young.counts <- FeatureMatrix(
  fragments = frags.young,
  features = combined.peaks,
  cells = rownames(md.young)
)
aged.counts <- FeatureMatrix(
  fragments = frags.aged,
  features = combined.peaks,
  cells = rownames(md.aged)
)

young_assay <- CreateChromatinAssay(young.counts, fragments = frags.young)
young <- CreateSeuratObject(young_assay, assay = "ATAC", meta.data=md.young)
aged_assay <- CreateChromatinAssay(aged.counts, fragments = frags.aged)
aged <- CreateSeuratObject(aged_assay, assay = "ATAC", meta.data=md.aged)

young$dataset <- 'young'
aged$dataset <- 'aged'
combined <- merge(
  x = young,
  y = aged,
  add.cell.ids = c("young", "aged")
)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"
Annotation(combined) <- annotations
combined <- NucleosomeSignal(object = combined)
combined <- TSSEnrichment(object = combined, fast = FALSE)
combined$pct_reads_in_peaks <- combined$peak_region_fragments / combined$passed_filters * 100
combined$blacklist_ratio <- combined$blacklist_region_fragments / combined$peak_region_fragments
combined$high.tss <- ifelse(combined$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(combined, group.by = 'dataset') + NoLegend()
combined$nucleosome_group <- ifelse(combined$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
VlnPlot(
  object = combined,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
combined <- subset(
  x = combined,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:30, reduction = 'lsi')
combined <- FindNeighbors(object = combined, reduction = 'lsi', dims = 2:30)
combined <- FindClusters(object = combined, verbose = FALSE,resolution = 0.4, algorithm = 3)
DimPlot(combined, pt.size = 0.1)
DimPlot(combined, group.by = 'dataset', pt.size = 0.1)

lst <- read.table("data/DARopen_combinedRegion.txt")
lst <- list(lst[,1])
combined <- AddModuleScore(object = combined, features = lst, name = 'DARopen_', random.seed = 1)

lst <- read.table("data/DARclose_combinedRegion.txt")
lst <- list(lst[,1])
combined <- AddModuleScore(object = combined, features = lst, name = 'DARclose_', random.seed = 1)

FeaturePlot(
  object = combined,
  features = "DARclose_1",
)

VlnPlot(
  object = combined,
  features = "DARclose_1",
  split.by = "dataset"
)

pdf(paste0("seed",sd,"_scATAC.pdf"), width = 6, height = 5)
DimPlot(combined, pt.size = 0.1)
DimPlot(combined, group.by = 'dataset', pt.size = 0.1)
FeaturePlot(
  object = combined,
  features = "DARopen_1",
)
VlnPlot(
  object = combined,
  features = "DARopen_1",
  split.by = "dataset"
)
dev.off()

saveRDS(combined,paste0("seed",sd,"_scATAC.rds"))
