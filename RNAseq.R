#------------------------------------------------------------------------
# 1. Load required packages
#------------------------------------------------------------------------
library(DESeq2)
library(edgeR)
library(pheatmap)
library(miodin)
library(ggplot2)

#------------------------------------------------------------------------
# 2. Define paths and parameters
#------------------------------------------------------------------------
data_dir <- "C:/Users/Dell/Documents/Assignment_1_a24ahmou/"  # Update to your actual path
file_name <- "a24ahmou_rnaseq.txt"     # Replace with your dataset name
count_path <- file.path(data_dir, file_name)

#------------------------------------------------------------------------
# 3. Import data
#------------------------------------------------------------------------
if (!file.exists(count_path)) stop("Count file not found at specified path")

# Read count table
countTable <- read.table(count_path,
                         header = TRUE, 
                         as.is = TRUE, 
                         row.names = 1, 
                         sep = "\t")

#------------------------------------------------------------------------
# 4. Dynamically extract sample groups from names
#------------------------------------------------------------------------
sample_names <- colnames(countTable)

# Extract group from sample names (e.g., "Control_1" -> "Control")
sample_groups <- sapply(strsplit(sample_names, "_"), `[`, 1)

# Create colData with inferred groups
col_data <- data.frame(
  row.names = sample_names,
  group = factor(sample_groups)
)

#------------------------------------------------------------------------
# 5. Filter low counts
#------------------------------------------------------------------------
keep <- rowMeans(countTable) >= 1
countTable <- countTable[keep, , drop = FALSE]

#------------------------------------------------------------------------
# 6. Prepare DESeq2 dataset (unpaired design)
#------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = countTable,
  colData = col_data,
  design = ~ group  # Unpaired design
)

#------------------------------------------------------------------------
# 7. Miodin Project Setup
#------------------------------------------------------------------------
project_dir <- file.path(".", "a24ahmou_RNAseqAnalysis")

# Remove existing project folder if it exists
if (dir.exists(project_dir)) {
  unlink(project_dir, recursive = TRUE)
}

# Create new project
mp <- MiodinProject(
  name = "a24ahmou_RNAseqAnalysis",
  author = "YourName",
  path = ".",
  overwrite = TRUE
)

#------------------------------------------------------------------------
# 8. Study Design (Case-Control)
#------------------------------------------------------------------------
# Get unique groups from sample names
unique_groups <- unique(sample_groups)
caseName <- unique_groups[1]
controlName <- unique_groups[2]

# Create sample and assay tables
sampleTable <- data.frame(
  SampleName = sample_names,
  SamplingPoint = paste0("sp", rep(1, length(sample_names))),
  Group = sample_groups
)

assayTable <- data.frame(
  SampleName = sample_names,
  DataFile = count_path,
  DataColumn = sample_names
)

# Define study design with unpaired configuration
ms <- studyDesignCaseControl(
  studyName = "CustomStudy",
  factorName = "Group",  # Matches column name in sampleTable
  caseName = caseName,
  controlName = controlName,
  contrastName = paste0(caseName, "_vs_", controlName),
  numCase = sum(sample_groups == caseName),
  numControl = sum(sample_groups == controlName),
  paired = "No",  # Fixed: Use "No" for unpaired samples
  sampleTable = sampleTable,
  assayTable = assayTable,
  assayTableName = "RNAseq"
)

# Insert into project
insert(ms, mp)

#------------------------------------------------------------------------
# 9. Workflow Definition
#------------------------------------------------------------------------
mw <- MiodinWorkflow("RNAseq workflow")

# Step 1: Import data
mw <- mw + 
  importProcessedData(
    name = "RNA importer",
    experiment = "sequencing",
    dataType = "rna",
    studyName = "CustomStudy",
    assayName = "RNAseq",
    datasetName = "CustomRNAseq",
    contrastName = paste0(caseName, "_vs_", controlName)
  )

# Step 2: Process data
mw <- mw + 
  processSequencingData(
    name = "RNA processor",
    contrastName = paste0(caseName, "_vs_", controlName),
    filterLowCount = TRUE
  )

# Insert workflow into project
mw <- insert(mw, mp)

#------------------------------------------------------------------------
# 10. Execute Analysis
#------------------------------------------------------------------------
mw <- execute(mw)

# Re-run variance-stabilizing transformation (VST) instead of rlog
normCounts <- varianceStabilizingTransformation(dds, blind = FALSE)

# Optional: Visualize sparsity issue (for diagnostics)
# plotSparsity(dds)

# Create directory for QC reports
qc_dir <- file.path(project_dir, "Exports", "Datasets", "CustomRNAseq", "qualityReports", "RNAseq", "Count data")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

# PCA Plot
pca_plot_path <- file.path(qc_dir, "PCA_plot.png")
png(pca_plot_path, width = 800, height = 600)
plotPCA(normCounts, intgroup = "group")
dev.off()

# Sample-to-sample distance heatmap
sampleDist <- cor(assay(normCounts), method = "spearman")
heatmap_plot_path <- file.path(qc_dir, "Heatmap.png")
png(heatmap_plot_path, width = 800, height = 600)
pheatmap::pheatmap(
  as.matrix(sampleDist),
  clustering_distance_rows = as.dist(1 - sampleDist),
  clustering_distance_cols = as.dist(1 - sampleDist),
  annotation_col = col_data
)
dev.off()

#------------------------------------------------------------------------
# 11. Save and Export
#------------------------------------------------------------------------
saveDataFile(mp)
export(mp, "dataset", "CustomRNAseq")