#!/bioinfo/users/tberthom/miniforge3/envs/newenvtest/bin/R

# -------------------------------------------------------------------
# Script Name: diff_analysis_peaks_LIBRA.R
# Author : Thomas Berthomier
# Description:
#   This script performs differential peak analysis using ArchR and the
#   Libra pseudobulk framework. It can optionally perform peak calling and
#   filter samples; it computes differential accessibility between groups,
#   and runs motif enrichment on significant peaks. Final results are
#   saved as an RDS file.
#
# ----------------------------- INPUTS ------------------------------
#
# Required arguments:
#   --archr_path             : Path to the ArchR project directory (string)
#   --out_directory          : Path to the output file directory (string)
#   --diff_annotation        : Column name in cellColData used for comparison (string)
#   --group_forward          : Label for the test (foreground) group (string)
#   --group_background       : Label for the background (control) group (string)
#
# Optional arguments:
#   --sample_included        : Vector of sample names to include (default = "All")
#   --cell_type_column       : Column with cell type annotations (default = "cell.type")
#   --replicate_column       : Column for replicate/patient IDs (default = "Patient.number")
#   --sample_column          : Column name specifying sample identifier (string, default = "Sample.number")
#   --control_column         : Column for batch correction (default = "No_control")
#   --paired                 : Whether to use paired analysis on replicates (default = FALSE)
#   --de_method              : DE method to use (e.g., "DESeq2") (default = "DESeq2")
#   --de_type                : Type of test used in DE (e.g., "Wald") (default = "Wald")
#   --binarization           : Whether to binarize matrix before DE (default = FALSE)
#   --min_cell_threshold     : Minimum number of cells per pseudobulk group (default = 150)
#   --max_cells              : Max number of cells allowed per group (default = 5000)
#   --reproducibility        : Strategy for reproducible peaks (default = "max(ceiling(n/3),2)")
#   --exclude_patterns       : Groups from cell.type or diff_annotations to be excluded from peak calling count of pseudobulk selection (see comment)
#   --seed                   : Random seed (default = 123)
#   --motifdb                : Motif database (e.g., "JASPAR2020") (default = "JASPAR2020").
#                              Can be set to "custom", see custom_enrichment_path
#   --threshold_FDR          : FDR threshold for motif enrichment (default = 0.05)
#   --threshold_Log2FC       : Log2 fold change threshold for peaks (default = 0.5)
#   --peak_call              : Whether to call peaks de novo (default = FALSE)
#   --use_peak               : Path to a precomputed peak set (RDS) to skip peak calling (optional)
#                              It should correspond to proj@peakSet in an ArchRProj object.
#   --custom_enrichment_path : Path to a GRangesList object (RDS) for each GR corresponding to an annotation of interest 
#                              (methylation sites, chip seq peaks...), with names corresponding to the name of the annotation.
#
# ----------------------------- OUTPUTS -----------------------------
#
# Output:
#   - Combined results saved as RDS file containing:
#       * Differential peaks (SummarizedExperiment)
#       * Motif enrichment (up and down)
#       * Libra pseudobulk matrices and differential results
#
#   - Optionally:
#       * Updated ArchR project with called peaks and motif annotations
#
# -------------------------------------------------------------------


print("###################################################################")
print('############### Starting Peak Libra diff analysis #################')
print("###################################################################")
# libraries
suppressMessages(suppressWarnings(library(ArchR)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(argparse)))
suppressMessages(suppressWarnings(library(BSgenome.Hsapiens.UCSC.hg38)))

# parser

parser <- ArgumentParser(description = "filter and annotate the project")

parser$add_argument("--archr_path", type = "character", help = "archr project directory path")
parser$add_argument("--out_directory", type = "character", help = "output directory path")

parser$add_argument("--sample_included", nargs = "+", default = c("All"), help = "list of samples to include (use 'All' to include all)")
parser$add_argument("--diff_annotation", type = "character", help = "annotation column to use for differential analysis")
parser$add_argument("--group_forward", type = "character", help = "forward comparison group name")
parser$add_argument("--group_background", type = "character", help = "background comparison group name")
parser$add_argument("--cell_type_column", type = "character", default = "cell.type", help = "column name for cell type")
parser$add_argument("--replicate_column", type = "character", default = "Patient.number", help = "column name for identifying replicates")
parser$add_argument("--sample_column", type = "character", default = "Sample.number", help = "column name for sample")
parser$add_argument("--control_column", type = "character", default = "No_control", help = "column name for controlling batch")
parser$add_argument("--paired", type = "logical", default = FALSE, help = "should i pair on replicate")
parser$add_argument("--de_method", type = "character", default = "DESeq2", help = "differential expression method")
parser$add_argument("--de_type", type = "character", default = "Wald", help = "differential expression type")
parser$add_argument("--binarization", type = "logical", default = FALSE, help = "whether to binarize input matrix")
parser$add_argument("--min_cell_threshold", type = "integer", default = 150, help = "minimum number of cells per group/sample")
parser$add_argument("--max_cells", type = "integer", default = 5000, help = "maximum number of cells per group")
parser$add_argument("--reproducibility", type = "character", default = "max(ceiling(n/3),2)", help = "peak reproducibility strategy")
parser$add_argument("--exclude_patterns", type = "character", nargs = "+", default = NULL, help = "see inside the peak calling comment")
parser$add_argument("--seed", type = "integer", default = 123, help = "random seed")

parser$add_argument("--motifdb", type = "character", default = "JASPAR2024", help = "Motif DB used in motif enrichment")
parser$add_argument("--threshold_FDR", type = "numeric", default = 0.05, help = "Threshold of adjusted pvalue for peaks")
parser$add_argument("--threshold_Log2FC", type = "numeric", default = 0.5, help = "Threshold of Log2FC for peaks")
parser$add_argument("--peak_call", type = "logical", default = FALSE, help = "whether I should do the peak calling")
parser$add_argument("--use_peak", type = "character", default = NULL, help = "path to peakset i should use")
parser$add_argument("--custom_enrichment_path", type = "character", default = NULL, help = "path to GRlist rds that I should use")

# Parse arguments
args <- parser$parse_args()

proj_global_dir <- args$archr
OUT_LIBRA <- args$out_directory

CODE_PATH <- "/mnt/beegfs/home/alyne/persistent/ArchR_functions_Thomas" #TO CHANGE
 
sample_included <- args$sample_included
diff_annotation <- args$diff_annotation
group_forward <- args$group_forward
group_background <- args$group_background
cell_type_column <- args$cell_type_column
replicate_column <- args$replicate_column
sample_column <- args$sample_column
control_column <-  args$control_column
paired <- args$paired
de_method <- args$de_method
de_type <- args$de_type
binarization <- args$binarization
minCellThreshold <- args$min_cell_threshold
maxCells <- args$max_cells
reproducibility_strategy <- args$reproducibility
SEED <- args$seed
motifdb <- args$motifdb
threshold_FDR <- args$threshold_FDR
threshold_LogFC <- args$threshold_Log2FC
peak_call <- args$peak_call
exclude_patterns <- args$exclude_patterns
use_peak <- args$use_peak

set.seed(SEED)
datetime <- strftime(Sys.time(), format = "%Y:%m:%d_%H:%M")

# functions scripts
source_all_r_files <- function(path_list) {
  origin_wd <- getwd()
  lapply(
    path_list, function(x) {
      setwd(x)
      file.sources <- list.files(pattern = "*.r$")
      sapply(file.sources, source, .GlobalEnv)
    }
  )
  setwd(origin_wd)
}
source_all_r_files(c(file.path(CODE_PATH, "Miscellaneous"), file.path(CODE_PATH, "Differential_analysis")))

# ---- PRE-CHECKS ----
if (!is.numeric(minCellThreshold) || minCellThreshold <= 0) {
  stop("minCellThreshold must be a positive number.")
}

# Vérifier fichier de pics si fourni
if (!is.null(use_peak) && !file.exists(use_peak)) {
  stop(paste("Provided peak set file not found:", use_peak))
}

# Vérifier enrichissement custom si utilisé
if (tolower(motifdb) == "custom" && !file.exists(custom_enrichment_path)) {
  stop(paste("Custom enrichment file not found:", custom_enrichment_path))
}

# Vérifier que le dossier de sortie existe (sinon le créer)
check_a_path(OUT_LIBRA)

# Vérifier que la base de motifs est connue
valid_motifdb <- c("cisbp", "encode", "homer", "jaspar2016", "jaspar2020", "jaspar2018", "jaspar2024", "vierstra", "custom")
if (!(tolower(motifdb) %in% valid_motifdb)) {
  stop(paste("Unknown motif database:", motifdb))
}
# ---- END PRE-CHECKS ----

# set up
addArchRThreads(ceiling(parallelly::availableCores() / 2))
addArchRGenome("hg38")
addArchRLocking(locking = T)

# load project
proj_global <- loadArchRProject(proj_global_dir)

# call peaks using method described in my report if peak_call == TRUE
if (peak_call) {

  # create a combined annotation of anno of interest (ex : treatment) + cell type
  proj_global <- addCombinedAnnotation(proj_global, diff_annotation, cell_type_column)
  combined_name <- paste(diff_annotation, "by", cell_type_column, sep = "_")

  # count number of cells per combined annotation per sample
  count_cells <- proj_global@cellColData[, c(sample_column, combined_name)] %>%
    as.data.frame() %>%
    group_by(.data[[combined_name]], .data[[sample_column]]) %>%
    summarize(count = n(), .groups = "drop")
  
  # count number of pseudobulk samples passing the minCell Threshold
  count_per_condition <- count_cells %>% filter(count >= minCellThreshold) %>%
    group_by(.data[[combined_name]]) %>%
    summarize(num_samples = n(), .groups = "drop")

  print("Number passing conditions :")
  print(count_per_condition)


  # This code sets the global min and max number of pseudobulk replicates that
  # ArchR can create per annotation group (e.g., State_CellType).
  # ArchR only accepts a single min/max value across all groups.
  # We compute these values based on how many biological replicates (samples)
  # are available per group.

  # To ensure quality, we exclude unreliable or uninformative groups using
  # `exclude_patterns` (e.g., "undetermined", "rec2").
  # You can customize `exclude_patterns` to skip annotation values that should
  # not influence the min/max estimation.

  # Importantly, include in `exclude_patterns` all groups with very few cells,
  # as they can distort the replicate count and degrade pseudobulk quality.

  # The smallest number of replicates across valid groups becomes `minReplicate`.
  # We then set `maxReplicate` at the max + 2 to allow ArchR flexibility.
  # This ensures each group has enough samples while still allowing large groups
  # to be split into multiple pseudobulks if enough cells are present.

  # default values
  if (is.null(exclude_patterns)) {
    if (diff_annotation == "State") {
      exclude_patterns <- c("undetermined", "rec2")
    } else {
      exclude_patterns <- c("undetermined")
    }
  }

  minReplicate <- count_per_condition %>%
    filter(
      !stringr::str_detect(.data[[combined_name]], str_c(exclude_patterns, collapse = "|"))
    ) %>%
    summarize(min(num_samples)) %>%
    as.numeric()

  maxReplicate <- count_per_condition %>%
    filter(
      !stringr::str_detect(.data[[combined_name]], str_c(exclude_patterns, collapse = "|"))
    ) %>%
    summarize(max(num_samples)) %>% as.numeric() + 2

  proj_global <- addGroupCoverages(
    ArchRProj = proj_global,
    groupBy = combined_name,
    useLabels = TRUE,
    sampleLabels = sample_column,
    minCells = minCellThreshold,
    maxCells = maxCells,
    minReplicates = minReplicate,
    maxReplicates = maxReplicate,
    force = TRUE
  )

  proj_global <- addReproduciblePeakSet(
    ArchRProj = proj_global,
    groupBy = combined_name,
    reproducibility = reproducibility_strategy,
    force = TRUE
  )

  proj_global <- addPeakMatrix(proj_global)

  proj_global <- saveArchRProject(proj_global)

} 
#delete this, peak set and matrix already in loaded archr object
# else if (!is.null(use_peak)) {
#   # If peak set is false and you provide a peakSet rds object
#   # will load it
#   proj_global@peakSet <- readRDS(use_peak)
#   proj_global <- addPeakMatrix(proj_global)
#   proj_global <- saveArchRProject(proj_global)
# }


# filters for given samples
if (!identical(sample_included, c("All"))) {
  idx_condition <- BiocGenerics::which(proj_global@cellColData[,sample_column] %in% sample_included)
  cellsSample <- proj_global$cellNames[idx_condition]
  proj_global <- proj_global[cellsSample, ]
}

# this is for output naming issues only
if (control_column == "No_control") {
  control_col <- NULL
} else {
  control_col <- control_column
}

# for numerically coded categories
if (is.numeric(proj_global@cellColData[, diff_annotation])) {
  warning(paste0(diff_annotation, " is a numerical variable. Converting to factor..."))
  proj_global@cellColData[, diff_annotation] <- as.factor(proj_global@cellColData[, diff_annotation])
}

# out_libra_pseudobulk <- getLibraPseudoBulkAnalysis(proj_global, diff_annotation, group_forward, group_background,
#                    cell_type = cell_type_column,
#                    replicate_col = replicate_column,
#                    control_col = control_col,
#                    paired_by_replicate = paired,
#                    useMatrix = "PeakMatrix",
#                    de_method = de_method,
#                    de_type = de_type,
#                    minCell= minCellThreshold, 
#                    n_threads = ceiling(parallelly::availableCores() / 2)
#                    )
out_libra_pseudobulk <- getLibraPseudoBulkAnalysis(proj_global, diff_annotation, cell_type_column, group_forward, group_background,
                                                   replicate_col = replicate_column,
                                                   control_col = control_col,
                                                   paired_by_replicate = paired,
                                                   useMatrix = "PeakMatrix",
                                                   de_method = de_method,
                                                   de_type = de_type,
                                                   minCell= minCellThreshold, 
                                                   n_threads = ceiling(parallelly::availableCores() / 2)
)

peak_global <- getPeakSet(proj_global)
results_se <- convertLibraOutToSE(out_libra_pseudobulk$results, peak_global)

if (tolower(motifdb) == "custom") {
  custom_enrich_grlist <- readRDS(custom_enrichment_path)
  proj_global <- addPeakAnnotations(ArchRProj = proj_global, regions = custom_enrich_grlist, name = "Motif")
} else {
  # motif analysis
  # need to load jaspar because of a bug in ArchR
  if (grepl("jaspar", tolower(motifdb), fixed = TRUE)) {
    library(motifdb, character.only = TRUE)
  }

  # need to add this line as JASPAR2024 is not natively used in archr
  # we provide it as a custom motif db
  if (tolower(motifdb) == "jaspar2024") {
    # jaspar <- JASPAR2024()
    # sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(jaspar))
    # motifs24 <- TFBSTools::getMatrixSet(sq24, list(species = "Homo sapiens", collection = "CORE"))
    # motif_mat <- TFBSTools::toPWM(motifs24)
    motif_mat = readRDS("/mnt/beegfs/home/alyne/persistent/motif_mat_JASPAR2024_human.rds")
    
    #change names so that we have the TF name instead of the ID
    adjust_labels <- function(vec) {
      dup_vals <- duplicated(vec) | duplicated(vec, fromLast = TRUE)
      vec[dup_vals] <- paste0(vec[dup_vals], "_", names(vec)[dup_vals])
      return(vec)
    }
    name_matrices <- sapply(motif_mat, function(pfm) TFBSTools::name(pfm))
    deduplicated_name_matrices <- adjust_labels(name_matrices)

    names(motif_mat) <- deduplicated_name_matrices[names(motif_mat)]
    
    proj_global <- addMotifAnnotations(proj_global, motifPWMs = motif_mat, annoName = "Motif", force = TRUE)
  } else {
    proj_global <- addMotifAnnotations(proj_global, motifSet = motifdb, annoName = "Motif", force = TRUE)
  }
}

threshold_up = paste('FDR <=', threshold_FDR, '& Log2FC >=', threshold_LogFC, sep = " ")

motifsUp <- peakAnnoEnrichment(
  seMarker = results_se,
  ArchRProj = proj_global,
  peakAnnotation = "Motif",
  cutOff = threshold_up
)

threshold_down = paste('FDR <= ', threshold_FDR, ' & Log2FC <= -', threshold_LogFC, sep = "")

motifsDown <- peakAnnoEnrichment(
  seMarker = results_se,
  ArchRProj = proj_global,
  peakAnnotation = "Motif",
  cutOff = threshold_down
)

#combining results into one rds
results_fin <- list(Peaks =  results_se, MotifsUp = motifsUp, MotifsDown = motifsDown, matrices = out_libra_pseudobulk$matrices, Libra = out_libra_pseudobulk$results)

saveRDS(results_fin, file = file.path(OUT_LIBRA, paste0(paste(basename(proj_global_dir),"libra_diff", 
                                                              diff_annotation, group_background, group_forward,
                                                              control_column, de_method, de_type, sep = "_"),".rds")))

print("###################################################################")
print('############# Successful analysis for peaks analysis ##############')
print("###################################################################")
