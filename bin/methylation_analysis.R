#!/usr/bin/env Rscript

# Load necessary libraries
suppressPackageStartupMessages({
  library(argparse)
  library(methylKit)
  library(dplyr)
  library(genomation)
  library(GenomicRanges)
  library(txdbmaker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(ChIPpeakAnno)
  library(org.Hs.eg.db)
  library(enrichR)
  library(purrr)
  library(tidyr)
  library(readr)
})

# Function to parse command-line arguments
get_args <- function() {
  parser <- ArgumentParser(description = "Run methylation analysis pipeline.")
  
  parser$add_argument("--metadata_csv", type = "character", default = "samplesheet.csv", help = "Path to sample sheet CSV file")
  parser$add_argument("--base_dir", type = "character", default = NULL, help = "Optional base directory to prepend to relative file paths")
  parser$add_argument("--gtf_file", type = "character", default = "genes.gtf", help = "Path to GTF annotation file")
  parser$add_argument("--file_type", type = "character", default = "bismark", choices = c("bismark", "methyldackel"), help = "Type of input files")
  parser$add_argument("--min_coverage", type = "integer", default = 10, help = "Minimum coverage threshold")
  parser$add_argument("--high_percentile", type = "double", default = 99.9, help = "High percentile for coverage filtering")
  parser$add_argument("--meth_diff_threshold", type = "double", default = 25, help = "Methylation difference threshold")
  parser$add_argument("--qval_threshold", type = "double", default = 0.05, help = "Q-value threshold")
  parser$add_argument("--top_gene_count", type = "integer", default = 10, help = "Number of top genes to select")
  parser$add_argument("--min_genes_for_enrichment", type = "integer", default = 3, help = "Minimum genes required for enrichment")
  parser$add_argument("--output_dir", type = "character", default = "./", help = "Output directory")
  parser$add_argument("--enrichr_databases", type = "character", default = "GO_Molecular_Function_2023,GO_Biological_Process_2023,GO_Cellular_Component_2023,Reactome_2022,WikiPathway_2023_Human,DGIdb_Drug_Targets_2024", help = "Comma-separated list of EnrichR databases")
  
  parser$parse_args()
}

# Function to load and process metadata
load_metadata <- function(args) {
  metadata <- read.csv(args$metadata_csv, stringsAsFactors = FALSE)
  
  # Check for file path column (flexible naming)
  filepath_col <- NULL
  if ("path" %in% names(metadata)) {
    filepath_col <- "path"
  } else if ("filename" %in% names(metadata)) {
    filepath_col <- "filename"
  } else if ("filepath" %in% names(metadata)) {
    filepath_col <- "filepath"
  } else {
    stop("Samplesheet must contain either 'path', 'filename', or 'filepath' column with individual file paths")
  }
  
  # Check for cohort column (flexible naming)
  cohort_col <- NULL
  if ("cohort" %in% names(metadata)) {
    cohort_col <- "cohort"
  } else if ("group" %in% names(metadata)) {
    cohort_col <- "group"
  } else if ("phenotype" %in% names(metadata)) {
    cohort_col <- "phenotype"
  } else {
    stop("Samplesheet must contain either 'cohort', 'group', or 'phenotype' column for treatment grouping")
  }
  
  # Ensure genome column has values
  if ("genome" %in% names(metadata) && any(is.na(metadata$genome) | metadata$genome == "")) {
    warning("Empty genome values found. Setting to 'hg38' as default.")
    metadata$genome[is.na(metadata$genome) | metadata$genome == ""] <- "hg38"
  } else if (!("genome" %in% names(metadata))) {
    warning("No genome column found. Adding 'hg38' as default.")
    metadata$genome <- "hg38"
  }
  
  metadata <- metadata %>%
    mutate(
      # Use the file path from the samplesheet
      filename = get(filepath_col),
      # If base_dir is provided and paths are relative, prepend it
      filename = if (!is.null(args$base_dir) && !all(grepl("^/", filename))) {
        file.path(args$base_dir, filename)
      } else {
        filename
      },
      # Use the cohort column (flexible naming)
      cohort = get(cohort_col),
      condition_inferred = as.numeric(as.factor(cohort)) - 1
    ) %>%
    group_by(cohort) %>%
    mutate(sample.id = paste0(cohort, "_", row_number())) %>%
    ungroup()
  
  # Save cohort mapping
  unique_mapping <- metadata %>%
    dplyr::select(cohort, condition_inferred) %>%
    distinct() %>%
    arrange(cohort)
  write.csv(unique_mapping, file.path(args$output_dir, "cohort_to_condition_mapping.csv"), row.names = FALSE)
  
  # Compare with original condition column if it exists
  if ("condition" %in% names(metadata)) {
    comparison <- metadata %>%
      dplyr::select(sample, cohort, condition, condition_inferred) %>%
      mutate(match = condition == condition_inferred)
    write.csv(comparison, file.path(args$output_dir, "condition_comparison.csv"), row.names = FALSE)
  }
  
  # Check if all files exist
  missing_files <- metadata$filename[!file.exists(metadata$filename)]
  if (length(missing_files) > 0) {
    stop(paste("The following files do not exist:", paste(missing_files, collapse = ", ")))
  }
  
  metadata
}

# Function to process methylation data
process_methylation_data <- function(metadata, args) {
  pipeline_type <- if (args$file_type == "bismark") "bismarkCoverage" else "amp"
  
  methyl_objs <- mapply(
    FUN = methRead,
    location = metadata$filename,
    sample.id = metadata$sample.id,
    assembly = metadata$genome,
    treatment = metadata$condition_inferred,
    MoreArgs = list(pipeline = pipeline_type, mincov = args$min_coverage),
    SIMPLIFY = FALSE
  )
  
  myobj <- methylRawList(methyl_objs, treatment = metadata$condition_inferred)
  filtered <- filterByCoverage(myobj, lo.count = args$min_coverage, hi.perc = args$high_percentile)
  normed <- normalizeCoverage(filtered, method = "median")
  methylKit::unite(normed, destrand = FALSE)
}

# Function to perform PCA
perform_pca <- function(meth, output_dir) {
  pca <- PCASamples(meth, obj.return = TRUE)
  pca_df <- as.data.frame(pca$x) %>%
    mutate(
      sample = rownames(.),
      group = sub("_(\\d+)$", "", sample)
    )
  write_csv(pca_df, file.path(output_dir, "pca.csv"))
}

# Function to perform differential methylation analysis
perform_differential_methylation <- function(meth) {
  myDiff <- calculateDiffMeth(meth, overdispersion = "MN", adjust = "BH")
  getData(myDiff)
}

# Function to annotate genomic regions
annotate_regions <- function(diff_df, gtf_file) {
  gr <- GRanges(
    seqnames = diff_df$chr,
    ranges = IRanges(start = diff_df$start, end = diff_df$end),
    strand = diff_df$strand
  )
  
  txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
  genes_tx <- genes(txdb)
  
  peakAnno <- annotatePeakInBatch(gr, AnnotationData = genes_tx)
  gene_symbols <- mapIds(
    org.Hs.eg.db,
    keys = as.character(peakAnno$feature),
    column = "SYMBOL",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  
  diff_df$geneSymbol <- NA_character_
  hits <- findOverlaps(gr, peakAnno)
  diff_df$geneSymbol[queryHits(hits)] <- gene_symbols[as.character(peakAnno$feature[subjectHits(hits)])]
  diff_df$geneSymbol <- sapply(diff_df$geneSymbol, paste, collapse = ",")
  diff_df$geneSymbol[diff_df$geneSymbol == ""] <- NA_character_
  
  diff_df
}

# Function to select top differentially methylated regions
select_top_regions <- function(diff_df, args) {
  sig_df <- diff_df %>%
    filter(abs(meth.diff) >= args$meth_diff_threshold, qvalue < args$qval_threshold)
  
  top_hyper <- sig_df %>%
    arrange(desc(meth.diff)) %>%
    head(args$top_gene_count)
  top_hypo <- sig_df %>%
    arrange(meth.diff) %>%
    head(args$top_gene_count)
  
  list(
    hyper = top_hyper,
    hypo = top_hypo,
    combined = bind_rows(top_hyper, top_hypo)
  )
}

# Function to extract methylation data for heatmap
extract_methylation_for_heatmap <- function(meth, top_regions, output_dir) {
  sig_regions <- with(top_regions$combined, paste0(chr, "_", start, "_", end))
  meth_clean <- meth[complete.cases(meth), ]
  pm <- percMethylation(meth_clean, rowids = TRUE)
  pm_ids <- gsub("\\.", "_", rownames(pm))
  subset_pm <- pm[pm_ids %in% sig_regions, ]
  rownames(subset_pm) <- sub("^(([^.]+\\.[^.]+))\\..*", "\\1", rownames(subset_pm))
  subset_pm_df <- data.frame(Symbol = rownames(subset_pm), subset_pm, row.names = NULL)
  write.csv(subset_pm_df, file.path(output_dir, "grapher.csv"), row.names = FALSE)
}

# Function to run enrichment analysis
run_enrichment_analysis <- function(genes, label, db_list, out_dir, min_genes) {
  if (length(genes) < min_genes) {
    message(sprintf("Skipping enrichment for '%s': only %d gene(s) (minimum required: %d)", label, length(genes), min_genes))
    return(tibble())
  }
  
  enrich_results <- enrichr(genes, db_list)
  
  combined_list <- map(names(enrich_results), function(db_name) {
    result <- enrich_results[[db_name]]
    
    if (is.null(result) || nrow(result) == 0 || all(is.na(result$Term))) {
      return(tibble())
    }
    
    result %>%
      as_tibble() %>%
      filter(!is.na(Term), Term != "", !is.na(Genes), Genes != "") %>%
      separate_rows(Genes, sep = ";") %>%
      distinct() %>%
      dplyr::select(Genes, Term, P.value) %>%
      dplyr::rename(pvalue = P.value) %>%
      mutate(Database = db_name)
  })
  
  combined <- bind_rows(combined_list)
  
  if (nrow(combined) > 0) {
    write_csv(combined, file.path(out_dir, paste0("enrichment_", label, "_combined.csv")))
  }
  
  combined
}

# Main function to orchestrate the pipeline
main <- function() {
  args <- get_args()
  
  # Create output directory
  if (!dir.exists(args$output_dir)) {
    dir.create(args$output_dir, recursive = TRUE)
  }
  
  # Run the pipeline
  metadata <- load_metadata(args)
  meth <- process_methylation_data(metadata, args)
  perform_pca(meth, args$output_dir)
  diff_df <- perform_differential_methylation(meth)
  annotated_df <- annotate_regions(diff_df, args$gtf_file)
  write_csv(annotated_df, file.path(args$output_dir, "volcano.csv"))
  
  top_regions <- select_top_regions(annotated_df, args)
  extract_methylation_for_heatmap(meth, top_regions, args$output_dir)
  
  db_vector <- strsplit(args$enrichr_databases, ",")[[1]]
  
  run_enrichment_analysis(
    genes = na.omit(top_regions$hyper$geneSymbol),
    label = "hyper",
    db_list = db_vector,
    out_dir = args$output_dir,
    min_genes = args$min_genes_for_enrichment
  )
  
  run_enrichment_analysis(
    genes = na.omit(top_regions$hypo$geneSymbol),
    label = "hypo",
    db_list = db_vector,
    out_dir = args$output_dir,
    min_genes = args$min_genes_for_enrichment
  )
}

# Run the main function
if (!interactive()) {
  main()
}