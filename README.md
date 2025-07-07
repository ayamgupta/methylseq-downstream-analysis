# Methylation Analysis Pipeline

This Nextflow pipeline performs differential methylation analysis from Bismark or methyldackel output. It identifies differentially methylated regions, annotates them to genes, and performs enrichment analysis.

# Input Files

samplesheet - present in the root dir for testing. This will work with both the file types.

genes.gtf s3 path - s3://refdata-1/nfcore-pipelines/bio-data/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf

input-bismark-path - s3://quark-dev-data/quarkdata/home/quark/ayam-invisibl-io/myresults/jobs/methylseq-bismark-b1-tt294/

input-methyldackel-path - s3://quark-dev-data/quarkdata/home/quark/ayam-invisibl-io/myresults/jobs/methylseq-bwa-b1-9vrg9/

## Quickstart

1.  **Build Docker Image:**
    ```bash
    docker build -t methylseq .
    ```

2.  **Run the Pipeline:**
    ```bash
    nextflow run main.nf --samplesheet <path_to_samplesheet> --input_dir <path_to_input_data> --gtf <path_to_gtf>
    ```

## Parameters

### Mandatory

*   `--samplesheet`: Path to the samplesheet CSV file. The file should contain the following columns:
    *   `sample`: Sample name
    *   `cohort`: Cohort or group for each sample
    *   `genome`: Genome assembly (e.g., `hg38`)
*   `--gtf`: Path to the gene transfer format (GTF) file for gene annotation.

### Optional
*   `--input_dir`: Path to the directory containing methylation data (Bismark coverage or methyldackel output).
*   `--file_type`: Type of input files. Choose from `bismark` or `methyldackel`. (Default: `bismark`)
*   `--outdir`: Output directory for results. (Default: `./results`)
*   `--min_coverage`: Minimum read coverage to include a CpG site. (Default: `10`)
*   `--meth_diff_threshold`: Absolute methylation difference threshold for calling significant regions. (Default: `25`)
*   `--qval_threshold`: Q-value (adjusted p-value) threshold for significance. (Default: `0.05`)
*   `--top_gene_count`: Number of top hyper- and hypo-methylated genes to select for enrichment analysis. (Default: `10`)
*   `--enrichr_databases`: Comma-separated list of Enrichr databases for enrichment analysis. (Default: See `nextflow.config`)

## Examples

### Bismark Analysis

This example runs the pipeline on Bismark output, with a minimum coverage of 10 and a q-value threshold of 0.01.

```bash
nextflow run main.nf \
    --samplesheet samplesheet.csv \
    --input_dir bismark/ \
    --gtf genes.gtf \
    --file_type bismark \
    --min_coverage 10 \
    --qval_threshold 0.01 \
    --outdir results_bismark/
```

### Methyldackel Analysis

This example runs the pipeline on methyldackel output, with a higher minimum coverage and a custom output directory.

```bash
nextflow run main.nf \
    --samplesheet samplesheet.csv \
    --input_dir methyldackel/ \
    --gtf genes.gtf \
    --file_type methyldackel \
    --min_coverage 25 \
    --outdir results_methyldackel/
```

Tested With the following command
nextflow run main.nf \
    --samplesheet methylseq-samplesheet.csv \
    --gtf genes.gtf \
    --file_type bismark \
    --min_coverage 10 \
    --qval_threshold 0.01 \
    --outdir results_bismark/