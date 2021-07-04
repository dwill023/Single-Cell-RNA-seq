# Single-Cell-RNA-seq

One of the most popular single cell isolation techniques is with the Drop-seq approach, where a single cell is encapsulated into a micordroplet containing a bead with unique barcodes, primers and enzymes where cDNA synthesis and library generation is performed.

![](https://github.com/dwill023/Single-Cell-RNA-seq/blob/main/figures/drop-seq.png) <sub>[Macosko et al., 2015](https://www.cell.com/cell/fulltext/S0092-8674(15)00549-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867415005498%3Fshowall%3Dtrue)</sub>

One of the most widely used technologies for the Drop-seq approach is the 10x Genomics Chromium platform. There has been a study by [Svensson et al,2017](https://www.nature.com/articles/nmeth.4220) comparing the different scRNA methods. If transcript level quantificaion and detection is what your experiment is aiming for, then the drop-seq (10x chromium) is fine. If transcript isoforms are part of the study then the full-length method (Smart-seq2) is required. 

Since my lab focuses on comparing transcript levels I will only outline the experimental steps for the 10x chromium protocol. Below is a schematic of the chromium chip where the cells that have been enzymatically broken up and mixed with reagent into the 2nd row of the chip are added. The microfluidics will mix these cells with the beads and oil to form the gel emulsion (GEMs) at the top row. The GEMs from this row are taken off the chip and put on a thermocycler for cDNA synthesis which will also incorporate a unique barcode. Then the GEMs are broken open, amplified, and cleaned up for sequencing. The remainder of the schematic outlines the data processing and visualization steps. 

![](https://github.com/dwill023/Single-Cell-RNA-seq/blob/main/figures/workflow.png)

When setting up the experiment it is **important** that samples should be balanced, evenly distributed, across all stages of the experiment. This will reduce sources of technical variation in the experiment. 

For example, you have samples on Day 0 sequenced on one flowcell and Day 7 samples run on another flowcell. The variation you observe can not be determined to be from the stage, or from technical variation in the sequencing run. Therefore, samples from day 0 and day 7 should be run on both flowcells. 

There are many similarities to the processing of scRNA-seq and the traditional RNA-seq. Both of which go through the same initial processing of read quality assessment, alignment and mapping quality assessment. 

The major difference is each library in a scRNA-seq represents one cell instead of a population. This creates some considerations on how the data can be filtered before analyzing. Things to consider are:

- library amplication depth: each cell can have differing number of reads between them.
- Gene 'dropouts': A gene may have a moderate amount reads in one cell but not the other.

The above can be introduced due to low starting material. A way to mediate this is to have lots of cells. It's difficult to estimate how may different types of cells you have in a sample but there's a [calculator](https://satijalab.org/howmanycells/) that can estimate this based on some assumptions. Based on the [10x literature](https://support.10xgenomics.com/single-cell-gene-expression/library-prep/doc/user-guide-chromium-single-cell-3-reagent-kits-user-guide-v3-chemistry) and reading papers (PMID: 32795399, PMID: 32302522) I would say 5,000 cells per condition per replicate is a good enough for complex samples (20 different cell types).

The **sequencing depth** recommended is a minimum of 20,000 read pairs per cell.  Paired-ended 50 bp reads.

However, it has been reported that the optimal allocation is to sequence one read per cell per gene = 
- For humans: ~21,000 read pairs per cell
- For mouse:  ~25,000 read pairs per cell

An example of a sequencing run using the 10x chromium chip is outlined below. In this example you have multiple samples that are processed through multiple GEM wells which generate multiple libraries that are pooled into one flowcell. After demultiplexing, the 10x software cellranger performs a count separately for each GEM well. For example, if you have 6 samples (6 GEM wells) you have to run `cellranger count` six times. Then you can aggregate them with a single instance of `cellranger aggr`. 

![](https://github.com/dwill023/Single-Cell-RNA-seq/blob/main/figures/cellranger.png)

The [cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) software performs the following:
- cellranger mkfastq Demultiplexes raw base call (BCL) files and  outputs fastq. If sequencing facility hasn't already done so.
- cellranger count Performs alignment, filtering, barcode counting, and UMI counting.
- aggregates outputs from multiple runs of cellranger count, normalizing those runs to the same sequencing depth. See [Multi-Library Aggregation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate).
	
The output of the above is a UMI count matrix. The values in this matrix represent the number of molecules for each feature (i.e. gene; row) that are detected in each cell (column). This file is used with the [Seurat R library](https://satijalab.org/seurat/archive/v3.2/pbmc3k_tutorial.html) to select and filter cells based on QC metrics, data normalization & scaling, and the detection of highly variable features. 

```diff
@@ Below are the steps to analyze the data starting with fastq files. @@
```

# Read quality assessment
Run fastqc and MultiQC for read quality and quantitiy.
```bash
fastqc *.fastq.gz
```

```bash
multiqc .
```

# Trimming
If trimming of the low quality bases and adapters are required there are different trimming tools available.

[Trim Galore](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md), that uses Cutadapt and FastQC, which will do this and automatically generate fastQ files afterwards. It also has built in adapter auto-detection. If you want to use trim galore, Cutadapt and FastQC have to be preinstalled.

```bash
# Check that cutadapt is installed
cutadapt --version
# Check that FastQC is installed
fastqc -v
# Install Trim Galore
curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz -o trim_galore.tar.gz
tar xvzf trim_galore.tar.gz
# Run Trim Galore
~/TrimGalore-0.6.6/trim_galore -o fastqc_trimmed_results Sample_1.fastq Sample_2.fastq
```
# Demultiplexing
This step is performed to demultiplex the UMI (12 nt unique molecular identifier) which is used to identify the single cell.
For Smartseq2 or other paired-end full transcript protocols the data will usually already be demultiplexed.
If the sequencing facility hasn't already demultiplexed you have to do it yourself.
Also in the droplet-based protocol  the cell-barcode will be attached to the read name and demultiplexing will happen during the quantification step.

# Cell Counting
Most experiments will have multiple samples, processed from multiple wells of a chromium chip (8 max). Each well represents a treatment, timepoint, or condition and has its own barcode to distinguish it from the other wells. 
You must perform the count step separately for each well, see the figure above.

The cellranger software requires lots of memory (64 GB) and 16 cores (CPU). Use a HPCC to run it.

```bash
cellranger count --id=day1 \
--fastqs=/scRNA_seq/d1_N2B27/ \
--transcriptome=refdata-gex-mm10-2020-A
```
The output files are listed below:
![](https://github.com/dwill023/Single-Cell-RNA-seq/blob/main/figures/output1.png)

Summary of what each file is:

|             File Name            |                                                                             Description                                                                            |
|:--------------------------------:|:------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
| web_summary.html                 | Run summary metrics and charts in HTML format                                                                                                                      |
| metrics_summary.csv              | Run summary metrics in CSV format                                                                                                                                  |
| possorted_genome_bam.bam         | Reads aligned to the genome and transcriptome annotated with barcode information                                                                                   |
| possorted_genome_bam.bam.bai     | Index for possorted_genome_bam.bam                                                                                                                                 |
| filtered_feature_bc_matrix       | Filtered feature-barcode matrices containing only cellular barcodes in MEX format. (In Targeted Gene Expression samples, the non-targeted genes are not present.)  |
| filtered_feature_bc_matrix_h5.h5 | Filtered feature-barcode matrices containing only cellular barcodes in HDF5 format. (In Targeted Gene Expression samples, the non-targeted genes are not present.) |
| raw_feature_bc_matrices          | Unfiltered feature-barcode matrices containing all barcodes in MEX format                                                                                          |
| raw_feature_bc_matrix_h5.h5      | Unfiltered feature-barcode matrices containing all barcodes in HDF5 format                                                                                         |
| analysis                         | Secondary analysis data including dimensionality reduction, cell clustering, and differential expression                                                           |
| molecule_info.h5                 | Molecule-level information used by cellranger aggr to aggregate samples into larger datasets                                                                       |
| cloupe.cloupe                    | Loupe Browser visualization and analysis file                                                                                                                      |
| feature_reference.csv            | (Feature Barcode only) Feature Reference CSV file                                                                                                                  |
| target_panel.csv                 | (Targeted GEX only) Targed panel CSV file                                                                                                                          |

# Aggregating Multiple GEM Wells

A CSV file must first be created with the following columns:
- library_id: Unique identifier for this input GEM well. This will be used for labeling purposes only; it doesn't need to match any previous ID you've assigned to the GEM well.
- molecule_h5: Path to the molecule_info.h5 file produced by `cellranger count`.

| library_id | molecule_h5                                    |
|------------|------------------------------------------------|
| day0_naive | /scRNA_seq/d0_2iL/d0_2iL/outs/molecule_info.h5 |
| day1_epi   | /scRNA_seq/d1_N2B27/day1/outs/molecule_info.h5 |
| day2_epi   | /scRNA_seq/d2_N2B27/day2/outs/molecule_info.h5 |

if you have different batches such as samples run with different 10X kits ( like version 2 or version 3) you can set up the above table with another column labeled "batch" and the software will perform batch correction.

Run cellranger aggr without normalization (none) this is because we will use Seurat to perform the normalization and it requires un-normalized counts.
```bash
cellranger aggr --id=mouseNaive_epiblast \
--csv=data.csv \
--normalize=none
```

# Seurat
Use the file filtered_feature_bc_matrix, generated from above step, as input into seurat. The filtered feature file is used as it contains only the detected cell barcodes without the background.

Follow the steps in the RMD file: [seurat_scRNAseq.Rmd](seurat_scRNAseq.Rmd).


