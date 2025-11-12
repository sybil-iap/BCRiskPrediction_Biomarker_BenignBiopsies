# Breast Cancer Biomarker Development for Benign Biopsies

## RNAseq Data files on HGCC
* `LBI12534-70650_Batch_1_2021/` : First batch of ~90 samples
* `LBI13454-118133_Batch_2_2024/` : Second batch of ~50 samples
	* `LBI**/RawFastq/` : Raw fastq files
	* `LBI**/**QC/` : QC files of mRNA biological samples and raw read files

### Mapped output files by STAR under `LBI**/`

1. `./Sample_ID.txt` : List of sample IDs.
2. `./filenames_**.txt` : Match sequence ID to sample IDs.
2. `/MAP_OUT/*Log.out` : Tool variables used for mapping
3. `/MAP_OUT/*Log.final.out` : Mapping quality metrics
4. `/MAP_OUT/*.SJ.out.tab` : Splicing junction output files
5. **`/MAP_OUT/*.ReadsPerGene.out.tab` : Read counts per gene**
	 * column 1: gene ID column 
	 * **column 2: counts for unstranded RNA-seq column (Used for DGE Analysis)**
	 * column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes) 
	 * column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)
	 * Note, if you have stranded data and choose one of the columns 3 or 4, the other column (4 or 3) will give you the count of antisense reads.
	 * Count the number of reads mapped to each strand by using a simple awk script
`grep -v "N_" p21035-s024_bng_bx_25_S24ReadsPerGene.out.tab | awk '{unst+=$2;forw+=$3;rev+=$4}END{print unst,forw,rev}'`
6. `/LBI**_BAM/*.Aligned.toTranscriptome.out.bam` : BAM files for reads aligned to transcriptome
7. `/LBI**_BAM/*.sortedByCoord.out.bam` : Sorted BAM files for aligned reads
8. `/LBI**_CRAM/*.sortedByCoord.out.cram` : Sorted CRAM files for aligned reads (converted from BAM files by samtools)

### Additional meta data under `LBI**/Misc/`

## Raw data and Processed Data under `BreastCancer_PredictiveBiomarker/data`
- `rawRNA`
- `rawPhenotype`
- `ProcessedRNa`
- `ProcessedPhenotype`

## Scripts under `BreastCancer_PredictiveBiomarker/scripts`
- `Process_Phenotype.rmd`
- `Process_RNAseq.rmd`
- `Process_DEGs.rmd`
- `pam50.r`
- `Rfunc.r`: self-defined R functions