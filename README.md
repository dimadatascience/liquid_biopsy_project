# Liquid Biopsy Project
edited by [Marcello Del Corvo](mailto:marcello.delcorvo@gmail.com)

<b>This repo provides details to understand and  reproduce all the analyses done for the upcoming paper (<i> Pelicci G. et al. </i>) from liquid biopsy project</b>

<blockquote>
The aim of the following work is to test the validity of this approach (cfDNA + evDNA sequencing for GBM cases) to identify somatic copy number
aberrations (sCNAs) in plasma samples in comparison with their matched tumor samples. We also want to benchmark sCNAs calling between WGS and WES tumor samples.
</blockquote>

## Contents
- [Contents](#contents)
- [Sample collection](#sample-collection)
- [Code](#code)
- [Analysis](#analysis)

## Sample collection

We performed low-pass WGS analysis of plasma (cfDNA + evDNA) from [***6 GBM patients***](./metafile/keys_sample_sequencing_ID.xlsx)  and their matched tumor samples and healty control (buffy coat).

With exception of ***C2-57*** patient, all samples were also sequenced with WES for tumor and control tissue.

<img src="img/samples_table.png" width="500" />

<i>Raw data used for the analyses are stored in the cluster.</i>

**1. Low-pass WGS (Plasma-Tumor-Control samples)**
```
cd /hpcnfs/techunits/genomics/PublicData/PelicciG/sfaletti/FASTQ/230919_A00302_0562_BHJ7C3DRX3
```
**2. WES (Tumor-Control samples)**

[Raw_data](./data/wes/raw_data_WES.txt) paths

## Code
All code used for analysis is provided at:

* [https://github.com/mdelcorvo/DeSeq-Free](https://github.com/mdelcorvo/DeSeq-Free) (branches: master, see README for usage details)

* [scripts](./scripts/) directory

## Analysis

At the surface level, the analyses can be broadly grouped into these sections:

<h3>Low-pass WGS</h3>
  
* Input: [metafile](./metafile/GBM_low-pass_WGS_samples.xlsx)
  
  a samplesheet  (can be .xlsx or .csv) with raw fastq.gz data that looks as follows:
  ```
  sample, lane, fq1, fq2, type
  
  E06, lane1, S52505_B-E06_S16_L001_R1_001.fastq.gz, S52505_B-E06_S16_L001_R2_001.fastq.gz, 0
  ```
  Each row represents a single-end fastq file. Rows with the same sample identifier are considered technical replicates and will be automatically merged. ``` type ``` refers to sample type (0=
  buffy coat, 1= plasma, 2=tumor).

  - Reference genome<br />
   
    Before starting, a user need to download reference genome. 

    Download from [NCBI](https://www.ncbi.nlm.nih.gov/genome/guide/human/), [Ensembl](https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/), or any other autorities
    ```
    wget https://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
    ```
    
    - Index reference genome for bwa-mem2<br />
  
      Prepare indexed genome for bwa-mem2 to boost mapping.  Refer to the [bwa-mem2 instruction](https://github.com/bwa-mem2/bwa-mem2).<br />
      
      Example code:
      ```
      ./bwa-mem2 index <in.fasta>
      Where 
      <in.fasta> is the path to reference sequence fasta file and 
      ```
      
* Code:

  ```
  git clone https://github.com/mdelcorvo/DeSeq-Free.git
  cd DeSeq-Free && conda env create -f envs/workflow.yaml
  conda activate DeSeq-Free_workflow

  snakemake --use-conda \
  --config \
  input=inputfile.xlsx \
  output=output_directory \
  genome=genome.fasta
  ```

  - <b>plasma - tumor CNAs comparions</b>
 
    file: [low_pass_wgs_analysis.R](./scripts/cna/low_pass_wgs_analysis.R)
 
    
    ```
    /hpcnfs/scratch/DIMA/delcorvo/liquid_biobsy_project/low_pass_wgs/scripts/low_pass_wgs_analysis.R
    
    # example of usage
    cd /hpcnfs/scratch/DIMA/delcorvo/liquid_biobsy_project/low_pass_wgs
    
    Rscript --vanilla scripts/low_pass_wgs_analysis.R \
    ./derived_data/cna/seg \ # input directory with IchorCNA output (.seg file)
    ./output \ # output directory
    ./derived_data/cna/canonical_exon_transcripts_hg38.bed # Gene annotation (.bed)
    ```
    
  - <b>fragment analysis</b>

    file: [fragment_size_distribution.py](./scripts/fragmentomics/fragment_size_distribution.py)
    
    ```
    /hpcnfs/scratch/DIMA/delcorvo/liquid_biobsy_project/low_pass_wgs/scripts/fragment_size_distribution.py
    
    # example of usage
    cd /hpcnfs/scratch/DIMA/delcorvo/liquid_biobsy_project/low_pass_wgs
    
    python scripts/fragment_size_distribution.py \
    --plasma_bam_list ./derived_data/recal/plasma_list.txt \
    --tumor_bam_list ./derived_data/recal/tumor_list.txt  \
    --output_csv derived_data/fragmentomics/gbm_fr.csv \
    --output_plot derived_data/fragmentomics/gbm_fr.png
    ```
  
- Results:

  All results are always stored also on cluster:
  ```
  /hpcnfs/scratch/DIMA/delcorvo/liquid_biobsy_project/low_pass_wgs
  ```

  - Quality control
  
      - coverage:

      <img src="results/qc/coverage/mosdepth.dist.png" width="700"/>
      
      - sequencing stats:
        [samtools_stats](./results/qc/samtools_stats/) directory
      
        [fastp_stats](./results/qc/fastp/) directory

      - fragmentomics metrics:
   
         [cfDNAPro](./results/qc/fragmentomics/) directory
      
         <img src="results/qc/fragmentomics/fragment_size_dist.png" width="700"/>

  - Alignment files
 
  ```
  /hpcnfs/scratch/DIMA/delcorvo/liquid_biobsy_project/low_pass_wgs/derived_data/recal
  
  ```
      
  - [CNAs](./results/cna/)
 
  ```
  /hpcnfs/scratch/DIMA/delcorvo/liquid_biobsy_project/low_pass_wgs/derived_data/cna/seg
  ```
    
  - [CNAs comparison](./results/cna/comparison/)
    
  ```
  /hpcnfs/scratch/DIMA/delcorvo/liquid_biobsy_project/low_pass_wgs/cna_comparison
  ```

  - [Fragmentomics](./results/fragmentomics/)
 
  ```
  /hpcnfs/scratch/DIMA/delcorvo/liquid_biobsy_project/low_pass_wgs/derived_data/fragmentomics
  ```
     
<h3>WES / WGS Benchmarking</h3>

Number of CNAs in WES confirmed by low-pass WGS
  
  CNAs were considered valid if:
  
    - overlap ≥ 50% between wgs and wes call
    - must be same type (e.g. duplication- duplication, deletion - deletion)


* Input data:

  [WES_CNAs call](./data/wes/): Dragen call on [WES Raw_data](./data/wes/raw_data_WES.txt)

  [WGS tumor - plasma comparison](./results/cna/comparison/CNAs_calls.txt)
  
  
* Code:

Filter for WES target regions

```
import pandas as pd
# File paths
tumor_plasma_file = "CNAs_calls.txt"
twist_bed_file = "Twist_hg38_CORRECT.bed"
output_file = "CNAs_calls.FILT.txt"

# Load the Twist BED file
twist_bed_df = pd.read_csv(twist_bed_file, sep="\t", header=None, names=["Chr", "Start", "End"])

# Ensure Twist BED columns are numeric for efficient comparisons
twist_bed_df["Start"] = pd.to_numeric(twist_bed_df["Start"])
twist_bed_df["End"] = pd.to_numeric(twist_bed_df["End"])

# Initialize an empty list to store filtered results
filtered_results = []

# Function to check overlap
def find_overlaps(chunk, twist_bed):
    overlaps = []
    for _, row in chunk.iterrows():
        overlapping = twist_bed[
            (twist_bed["Chr"] == row["Chr"]) &
            (twist_bed["Start"] < row["End"]) &
            (twist_bed["End"] > row["Start"])
        ]
        if not overlapping.empty:
            overlaps.append(row)
    return overlaps

# Read tumor plasma file in chunks for memory efficiency
chunk_size = 5000
with pd.read_csv(tumor_plasma_file, sep="\t", chunksize=chunk_size) as reader:
    for chunk in reader:
        # Ensure necessary columns are numeric
        chunk["Start"] = pd.to_numeric(chunk["Start"], errors="coerce")
        chunk["End"] = pd.to_numeric(chunk["End"], errors="coerce")
        
        # Filter overlapping regions
        filtered_chunk = find_overlaps(chunk, twist_bed_df)
        filtered_results.extend(filtered_chunk)

# Convert results to a DataFrame
filtered_df = pd.DataFrame(filtered_results)

# Save the filtered DataFrame
filtered_df.to_csv(output_file, index=False)

print(f"Filtered overlapping regions saved to {output_file}")
```

file: [wgs_wes_bench.R](./scripts/cna/wgs_wes_bench.R)

```
/hpcnfs/scratch/DIMA/delcorvo/liquid_biobsy_project/wes/scripts/wgs_wes_bench.R

# example of usage

cd /hpcnfs/scratch/DIMA/delcorvo/liquid_biobsy_project

Rscript --vanilla wes/scripts/wgs_wes_bench.R \
wes/dragen_call \
low_pass_wgs/cna_comparison/CNAs_calls.FILT.txt \
wes_wgs_comparison
```

* Results:

  - [wes_wgs_benchmarking](./results/cna/wes_bench/)
 
  ```
  /hpcnfs/scratch/DIMA/delcorvo/liquid_biobsy_project/wes_wgs_comparison
  ```
  
