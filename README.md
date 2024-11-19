# Liquid Biopsy Project
edited by Marcello Del Corvo

This repo provides details to understand and  reproduce all the analyses done for the <i> Pelicci G. et al. </i>  liquid biopsy project

## Contents
- [Contents](#contents)
- [Data](#data)
- [Code](#code)
- [Analysis](#analysis)

## Data
Raw data used for the analyses are stored on the cluster. 

We performed low-pass WGS analysis of plasma (cfDNA + evDNA) from ***6 GBM patients*** and their corresponding tumor samples and healty control (buffy coat).

With exception of ***C2-57*** patient, all samples have been also sequencing with WES for tumor and control tissue

<img src="img/samples_table.png" width="500" />

**1. Low-pass WGS (Plasma-Tumor-Control samples)**
```
cd /hpcnfs/techunits/genomics/PublicData/PelicciG/sfaletti/FASTQ/230919_A00302_0562_BHJ7C3DRX3
```
**2. WES (Tumor-Control samples)**
```

```
Results (both data and plots) are also stored on the cluster, or can be easily generated by running the appropriate scripts.

## Code
All code used for the analyses is provided at: 

* [https://github.com/mdelcorvo/DeSeq-Free](https://github.com/mdelcorvo/DeSeq-Free) (branches: master, see README for usage details)

* [scripts](./scripts/) directory

## Analysis

At the surface level, the analyses can be broadly grouped into these sections:

1) Low-pass WGS analysis
* Input data:
* Code:
* Results:

2) WES analysis
* Input data:
* Code:
* Results:
