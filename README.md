
# Mouse Liver Microarray Analysis — GSE54723

## Overview

This project presents an exploratory analysis of public mouse liver microarray data from GEO accession GSE54723. It compares gene-expression profiles from mice fed sucrose-containing and palatinose-containing diets for 22 weeks.

The workflow uses R and Bioconductor to retrieve the dataset, inspect sample-level patterns, perform differential-expression analysis, annotate probes, and generate reproducible tables and figures.

## Research question

How do liver gene-expression profiles differ between the sucrose and palatinose dietary groups?

## Dataset

- GEO accession: GSE54723
- Organism: *Mus musculus*
- Platform: GPL6885 — Illumina MouseRef-8 v2.0
- Samples: 14 total (7 sucrose, 7 palatinose)
- Expression features: 25,697 probes

Data source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54723

## Methods

1. Retrieve the GEO Series Matrix using `GEOquery`.
2. Inspect sample metadata and assign dietary groups.
3. Explore expression distributions and sample-level patterns.
4. Apply a log2 transformation to the positive Series Matrix values, consistent with GEO2R's automatic transformation decision for this dataset.
5. Visualize samples using principal component analysis (PCA) and a sample-distance heatmap.
6. Test differential expression using `limma`.
7. Adjust p-values using the Benjamini–Hochberg procedure.
8. Annotate probes using `AnnotationDbi` and `illuminaMousev2.db`.
9. Export results, figures, and session information.

## Results

- Samples analyzed: 14
- Probes tested: 25,697
- Probes with adjusted p-value < 0.05: 0
- SOCS2 probe: `ILMN_2628178`
- SOCS2 log2 fold change (palatinose vs sucrose): −1.134
- SOCS2 raw p-value: 0.00516
- SOCS2 adjusted p-value: 0.7445

No probe passed the chosen false discovery rate (FDR) threshold of 0.05. SOCS2 showed a negative estimated log-fold change, but this result was not statistically significant after multiple-testing correction. It should be treated as an exploratory observation, not a confirmed biological finding.

## Quality control and sensitivity analysis

The first two principal components explained approximately 25.06% of the variance. The PCA showed substantial overlap between the dietary groups.

Sample `GSM1322815` had an unusual position in the PCA plot. Its overall expression summaries were similar to those of the other samples, and its expression profile remained highly correlated with the other sucrose samples.

A sensitivity analysis excluding this sample also found zero probes passing FDR < 0.05. The smallest adjusted p-value was 0.7445 in the full analysis and 0.3092 in the reduced analysis.

The full 14-sample analysis remains the primary analysis. The sensitivity analysis does not, by itself, justify excluding the sample.

## Annotation handling

Probe annotations were obtained using `illuminaMousev2.db`.

The annotated differential-expression table contains one row per tested probe. A separate file, `probe_annotation_all_mappings.csv`, preserves the alternative annotation records returned by the annotation database.

Some probes have multiple possible gene-symbol mappings. The one-row-per-probe table uses a single annotation record for convenient reporting; this does not establish that the selected mapping is definitive for every ambiguous probe. Interpret gene-level results with this limitation in mind.

## Output files

### Figures — `results/figures/`

- `01_expression_distribution.png`
- `02_log2_expression_distribution.png`
- `03_pca_samples.png`
- `04_sample_distance_heatmap.png`
- `05_volcano_plot.png`

### Tables — `results/tables/`

- `analysis_summary.csv`
- `sample_metadata.csv`
- `differential_expression_all_probes.csv`
- `differential_expression_annotated.csv`
- `probe_annotation_all_mappings.csv`

### Other files

- `scripts/analysis.R` — analysis script
- `sessionInfo.txt` — R and package session information

## Reproducibility

Use R with the required packages installed. Set the working directory to the project root, then run:

```r
source("scripts/analysis.R")
```

The script retrieves the GEO Series Matrix, performs the analysis, and regenerates the figures and result tables. An internet connection may be required to retrieve data if no cached copy is available.

The analysis was tested with R 4.6.1 and Bioconductor 3.23 in the current working environment. Package versions may differ in other environments; consult `sessionInfo.txt`.

## Limitations

- GEO metadata reports background subtraction, offset adjustment, and quantile normalization, but the complete original preprocessing history has not been independently established.
- The log2 transformation matches GEO2R's automatic decision. This agreement does not independently validate every aspect of the original processing.
- The sample size is small: seven samples per dietary group.
- No probes passed the FDR significance threshold of 0.05.
- The PCA and sensitivity analysis are exploratory and do not establish the biological cause of sample-level differences.
- Some probe-to-gene annotations are ambiguous.

These results should not be presented as validated biological discoveries or as proof that the diets have no biological effects.

## Tools

R, GEOquery, Biobase, limma, AnnotationDbi, illuminaMousev2.db, ggplot2, and pheatmap.

## License

Add an appropriate license before redistributing this project. Confirm the terms for any included data, annotations, and other third-party materials separately.
