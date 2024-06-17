# About 
---
# short_ds_cfDNA: DNA footprinting using short double-stranded cell-free DNA from plasma

## Introduction
This repository features custom scripts and workflows for the study described in the manuscript: [A novel approach for DNA footprinting using short double-stranded cell-free DNA from plasma](https://www.biorxiv.org/content/10.1101/2024.02.09.579588v1.full). It provides tools for analyzing short double-stranded cfDNA fragments enriched from human plasma and reproduction of results in the manuscript.

## Abstract
We introduce a method for enriching short double-stranded cfDNA, approximately 40 base pairs in length, from cfDNA for high-throughput DNA sequencing. These cfDNA fragments, enriched at regulatory genomic loci and transcription factor binding sites, enable genome-wide DNA footprinting from liquid biopsies. Our findings reveal significant enrichment of transcription factor motifs and correlations with DNA methylation, histone modifications, and gene transcription, as well as the basic demonstration of the diagnostic potential of short double-stranded cfDNA in disease differentiation via liquid biopsies.

## Contents
The repository is organized as follows:
1. **Data Processing:** Scripts for preprocessing, mapping, and postprocessing of short cfDNA sequencing data. [View](https://github.com/janmueller17/short_ds_cfDNA/tree/main/map_and_process) - "Sequencing data processing"
2. **Peak Calling:** Methods for identifying significant peaks in DNA footprints. [View](https://github.com/janmueller17/short_ds_cfDNA/tree/main/peak_calling) - "Peak calling"
3. **Motif Enrichment:** Analysis of transcription factor motifs within short double-stranded cfDNA. [View](https://github.com/janmueller17/short_ds_cfDNA/tree/main/transcription_factor_motif_analyses) - "Transcription factor motif enrichment analysis"
4. **NFR Calling:** Identification of nucleosome-free regions from regular cfDNA. [View](https://github.com/janmueller17/short_ds_cfDNA/tree/main/nucleosome_free_regions) - "Nucleosome-free region calling"
5. **Visualization:** Scripts for generating cumulative coverage profiles and heatmaps. [View](https://github.com/janmueller17/short_ds_cfDNA/tree/main/average_coverage_analyses) - "Average coverage profiles and heatmaps"
6. **Composite Analysis:** Examining the effects of DNA methylation and short double-stranded cfDNA signals on gene transcription. [View](https://github.com/janmueller17/short_ds_cfDNA/tree/main/composite_effects_short_methyl) - "Composite effect of DNA methylation and DNA footprint signals on gene transcription"

## Usage
Ensure all dependencies are installed as specified in each script's header. All analyses are designed for execution on Linux systems.

## License
This project is available under the Apache-2.0 License. For more information, see the [LICENSE](./LICENSE) file.

## Contact
For further inquiries, please contact us via github or our [website](https://www.igb.fraunhofer.de/en/research/in-vitro-diagnostics.html).

