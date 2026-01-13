## Restoration of Ancestral Plasticity Contributes to Plastic Heterosis in the Fatty Liver of Hybrid Ducks

### Overview

This repository contains scripts and analysis pipelines used to investigate transcriptomic plasticity and heterosis in hybrid ducks

---

### Workflow

#### **01. FastQC**

FastQC quality control results for all RNA-seq samples.

---

#### **02. Parental Transcriptomic Analysis**

Analyses applied to the parental duck species, including:

- **01.gene_expression_quantification.sh**  
  Quantification of gene expression in parental species under two feeding conditions.

- **02.orthologous_genes_identification.sh**  
  Identification of orthologous genes among Peking duck, Muscovy duck, and goose.

- **03.ancestral_plasticity_inferences.sh**  
  Estimation of evolutionary rates of the hepatic transcriptome and inference of ancestral gene expression levels.

- **04.adaptive_nonadaptive_plasticity_inferences.R**  
  Identification of genes exhibiting adaptive and nonadaptive plasticity based on gene expression patterns.

---

#### **03. Hybrid Transcriptomic Analysis**

Analyses applied to hybrid ducks, including:

- **01.gene_expression_quantification.sh**  
  Quantification of gene expression in hybrid ducks under two feeding conditions.

- **02.calculate_heterosis_percentage.R**  
  Estimation of the degree of heterosis associated with fatty liver development in hybrid ducks.

- **03.inheritance_modes_classification.R**  
  Classification of gene expression inheritance modes in hybrid ducks.

- **04.WGCNA.R**  
  Weighted Gene Co-expression Network Analysis (WGCNA) based on allele expression in hybrid ducks.

---

#### **04. Simulation**

Simulation-based analyses, including:

- **01.build_pseudo_hybrid.sh**  
  Construction of pseudo-hybrid ducks by randomly downsampling parental RNA-seq reads to 50%.

- **02.pseudo_hybrid_align.sh**  
  Alignment of pseudo-hybrid RNA-seq reads to the pseudo-haplotype-resolved hybrid genome.

- **03.haplotype_levels_stat.sh**  
  Estimation of the proportion of reads correctly assigned to parental haplotypes.

- **04.allele_levels_stat.sh**  
  Estimation of the proportion of reads correctly assigned to parental alleles.
