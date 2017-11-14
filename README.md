# CODEX2
Full-spectrum copy number variation detection by high-throughput DNA sequencing

## Author
Yuchao Jiang, Nancy R. Zhang

## Maintainer
Yuchao Jiang <yuchaoj@email.unc.edu>

## Description
High-throughput DNA sequencing enables detection of copy number variations (CNVs) on the genome-wide scale with finer resolution compared to array-based methods, but suffers from biases and artifacts that lead to false discoveries and low sensitivity. We describe CODEX2, a statistical framework for full-spectrum CNV profiling that is sensitive for variants with both common and rare population frequencies and that is applicable to study designs with and without negative control samples. We demonstrate and evaluate CODEX2 on whole-exome and targeted sequencing data, where biases are the most prominent. CODEX2 outperforms existing methods and, in particular, significantly improves sensitivity for common CNVs.


## Installation
```r
# Install dependent packages first
source("https://bioconductor.org/biocLite.R")
biocLite("GenomeInfoDb")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
biocLite("WES.1KG.WUGSC")

install.packages("devtools")
library(devtools)
install_github("yuchaojiang/CODEX/package")
install_github("yuchaojiang/CODEX2/package")
```


## Questions?
If you've encountered problems or bugs running CODEX2, please report directly [here](https://github.com/yuchaojiang/CODEX2/issues) using the Issues tab by GitHub.


## Citation
Yuchao Jiang, Ioannis N. Anastopoulos, Katherine L. Nathanson, Nancy R. Zhang, 2017. Full-spectrum copy number variation detection by high-throughput DNA sequencing. *Submitted* ([bioRxiv](https://www.biorxiv.org/content/early/2017/10/30/211698)).


## Running CODEX2
The figure below illustrates the two experimental designs for which CODEX2 can be applied: (i) case-control design with a group of negative control samples, where the goal is to detect CNVs disproportionately present in the ‘cases’ versus the ‘controls’; and (ii) detection of all CNVs present in all samples design, such as in the Exome Aggregation Consortium. The key innovation in CODEX2 is the usage of *negative control genome regions* in a genome-wide latent factor model for sample- and position-specific background correction, and the utilization of *negative control samples*, under a case-control design, to further improve background bias estimation under this model. The negative control genome regions defined by CODEX2 are regions that do not harbor common CNVs, but that are still allowed to harbor rare CNVs, and can be constructed from existing studies or learned from data.

<p align="center">
  <img src='https://github.com/yuchaojiang/CODEX2/blob/master/demo/Figure1.png' width='450' height='350'>
</p>



## IMPORTANT: CODEX2 for cancer genomics
* In segmentation step, use **fractional** mode for somatic CNA detection (cancer is heterogenous) and **interger** mode for germline CNV detection (you will get CNV calls in your blood samples, which are germline).
* For segmentation with paired tumor-normal experimental design, a modified CBS (circular binary segmentation) algorithm can be adopted, which ultilizes the pair information. Refer to the paired_tumor_normal_segmentation folder for code (not actively updated/maintained). Note that, from our experience, the default segmentation by CODEX2 (not using the pair information) does not make much difference. Normalization is the first order effect in WES study design.

## CODEX2 for targeted sequencing
We've adapted CODEX2 for targeted sequencing. Instead of normalizing and segmenting each chromosome separately, for targeted sequencing, we combine all targets across the genome to perform normalization, followed by segmentation within each gene. Refer to codes below (need to source segment_targeted.R for gene-based segmentation).
* [codex2_targeted.R](https://github.com/yuchaojiang/CODEX2/blob/master/targeted_sequencing/codex2_targeted.R)
* [segment_targeted.R](https://github.com/yuchaojiang/CODEX2/blob/master/targeted_sequencing/segment_targeted.R)

## Visualization by IGV
One can load CODEX2's CNV calling results into [IGV](http://www.broadinstitute.org/igv/) for visualization by generating a tab-delimited seg file for each sample. Below is a sample code that we use in our daily practice -- for each sample, a *.seg.txt file is generated with six columns and header 'Sample', 'Chromosome','Start','End','Num_Probes','Segment_Mean', which correspond to sample name, chromosome, CNV start bp, CNV end bp, number of exonic targets, and log ratio of raw (i.e. observed) depths of coverage versus normalized (i.e. expected) coverage (deletion has a negative log ratio, duplication has a positive log ratio, copy-neutral region has a log ratio around 0).
* [CODEX2_IGV.R](https://github.com/yuchaojiang/CODEX2/blob/master/IGV_visualization/CODEX2_IGV.R)

## CODEX2 for hg38?
CODEX2 by default is for hg19 reference. It can be adapted to hg38: only the calculations of GC content and mappability need to be changed; to get coverage for exons across samples stays the same (make sure that the exonic targets in the bed file are also in hg38 coordinates). To calculte GC content in hg38, you need to download the hg38 reference from [Bioconductor](https://bioconductor.org/packages/BSgenome.Hsapiens.UCSC.hg38/). Then, after loading CODEX2, load the hg38 reference package so that the Hsapiens (hg19 by CODEX2's default) is masked to hg38. Note that the getgc() function needs to be sourced again so that the correct version of Hsapiens is used.

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg38")

library(CODEX2)
library(BSgenome.Hsapiens.UCSC.hg38)
# The following object is masked from ‘package:BSgenome.Hsapiens.UCSC.hg19’:  Hsapiens

# Source getgc() function again so it uses the right Hsapiens
getgc = function (chr, ref) {
  if (chr == "X" | chr == "x" | chr == "chrX" | chr == "chrx") {
    chrtemp <- 23
  } else if (chr == "Y" | chr == "y" | chr == "chrY" | chr == "chry") {
    chrtemp <- 24
  } else {
    chrtemp <- as.numeric(mapSeqlevels(as.character(chr), "NCBI")[1])
  }
  if (length(chrtemp) == 0) message("Chromosome cannot be found in NCBI Homo sapiens database!")
  chrm <- unmasked(Hsapiens[[chrtemp]])
  seqs <- Views(chrm, ref)
  af <- alphabetFrequency(seqs, baseOnly = TRUE, as.prob = TRUE)
  gc <- round((af[, "G"] + af[, "C"]) * 100, 2)
  gc
}
```
To calculate mappability for hg38 is a bit more complicated and time-consuming. For CODEX2, we pre-compute mappabilities for all hg19 exons and store them as part of the package. For hg38, there are two workarounds: 1) set all mappability to 1 using mapp=rep(1,length(gc)) since mappability is only used in the QC step to filter out exons with low mappability and thus should not affect the final output too much; 2) adopt QC procedures based on annotation results, e.g., filter out all exons within segmental duplication regions, which generally have low mappability.

Note that CODEX2 can also be adapted to the mouse genome, see below.

## CODEX2 for mouse genome
CODEX2 can be applied to WES of the mouse genome. The library for the mm10 mouse genome sequencing needs to be loaded: 
* [BSgenome.Mmusculus.UCSC.mm10](http://bioconductor.org/packages/release/data/annotation/html/BSgenome.Mmusculus.UCSC.mm10.html).

The GC content and the mappability can be obtained from the code below with minor adaptations for the mouse genome:
* [GC content](https://github.com/yuchaojiang/CODEX2/blob/master/mouse/getgc.R)
* [Mappability pre-calculation](https://github.com/yuchaojiang/CODEX2/blob/master/mouse/mapp.R) (Note: This step can be computationally extensive and thus parallel computing is recommended. For CODEX2 in its default setting, the mappability for exonic targets in human h19 assembly is pre-computed and stored as part of the package).
* [Mappability](https://github.com/yuchaojiang/CODEX2/blob/master/mouse/getmapp.R)

## Common questions
  * How many samples does CODEX2 need? Should I separately run samples from different batches?
    
    We have applied CODEX2 to data sets of sample size ranging from 30 to 500. Yes, samples from different batches are highly recommended to run separately. The Poisson latent factor can presumably capture the batch effects but if additional knowledge is available beforehand, it should be ultilized. If batch information is not available, sometimes we refer to the header within the bam files.
    
  * Error in glm.fit?
    
    Yes, we are aware that sometimes the normalize() function leads to error in glm.fit. CODEX2 adopts an iterative estimation procedure to estimate the Poisson latent factors via Poisson glm, the exon-specific bias by taking the median across all samples, and the GC content bias by fitting a non-parametric smooth.spline. We did our best to make sure that the iteration/estimation runs properly, yet sometimes the Poisson glm function in R still fails to converge due to: (1) extreme heterogeneity in the data (i.e., the data is just too noisy or the samples are from multiple batches); (2) the number of Poisson latent factors to estimate is too large. See question below.
    
  * What is the range of K? Which one is optimal?
  
    Based on our experience, very rarely do we run CODEX2 with greater than 10 latent factors (i.e., K = 1:10 suffices for most, if not all, datasets we have). The larger the K is, the longer the estimation takes. Also, refer to the previous question regarding potential pitfalls with a large value of K.
    
    CODEX2 includes three statistical metrics to help the users more wisely choose the optimal K: AIC, BIC, and residual variance. A pdf plot is automatically generated by the choiceofK() function. Sometimes the optimal value based on the metrics are not clear. In this case, we recommend a sanity check by focusing on known positive/negative controls and visualizing the normalization/segmentation results. For normalize2() which specifies the normal samples, the effect on different optimal K values diminishes since only the normal samples are used to estimate the exon-wise biases and latent factors.


