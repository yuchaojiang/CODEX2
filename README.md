<img border="0" src="http://bioconductor.org/shields/availability/release/CODEX.svg"/> <img border="0" src="http://bioconductor.org/shields/downloads/CODEX.svg"/> <img border="0" src="http://bioconductor.org/shields/build/devel/bioc/CODEX.svg"/> <img border="0" src="http://bioconductor.org/shields/years-in-bioc/CODEX.svg"/>


# CODEX
A Normalization and Copy Number Variation Detection Method for Whole Exome Sequencing


## Author
Yuchao Jiang, Nancy R. Zhang

## Maintainer
Yuchao Jiang <yuchaoj@upenn.edu>

## Description
A normalization and copy number variation calling procedure for
whole exome DNA sequencing data. CODEX relies on the availability of 
multiple samples processed using the same sequencing pipeline for 
normalization, and does not require matched controls. The normalization 
model in CODEX includes terms that specifically remove biases due to GC 
content, exon length and targeting and amplification efficiency, and latent
systemic artifacts. CODEX also includes a Poisson likelihood-based recursive
segmentation procedure that explicitly models the count-based exome 
sequencing data.


## Installation
* Install the current release from Bioconductor
```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("CODEX")
biocLite("WES.1KG.WUGSC") # toy dataset
```

* Install the devel version from GitHub
```r
install.packages("devtools")
library(devtools)
install_github("yuchaojiang/CODEX/package")
```

## Demo code & Vignettes
* [Demo code](http://www.bioconductor.org/packages/devel/bioc/vignettes/CODEX/inst/doc/CODEX_vignettes.R)
* [Vignettes](http://www.bioconductor.org/packages/devel/bioc/vignettes/CODEX/inst/doc/CODEX_vignettes.pdf)
* [Bioconductor page](http://www.bioconductor.org/packages/release/bioc/html/CODEX.html)
* [Bioconductor download stats](http://bioconductor.org/packages/stats/bioc/CODEX/)


## Google user group (Q&A)
https://groups.google.com/d/forum/codex_wes_cnv


## Citation
Jiang, Y., Oldridge, D.A., Diskin, S.J. and Zhang, N.R., 2015. CODEX: a normalization and copy number variation detection method for whole exome sequencing. *Nucleic acids research*, 43(6), pp.e39-e39. [[html](http://nar.oxfordjournals.org/content/43/6/e39), [pdf](http://nar.oxfordjournals.org/content/43/6/e39.full.pdf+html)]


## IMPORTANT: CODEX for cancer genomics
When apply CODEX to whole-exome sequencing and targeted sequencing in cancer genomics:
* In normalization step, use *normalize2(...)* function if there are normal samples and specify the index of the normal samples as a numerical vector in the *normal_index* argument;
* In segmentation step, use **fractional** mode for somatic CNA detection (cancer is heterogenous) and **interger** mode for germline CNV detection (you will get CNV calls in your blood samples, which are germline).
* For segmentation with paired tumor-normal experimental design, a modified CBS (circular binary segmentation) algorithm can be adopted, which ultilizes the pair information. Refer to the paired_tumor_normal_segmentation folder for code (not actively updated/maintained). Note that, from our experience, the default segmentation by CODEX (not using the pair information) does not make much difference. Normalization is the first order effect in WES study design.

## CODEX for targeted sequencing
We've adapted CODEX for targeted sequencing. Instead of normalizing and segmenting each chromosome separately, for targeted sequencing, we combine all targets across the genome to perform normalization, followed by segmentation within each gene. Refer to codes below (need to source segment_targeted.R for gene-based segmentation).
* [codex_targeted.R](https://github.com/yuchaojiang/CODEX/blob/master/targeted_sequencing/codex_targeted.R)
* [segment_targeted.R](https://github.com/yuchaojiang/CODEX/blob/master/targeted_sequencing/segment_targeted.R)

## Visualization by IGV
One can load CODEX's CNV calling results into [IGV](http://www.broadinstitute.org/igv/) for visualization by generating a tab-delimited seg file for each sample. Below is a sample code that we use in our daily practice -- for each sample, a *.seg.txt file is generated with six columns and header 'Sample', 'Chromosome','Start','End','Num_Probes','Segment_Mean', which correspond to sample name, chromosome, CNV start bp, CNV end bp, number of exonic targets, and log ratio of raw (i.e. observed) depths of coverage versus normalized (i.e. expected) coverage (deletion has a negative log ratio, duplication has a positive log ratio, copy-neutral region has a log ratio around 0).
* [CODEX_IGV.R](https://github.com/yuchaojiang/CODEX/blob/master/IGV_visualization/CODEX_IGV.R)

## CODEX for hg38?
CODEX by default is for hg19 reference. It can be adapted to hg38: only the calculations of GC content and mappability need to be changed; to get coverage for exons across samples stays the same (make sure that the exonic targets in the bed file are also in hg38 coordinates). To calculte GC content in hg38, you need to download the hg38 reference from [Bioconductor](https://bioconductor.org/packages/BSgenome.Hsapiens.UCSC.hg38/). Then, after loading CODEX, load the hg38 reference package so that the Hsapiens (hg19 by CODEX's default) is masked to hg38.

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg38")

library(CODEX)
library(BSgenome.Hsapiens.UCSC.hg38)
# The following object is masked from ‘package:BSgenome.Hsapiens.UCSC.hg19’:  Hsapiens
```
To calculate mappability for hg38 is a bit more complicated and time-consuming. For CODEX, we pre-compute mappabilities for all hg19 exons and store them as part of the package. For hg38, there are two workarounds: 1) set all mappability to 1 using mapp=rep(1,length(gc)) since mappability is only used in the QC step to filter out exons with low mappability and thus should not affect the final output too much; 2) adopt QC procedures based on annotation results, e.g., filter out all exons within segmental duplication regions, which generally have low mappability.

Note that CODEX can also be adapted to the mouse genome, see below.

## CODEX for mouse genome
CODEX can be applied to WES of the mouse genome. The library for the mm10 mouse genome sequencing needs to be loaded: 
* [BSgenome.Mmusculus.UCSC.mm10](http://bioconductor.org/packages/release/data/annotation/html/BSgenome.Mmusculus.UCSC.mm10.html).

The GC content and the mappability can be obtained from the code below with minor adaptations for the mouse genome:
* [GC content](https://github.com/yuchaojiang/CODEX/blob/master/mouse/getgc.R)
* [Mappability pre-calculation](https://github.com/yuchaojiang/CODEX/blob/master/mouse/mapp.R) (Note: This step can be computationally extensive and thus parallel computing is recommended. For CODEX in its default setting, the mappability for exonic targets in human h19 assembly is pre-computed and stored as part of the package).
* [Mappability](https://github.com/yuchaojiang/CODEX/blob/master/mouse/getmapp.R)

## Common questions
  * How many samples does CODEX need? Should I separately run samples from different batches?
    
    We have applied CODEX to data sets of sample size ranging from 30 to 500. Yes, samples from different batches are highly recommended to run separately. The Poisson latent factor can presumably capture the batch effects but if additional knowledge is available beforehand, it should be ultilized. If batch information is not available, sometimes we refer to the header within the bam files.
    
  * Error in glm.fit?
    
    Yes, we are aware that sometimes the normalize() function leads to error in glm.fit. CODEX adopts an iterative estimation procedure to estimate the Poisson latent factors via Poisson glm, the exon-specific bias by taking the median across all samples, and the GC content bias by fitting a non-parametric smooth.spline. We did our best to make sure that the iteration/estimation runs properly, yet sometimes the Poisson glm function in R still fails to converge due to: (1) extreme heterogeneity in the data (i.e., the data is just too noisy or the samples are from multiple batches); (2) the number of Poisson latent factors to estimate is too large. See question below.
    
  * What is the range of K? Which one is optimal?
  
    Based on our experience, very rarely do we run CODEX with greater than 10 latent factors (i.e., K = 1:10 suffices for most, if not all, datasets we have). The larger the K is, the longer the estimation takes. Also, refer to the previous question regarding potential pitfalls with a large value of K.
    
    CODEX includes three statistical metrics to help the users more wisely choose the optimal K: AIC, BIC, and residual variance. A pdf plot is automatically generated by the choiceofK() function. Sometimes the optimal value based on the metrics are not clear. In this case, we recommend a sanity check by focusing on known positive/negative controls and visualizing the normalization/segmentation results. For normalize2() which specifies the normal samples, the effect on different optimal K values diminishes since only the normal samples are used to estimate the exon-wise biases and latent factors.
    
    
    
    
    
