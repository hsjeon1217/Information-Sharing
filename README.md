# DEG Information Sharing

## Improving genomic findings through information sharing between a pilot and main study

This is a statistical framework designed to identify differentially expressed genes (DEGs) while controlling the false discovery rate (FDR) using two p-value vectors from two independent studies. First, it is assumed that study A has a lower detection power than study B, with study A serving as a pilot study for study B. In addition, we also assume the same DEG across studies. If the DEGs are not the same, the FDR is controlled in the circumstance of detecting genes that are DE in at least one study, within a meta-analysis context. The implementation of the method is illustrated in Tutorial.R. If you have questions, please contact the author.
