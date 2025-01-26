# Linear-modelling-using-LIMMA

**Understanding linear modelling in differential methylation analysis using the limma method**  

As a part of my Master Thesis project, I have analysed Methylation datasets, and a very common Research question that I have come across is “Which of the cpg sites are differentially methylated across the different phenotypic groups in the dataset?” In this blog I will be sharing my understanding of the limma method in R using the case of answering my Research question.  

LIMMA is a library for the analysis of microarray data, especially the use of linear models for analysing designed experiments and the assessment of differential methylation. Limma operates on a matrix of expression or methylation values. Therfore, in order to answer the above Research question, we need to navigate through a couple of data pre-processing steps (which will be the content of another write-up) and transform our dataset into the format of a matrix, each cell containing M-values corresponding to a specific cpg probe for a particular sample.

