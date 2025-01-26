# Linear-modelling-using-LIMMA

**Understanding linear modelling in differential methylation analysis using the limma method**  

As a part of my Master Thesis project, I have analysed Methylation datasets, and a very common Research question that I have come across is “Which of the cpg sites are differentially methylated across the different phenotypic groups in the dataset?” In this blog I will be sharing my understanding of the limma method in R using the case of answering my Research question.  

LIMMA is a library for the analysis of microarray data, especially the use of linear models for analysing designed experiments and the assessment of differential methylation. Limma operates on a matrix of expression or methylation values. Therfore, in order to answer the above Research question, we need to navigate through a couple of data pre-processing steps (which will be the content of another write-up) and transform our dataset into the format of a matrix, each cell containing M-values corresponding to a specific cpg probe for a particular sample.  

<img width="308" alt="Screenshot 2025-01-26 at 14 13 16" src="https://github.com/user-attachments/assets/0851e703-15b9-461d-87bc-fdf27202d368" />

Before continuing with the methodology, I will provide a bit of introduction about the M-values used by limma for statistical modelling.
The quantitative of the methylation at any cpg site is obtained using beta values and M values. Every cpg site used in the assay has one methylated and one unmethylated intensity recorded against them.  

M-values take the ratio of the methylated and unmethylated intensities at a cpg site and **apply a log transformation**. This log transformation results in the M-values being **unbounded**. A value of around 0 means that the methylated and unmethylated signals are close to equal. **A positive M-value means that the methylated signal is higher, while a negative M-value means that the unmethylated signal is higher.**  
<img width="260" alt="Screenshot 2025-01-26 at 14 48 54" src="https://github.com/user-attachments/assets/1f1e12b8-bdaf-49d6-b223-6fbc1936eff7" />





