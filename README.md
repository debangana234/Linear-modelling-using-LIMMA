# Linear-modelling-using-LIMMA

**Understanding linear modelling in differential methylation analysis using the limma method**  

As a part of my Master Thesis project, I have analysed Methylation datasets, and a very common Research question that I have come across is “Which of the cpg sites are differentially methylated across the different phenotypic groups in the dataset?” In this blog I will be sharing my understanding of the limma method in R using the case of answering my Research question.  

LIMMA is a library for the analysis of microarray data, especially the use of linear models for analysing designed experiments and the assessment of differential methylation. Limma operates on a matrix of expression or methylation values. Therfore, in order to answer the above Research question, we need to navigate through a couple of data pre-processing steps (which will be the content of another write-up) and transform our dataset into the format of a matrix, each cell containing M-values corresponding to a specific cpg probe for a particular sample.  

<img width="308" alt="Screenshot 2025-01-26 at 14 13 16" src="https://github.com/user-attachments/assets/0851e703-15b9-461d-87bc-fdf27202d368" />

Before continuing with the methodology, I will provide a bit of introduction about the M-values used by limma for statistical modelling.
The quantitative of the methylation at any cpg site is obtained using beta values and M values. Every cpg site used in the assay has one methylated and one unmethylated intensity recorded against them.  

M-values take the ratio of the methylated and unmethylated intensities at a cpg site and **apply a log transformation**. This log transformation results in the M-values being **unbounded**. A value of around 0 means that the methylated and unmethylated signals are close to equal. **A positive M-value means that the methylated signal is higher, while a negative M-value means that the unmethylated signal is higher.**  
<img width="260" alt="Screenshot 2025-01-26 at 14 48 54" src="https://github.com/user-attachments/assets/1f1e12b8-bdaf-49d6-b223-6fbc1936eff7" />  

**If the unmethylated signal approaches 0, the ratio becomes very large, and the M-value trends towards +∞. Similarty, the M-value trends towards −∞ when the methylated signal approaches 0.**  

**So why do we use these unbounded mathematical values that are visually not so interpretable when compared to their visually interpretable counterparts: β values?**  

Well the answer is simple: **Limma and similar statistical methods operate better on data without artificial bounds. Unbounded M-values are more suitable for linear modeling because they adhere to the assumptions of normality in statistical tests.**  

I will elaborate a bit more on this point:
Firstly, beta values are bounded within the interval 0 to 1. A value of 0 implies that the cpg site is unmethylated and a value of 1 implies that the cpg is methylated. So they tend to cluster near 0 or 1, especially for highly methylated or unmethylated regions. **This skewness creates a non-normal distribution which is not a problem in case of M-values.**  


Secondly, differences in methylation levels near the extremes (e.g., 0.95 vs. 0.90) are less pronounced compared to differences in the middle of the scale (e.g., 0.55 vs. 0.50), leading to unequal variances what we call **heteroscedasticity** in statistics. These reasons make β values unsuitable for linear models that assume homoscedasticity and normality.  


Now we have developed an understanding of the M values, let us get back to understanding how limma will help us to solve our Research Question.  

Any Research question starts with the formulation of a null hypothesis. Let the null hypothesis in this case be: **“At a particular CpG site, there is no difference in methylation levels between the compared groups (e.g., case vs. control).”**  


How does Limma test for the Null Hypothesis?  

**Limma uses linear models to test if the mean M-values at a CpG site differ significantly between the groups.**
It helps test this hypothesis by fitting linear models to the M-values at each CpG site, comparing the methylation levels across the groups .   
**The objective is to evaluate evidence against the null hypothesis by calculating a p-value for the t-statistic that quantifies the difference in methylation between groups at each site.**  
We need to break that down to understand in detail.

What is t-statistic?

In simple statistical terms, The t-statistic measures the number of standard errors the estimated coefficient is away from the hypothesized value.
Now we need to interpret this in our methylation context.

Mathematically the t-statistic at each CpG site is calculated as:







