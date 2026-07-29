R code of a predictive model for dominant influenza A serotype by quantifying genetic mutation activities


Overview
--------
This repository provides an R implementation of the dominant influenza A serotype prediction model described in:

  Lou J, Zhao S, Cao L, Chong MKC, Chan RWY, Chan PKS, Zee BCY, Yeoh EK, Wang MH. 
  Predicting the dominant influenza A serotype by quantifying mutation activities. 
  Int J Infect Dis. 2020;100:255-257. doi:10.1016/j.ijid.2020.08.053

The model quantifies genetic mutation activity (g-measure) from hemagglutinin (HA) protein sequences of A/H1N1 and A/H3N2, and infersthe dominant probability of each influenza A subtypesby in the upcoming influenza season using a log-odds linear regression framework. A positive relationship between relative genetic activity and relative epidemic prevalence is established. A 5-fold cross-validation is used to assess model discrimination performance.


Input and Output
--------
The prediction model requires data to be collected separately for each target region. Here, the R script is conducted using toydata in Singapore as example.

Main Input:
  1. HA protein genetic sequences (FASTA format)
   - A/H1N1 HA sequences from Oct 2009 to Sep 2020 in Singapore (file name: SGH1_HA_sequence.fasta)
   - A/H3N2 HA sequences (file name: SGH3_HA_sequence.fasta)
   - Source: GISAID flu database
  2. Laboratory-based influenza surveillance data (file name: SG_seropositivity_rate.csv)
     - total_speciman: number of specimens tested
     - h1n1pdm09_positive: number of A/H1N1pdm09 positive specimens
     - h3n2_positive: number of A/H3N2 positive specimens
     - Source: WHO GISRS FluNet

Main Output
 1. Regression model: log(H3_rate / H1_rate) = beta0 + beta1 * log(H3_g / H1_g)
 2. Model summary: coefficients, R-squared, p-value
 3. In-sample confusion matrix (predicted vs. observed dominant serotype)
 4. Scatter plot with regression line and 95% confidence interval (PDF)
 5. Intermediate outputs (optional): sequence tables, site-wise prevalence tables, g-measure tables


R Dependencies
------------
R (>= 4.1.3) or newer version is required. To download R, please see: https://www.r-project.org/

Required R packages:
  Biostrings (2.64.1) : https://bioconductor.org/packages/release/bioc/html/Biostrings.html
  lubridate  (1.8.0)  : https://www.rdocumentation.org/packages/lubridate/versions/1.8.0
  plyr       (1.8.7)  : https://www.rdocumentation.org/packages/plyr/versions/1.8.7
  dplyr      (1.0.0)  : https://www.rdocumentation.org/packages/dplyr
  tibble     (3.0.0)  : https://www.rdocumentation.org/packages/tibble

Install packages in R:
  install.packages(c("lubridate", "plyr", "dplyr", "tibble"))
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("Biostrings")


Step-by-step Description and Code Execution
------------------------
Step 1. Preparing protein sequences: the R function 'PrepareSeq' reads aligned HA protein sequences in FASTA format and converts the data to a csv-compatible data frame for further analysis. For each sequence, it parses the accession number, isolate name, collection date, and assigns a flu season label. Sequence nomenclature: flu season 2011/12 is labeled as 2012.

Step 2. Calculating site-wise amino acid prevalence: the function 'SitePrev' calculates the yearly proportion of each amino acid at each site position of the HA protein. This produces a site-wise prevalence table through time, which captures the evolutionary dynamics at each residue.

Step 3. Calculating g-measure (mutational activity metric): the g-measure quantifies the overall level of key mutational activity through time and reflects both the scale and the prevalence of key mutations.
Reference for g-measure: 
  Wang MH, Lou J, Cao L, Zhao S, Chan RW, Chan PK, Chan MC, Chong MK, Wu WK, Wei Y, Zhang H, Zee BC, Yeoh EK. 
  Characterization of key amino acid substitutions and dynamics of the influenza virus H3N2 hemagglutinin. 
  J Infect. 2021;83(6):671-677. doi:10.1016/j.jinf.2021.09.026

Step 4. Regression fitting: this step fits a linear regression model that associates the relative genetic activity with the relative epidemic prevalence (Equation S6.1 in the paper): log(Pr_H3 / Pr_H1) = beta0 + beta1 * log(g_H3 / g_H1), where:
    Pr_H3, Pr_H1 = positive rates of H3N2 and H1N1, respectively
    g_H3, g_H1 = g-measure values of H3N2 and H1N1, respectively
The predicted dominant serotype for each season is determined by the sign of the predicted log odds.

Step 5. (Optional) Scatter plot: generates a scatter plot of log(H3/H1 positive rate) versus log(H3/H1 g-measure) with the fitted regression line and its 95% confidence interval. Data points are colored by the observed dominant serotype, and labeled by season.

To run the toy example:
  1. Set the working directory to the toydata/ folder
  2. Source the main script: source("../code/flu subtype prediction model.R")
  3. All intermediate and final outputs will be generated in toydata/toydata_result/

Note: Because Singapore represents a single region with a limited number of seasons (approximately 8-10 valid seasons after filtering), the model fit from this toy example may not achieve optimal performance. The published study used data from 7 countries/regions (New York, California, United Kingdom, Hong Kong, Singapore, Thailand, and New Zealand) spanning the 2009-2019 flu seasons, yielding 76 region-season samples and substantially better fitting results (AUC = 0.78, sensitivity = 0.84, precision = 0.79, R-squared = 0.42).


Data Availability & Acknowledgements
-------------------------------------
Influenza genetic sequences were obtained from:
  Global Initiative on Sharing All Influenza Data (GISAID): https://www.gisaid.org/

Laboratory surveillance data were obtained from:
  WHO Global Influenza Surveillance and Response System (GISRS): https://www.who.int/initiatives/global-influenza-surveillance-and-response-system

Detailed data sources for each region are listed in Supplementary Materials S1 of the paper.
We thank the contributions of all health care workers, scientists, the GISAID team, and the submitting and originating laboratories.


Citation
--------
  Lou J, Zhao S, Cao L, Chong MKC, Chan RWY, Chan PKS, Zee BCY, Yeoh EK, Wang MH. 
  Predicting the dominant influenza A serotype by quantifying mutation activities. 
  Int J Infect Dis. 2020;100:255-257. doi:10.1016/j.ijid.2020.08.053


Contact
-------
For questions or feedback, please contact the corresponding author:
  Maggie H Wang (maggiew@cuhk.edu.hk)
  JC School of Public Health and Primary Care, Chinese University of Hong Kong
