# Supplementary-Material-for-BiRS
Supplementary codes and tables for BiRS

## GWAS results of BiRS and Q-SCAN
- *GWAS-BiRS*: GWAS results of BiRS
- *GWAS-Scan*: GWAS results of Q-SCAN

## Code for BiRS
### Overview

1. Package *BiRS* contains the main functions of BiRS algorithm and DCF test. One can install the package *BiRS* through the file ***BiRS.tar.gz***. The main functions are as follows:
- *SigDetDCF*: R complementation of BiRS-DCF
- *SigDetCL*: R complementation of BiRS-CL
- *SigScanQ*: R complementation of Q-SCAN under the framework of two sample mean test
- *SigScanDCF*: R complementation of DCF-SCAN
- *SigScanCL*: R complementation of CL-SCAN
- *XFTest*: R complementation of DCF test
- *DCFM*: R complementation for generating the test statistic vector and critical value matrix when there is missing value in sample matrix appears

2. Directory *Simulation_BiRS* contains the main functions for the simulations. Specifically, 
- The R script ***M-dependence/Normal_ESBA*** contains the code for generating data for ''M-dependence with equal covariance and non-decay signal'' setting, BiRS-DCF detection, BiRS-CL detection and Q-SCAN detection.
- The R script ***M-dependence/Normal_NSBA*** contains the code for generating data for ''M-dependence with unequal covariance and non-decay signal'' setting, BiRS-DCF detection, BiRS-CL detection and Q-SCAN detection. 
- The R script ***M-dependence/Normal_ESBA_Decay*** conntains the code for generating data for ''M-dependence with equal covariance and decay signal'' setting, BiRS-DCF detection, BiRS-CL detection and Q-SCAN detection.
- The R script ***M-dependence/Normal_NSBA_Decay*** contains the code for generating data for ''M-dependence with unequal covariance and decay signal'' setting, BiRS-DCF detection, BiRS-CL detection and Q-SCAN detection. 
- The R scripts ***M-dependence/Describe*** and ***M-dependence/Describe_Decay*** contain the code for generating description information of the detection results under M-dependence covariance structure. (The RDATA file need to generate using the previous 4 detection functions and the RDATA file ***Mulist_NDi, i = 1, ..., 5*** is the same as the file ***Mulisti, i = 16, ..., 20***)
- The R scripts ***M-dependence/Organize*** and ***M-dependence/Organize_Decay*** contain the code for generating figures of the detection results under M-dependence covariance structure

- The R script ***Weak-dependence/Normal_ESBA*** contains the code for generating data for ''Weak-dependence with equal covariance and non-decay signal'' setting, BiRS-DCF detection, BiRS-CL detection and Q-SCAN detection.
- The R script ***Weak-dependence/Normal_NSBA*** contains the code for generating data for ''Weak-dependence with unequal covariance and non-decay signal'' setting, BiRS-DCF detection, BiRS-CL detection and Q-SCAN detection. 
- The R script ***Weak-dependence/Normal_ESBA_Decay*** conntains the code for generating data for ''Weak-dependence with equal covariance and decay signal'' setting, BiRS-DCF detection, BiRS-CL detection and Q-SCAN detection.
- The R script ***Weak-dependence/Normal_NSBA_Decay*** contains the code for generating data for ''Weak-dependence with unequal covariance and decay signal'' setting, BiRS-DCF detection, BiRS-CL detection and Q-SCAN detection. 
- The R scripts ***Weak-dependence/Describe*** and ***Weak-dependence/Describe_Decay*** contain the code for generating description information of the detection results under Weak-dependence covariance structure. (The RDATA file need to generate using the previous 4 detection functions and the RDATA file ***Mulist_NDi, i = 1, ..., 5*** is the same as the file ***Mulisti, i = 16, ..., 20***)
- The R scripts ***Weak-dependence/Organize*** and ***Weak-dependence/Organize_Decay*** contain the code for generating figures of the detection results under Weak-dependence covariance structure
- The R script ***Speed-Test*** contains the code for calculating computational times of BiRS-DCF, DCF-SCAN and Q-SCAN

3. Directory *Application_Code* contains the main functions for quality control and performing GWAS on UK Biobank data. Specifically, 
- The R scripts ***CTtraits*** and ***covariate*** generate samples and covariates we needed from the main dataset of UK Biobank.
- The R scripts ***SampleQC*** and ***VariantQC*** perform the sample and variant quality control in SI Appendix after generating .vcf file for each chromosome using PLINK 2.0 (the quality control on HWE and MAF performed at the same time as the .vcf file be generated)
- The R script ***SampleGLM*** generate data for Q-SCAN to perform GWAS.
- The R script ***GenerateFinalVariant*** generate the information for the selected variants for performing GWAS.
- The R scripts ***BiRSGenomeSize***, ***QScanGenomeSize*** and ***BiCLGenomeSize*** perform permutation tests of three methods.
- The R scripts ***BiRSGenome*** and ***QScanGenome*** perform GWAS on C50 Malignant neoplasms of breast using UK Biobank data.
