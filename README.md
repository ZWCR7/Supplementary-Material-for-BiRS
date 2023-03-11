# Supplementary-Material-for-BiRS
Supplementary codes and tables for BiRS

## GWAS results of BiRS and Q-SCAN
- *GWAS-BiRS*: GWAS results of BiRS
- *GWAS-Scan*: GWAS results of Q-SCAN
- *GWAS-KS*: GWAS results of KnockoffScreen
- *GWAS-REGENIE*: GWAS results of regenie

## Code for BiRS
### Overview

1. Package *BiRS* contains the main functions of BiRS algorithm and DCF test. One can install the package *BiRS* through the file ***BiRS_0.10.tar.gz***.  


2. Directory *Simulation_BiRS* contains the main functions for the simulations. Specifically, 
- The R script ***KSDet*** contains the code for applying KnockoffScreen method to normal setting in simulation.
- The R script ***Speed-Test*** contains the code for calculating computational times of BiRS-DCF, DCF-SCAN, Q-SCAN and KnockoffScreen.

- The R script ***Size-Test/Size-Test*** contains the contains the code for generating data and conducting simulation for size validation of BiRS-DCF, BiRS-CL and Q-SCAN in all settings.
- The R scripts ***Size-Test/Describe_Size*** contains the code for generating tables of the size validation in all settings.

- The R script ***M-dependence/Normal_ESBA*** contains the code for generating data and conducting simulation for ''M-dependence with equal covariance and non-decay signal'' setting using BiRS-DCF, BiRS-CL, Q-SCAN and KnockoffScreen.
- The R script ***M-dependence/Normal_NSBA*** contains the code for generating data and conducting simulation for ''M-dependence with unequal covariance and non-decay signal'' setting using BiRS-DCF, BiRS-CL, Q-SCAN and KnockoffScreen.
- The R script ***M-dependence/Normal_ESBA_Decay*** conntains the code for generating data and conducting simulation for ''M-dependence with equal covariance and decay signal'' setting using BiRS-DCF, BiRS-CL, Q-SCAN and KnockoffScreen.
- The R script ***M-dependence/Normal_NSBA_Decay*** contains the code for generating data and conducting simulation for ''M-dependence with unequal covariance and decay signal'' setting using BiRS-DCF, BiRS-CL, Q-SCAN and KnockoffScreen.
- The R scripts ***M-dependence/Describe*** and ***M-dependence/Describe_Decay*** contain the code for generating description information of the detection results under M-dependence covariance structure.  
- The R scripts ***M-dependence/Organize*** and ***M-dependence/Organize_Decay*** contain the code for generating figures and tables of the detection results under M-dependence covariance structure.

- The R script ***Weak-dependence/Normal_ESBA*** contains the code for generating data and conducting simulation for ''Weak-dependence with equal covariance and non-decay signal'' setting using BiRS-DCF, BiRS-CL, Q-SCAN and KnockoffScreen.
- The R script ***Weak-dependence/Normal_NSBA*** contains the code for generating data and conducting simulation for ''Weak-dependence with unequal covariance and non-decay signal'' setting using BiRS-DCF, BiRS-CL, Q-SCAN and KnockoffScreen.
- The R script ***Weak-dependence/Normal_ESBA_Decay*** conntains the code for generating data and conducting simulation for ''Weak-dependence with equal covariance and decay signal'' setting using BiRS-DCF, BiRS-CL, Q-SCAN and KnockoffScreen.
- The R script ***Weak-dependence/Normal_NSBA_Decay*** contains the code for generating data and conducting simulation for ''Weak-dependence with unequal covariance and decay signal'' setting using BiRS-DCF, BiRS-CL, Q-SCAN and KnockoffScreen.
- The R scripts ***Weak-dependence/Describe*** and ***Weak-dependence/Describe_Decay*** contain the code for generating description information of the detection results under Weak-dependence covariance structure. 
- The R scripts ***Weak-dependence/Organize*** and ***Weak-dependence/Organize_Decay*** contain the code for generating figures and tables of the detection results under Weak-dependence covariance structure.

- The R script ***KSDetGenetic*** contains the code for applying KnockoffScreen method to genetic setting in simulation.
- The R script ***Genetic_Setting/Genetic_Setting*** contains the code for generating data and conducting simulation for genetic setting using BiRS-DCF, BiRS-CL, Q-SCAN and KnockoffScreen.
- The R scripts ***Weak-dependence/Describe*** and ***Weak-dependence/Describe_Decay*** contain the code for generating description information of the detection results under Weak-dependence covariance structure. 
- The R scripts ***Weak-dependence/Organize*** and ***Weak-dependence/Organize_Decay*** contain the code for generating figures and tables of the detection results under Weak-dependence covariance structure.

3. Directory *Application_Code* contains the main functions for quality control and performing GWAS on UK Biobank data. Specifically, 
- The R scripts ***Pretreatment/CTtraits*** and ***Pretreatment/covariate*** generate samples and covariates we needed from the main dataset of UK Biobank.
- The R scripts ***Pretreatment/SampleQC*** and ***Pretreatment/VariantQC*** perform the sample and variant quality control in SI Appendix after generating .vcf file for each chromosome using PLINK 2.0 (the quality control on HWE and MAF performed at the same time as the .vcf file be generated)
- The R script ***Pretreatment/SampleGLM*** generate data for Q-SCAN to perform GWAS.
- The R script ***Pretreatment/GenerateFinalVariant*** generate the information for the selected variants for performing GWAS.
- The R scripts ***BiRS/BiRSGenomeSize***, ***QSCAN/QScanGenomeSize*** and ***BiRS/BiCLGenomeSize*** perform permutation tests of BiRS-DCF, BiRS-CL and QSCAN.
- The R scripts ***BiRS/BiRSGenome***, ***QSCAN/QScanGenome*** and ***KnockoffScreen/KSGWAS*** perform GWAS on C50 Malignant neoplasms of breast using UK Biobank data by BiRS-DCF, QSCAN and KnockoffScreen.
- The R script ***sh_file_generator*** generates .sh file for performing GWAS on C50 Malignant neoplasms of breast using UK Biobank data using regenie.
- The .sh files ***fit.sh*** and ***test.sh*** perform GWAS on C50 Malignant neoplasms of breast using UK Biobank data using regenie.
