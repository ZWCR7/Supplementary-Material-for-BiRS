# Supplementary-Material-for-BiRS
Supplementary codes and tables for BiRS

## Data and software availability
- Previous data published in UK Biobank were used for this work and this research has been conducted using UK Biobank Resource under project 79237. One can refer to 
https://biobank.ctsu.ox.ac.uk/crystal/index.cgi for accessing and enabling data download or detailed description of UKB data.
- All WGS studies are conducted in the UK Biobank Research Analysis Platform (RAP), see https://ukbiobank.dnanexus.com/landing for details.
- PLINK2.0 was used to pre-treat UKB data, see https://www.cog-genomics.org/plink/2.0/ for detailed manuls of using PLINK2.0.
- cosi2 was used to generated simulation data, see https://github.com/broadinstitute/cosi2 for detaild manuls of using cosi2.

## GWAS and WGS results
- *GWAS-BiRS.csv*: GWAS results of BiRS
- *GWAS-Scan.csv*: GWAS results of Q-SCAN
- *GWAS-KS.csv*: GWAS results of KnockoffScreen
- *GWAS-4S.csv*: GWAS results of 4S
- *GWAS-REGENIE.csv*: GWAS results of regenie
- *WGS-BiRS-Region*: Region-wise WGS results of BiRS

## Code for BiRS
### Overview

1. Package *BiRS* contains the main functions of BiRS algorithm and DCF test. One can install the package *BiRS* through the file ***BiRS_0.1.0.tar.gz***.  


2. Directory *Simulation_BiRS* contains the main functions for the simulations. Specifically, 
- The R script ***KSDet*** contains the code for applying KnockoffScreen method to normal setting in simulation.
- The R script ***4S_Algorithm*** contains the code for applying 4S method to simulation.
- The R script ***LRS_Detect*** contains the code for applying LRS method to simulation.

- The R script ***M-dependence/Normal_ESBA*** contains the code for generating data and conducting simulation for ''M-dependence with equal covariance and non-decay signal (MES)'' setting using BiRS-DCF, BiRS-CL, Q-SCAN, KnockoffScreen, 4S and LRS.
- The R script ***M-dependence/Normal_NSBA*** contains the code for generating data and conducting simulation for ''M-dependence with unequal covariance and non-decay signal (MNS)'' setting using BiRS-DCF, BiRS-CL, Q-SCAN, KnockoffScreen, 4S and LRS.
- The R script ***M-dependence/Normal_ESBA_Decay*** conntains the code for generating data and conducting simulation for ''M-dependence with equal covariance and decay signal (WES)'' setting using BiRS-DCF, BiRS-CL, Q-SCAN, KnockoffScreen, 4S and LRS.
- The R script ***M-dependence/Normal_NSBA_Decay*** contains the code for generating data and conducting simulation for ''M-dependence with unequal covariance and decay signal (WNS)'' setting using BiRS-DCF, BiRS-CL, Q-SCAN, KnockoffScreen, 4S and LRS.
- The R scripts ***M-dependence/Describe*** and ***M-dependence/Describe_Decay*** contain the code for generating description information of the detection results under M-dependence covariance structure.  
- The R scripts ***M-dependence/Organize*** and ***M-dependence/Organize_Decay*** contain the code for generating figures and tables of the detection results under M-dependence covariance structure.

- The R script ***Weak-dependence/Normal_ESBA*** contains the code for generating data and conducting simulation for ''Weak-dependence with equal covariance and non-decay signal'' setting using BiRS-DCF, BiRS-CL, Q-SCAN, KnockoffScreen, 4S and LRS.
- The R script ***Weak-dependence/Normal_NSBA*** contains the code for generating data and conducting simulation for ''Weak-dependence with unequal covariance and non-decay signal'' setting using BiRS-DCF, BiRS-CL, Q-SCAN, KnockoffScreen, 4S and LRS.
- The R script ***Weak-dependence/Normal_ESBA_Decay*** conntains the code for generating data and conducting simulation for ''Weak-dependence with equal covariance and decay signal'' setting using BiRS-DCF, BiRS-CL, Q-SCAN, KnockoffScreen, 4S and LRS.
- The R script ***Weak-dependence/Normal_NSBA_Decay*** contains the code for generating data and conducting simulation for ''Weak-dependence with unequal covariance and decay signal'' setting using BiRS-DCF, BiRS-CL, Q-SCAN, KnockoffScreen, 4S and LRS.
- The R scripts ***Weak-dependence/Describe*** and ***Weak-dependence/Describe_Decay*** contain the code for generating description information of the detection results under Weak-dependence covariance structure. 
- The R scripts ***Weak-dependence/Organize*** and ***Weak-dependence/Organize_Decay*** contain the code for generating figures and tables of the detection results under Weak-dependence covariance structure.

- The scripts ***Genetic_Setting/parama_ukb*** and ***Genetic_Setting/recParams_ukb*** contains the parameters for generating the haplotypes by cosi-2.
- The R scripts ***Genetic_Setting/BinaryGenerator*** contains the code for generating the binary files (.bim, .bed, .fam) of simulated haplotypes from cosi-2.
- The R scripts ***Genetic_Setting/get_corMat*** and ***Genetic_Setting/get_PDcorMat*** contain the code for generating the LD matrix of simulated sequence data.
- The R script ***Genetic_Setting/KSDetGenetic*** contains the code for applying KnockoffScreen method to genetic setting in simulation.
- The R script ***Genetic_Setting/Block_BiRS*** contains the code for applying BiRS method to genetic setting in simulation.
- The R script ***Genetic_Setting/Genetic_Setting*** contains the code for generating data and conducting simulation for genetic setting using BiRS-DCF, BiRS-CL, Q-SCAN and KnockoffScreen.
- The R scripts ***Genetic_Setting/Describe_Genetic*** contains the code for generating description information of the detection results in genetic setting. 
- The R scripts ***Genetic_Setting/Organize_Genetic*** contains the code for generating figures and tables of the detection results in genetic setting.
- The R scripts ***Genetic_Setting/Speed-Test*** and ***Genetic_Setting/Summary_Speed*** contains the code for calculating computational times of these methods under genetic setting.

- The R script ***Size-Test/Size-Test*** contains the contains the code for generating data and conducting simulation for size validation of BiRS-DCF, BiRS-CL, Q-SCAN, 4S and LRS in normal settings.
- The R script ***Size-Test/Genetic_Size*** contains the contains the code for generating data and conducting simulation for size validation of BiRS-DCF, BiRS-CL, Q-SCAN, 4S and SCANG-STAAR in genetic setting.
- The R scripts ***Size-Test/Describe_Size*** contains the code for generating tables of the size validation in normal settings.
- The R scripts ***Size-Test/Describe_Genetic_Size*** contains the code for generating tables of the size validation in genetic setting.


3. Directory *Application_Code* contains the main functions for quality control and performing GWAS and WGS on UK Biobank data. Specifically, 
- The R scripts ***Pretreatment/CTtraits*** and ***Pretreatment/covariate*** generate samples and covariates we needed from the main dataset of UK Biobank.
- The R scripts ***Pretreatment/SampleQC*** and ***Pretreatment/VariantQC*** perform the sample and variant quality control after generating .vcf file for each chromosome using PLINK 2.0 (the quality control on HWE and MAF performed at the same time as the .vcf file be generated)
- The R script ***Pretreatment/SampleGLM*** generates data for Q-SCAN to perform GWAS.
- The R script ***Pretreatment/GenerateFinalVariant*** generates the information for the selected variants for performing GWAS.
- The R script ***Pretreatment/SizeSampleGenerator*** generates the samples for performing permutation test in each chromosome.
- The R scripts ***BiRS/BiRSGenomeSize***, ***QSCAN/QScanGenomeSize*** and ***BiRS/BiCLGenomeSize*** perform permutation tests of BiRS-DCF, BiRS-CL and QSCAN.
- The R scripts ***BiRS/BiRSGenome***, ***QSCAN/QScanGenome***, ***KnockoffScreen/KSGWAS*** and ***4S_GWAS*** perform GWAS on C50 Malignant neoplasms of breast using UK Biobank data by BiRS-DCF, QSCAN, KnockoffScreen and 4S.
- The R script ***sh_file_generator*** generates .sh file for performing GWAS on C50 Malignant neoplasms of breast using UK Biobank data using regenie.
- The .sh files ***fit.sh*** and ***test.sh*** perform GWAS on C50 Malignant neoplasms of breast using UK Biobank data using regenie.
- The R script ***plot_gene*** plots identified signal regions of the four methods.
  
- The R script ***WGS/Generate_Block_Index*** contains the code for generating blocks for detection.
- The R script ***WGS/BiRS-WGS-Chr1*** performs whole genome sequencing association studies in chromosome 1 in rap platform (for other chromosomes, change the parameter CHR).
- The R script ***WGS/BlockBiRS_WGS*** contains the code for applying BiRS method to WGS in rap platform.
- The R script ***WGS/Summary_WGS*** contains the code for summarizing the WGS results.


4. Directory *Supplementary_Code* contains the additional simulations in Supplementary Material. Specifically, 
- The R script ***BiRS-Cpts/SummaryBiRS*** and ***BiRS-Cpts/WildBiRS*** contain the code for applying WBB.
- The R script ***BiRS-Cpts/MES-BiRS-Cpts***,  ***BiRS-Cpts/MNS-BiRS-Cpts***,  ***BiRS-Cpts/WES-BiRS-Cpts*** and  ***BiRS-Cpts/MES-BiRS-Cpts*** contain the codes for generating data and conducting simulations for MES, MNS, WES and WNS setting.
- The R scripts ***BiRS-Cpts/Summary*** contains the code for generating description information of the detection results under all settings.  
- The R scripts ***BiRS-Cpts/plot*** contains the code for generating figures of the detection results.

- The R script ***Sensitivity/MES-Sensitivity***,  ***Sensitivity/MNS-Sensitivity***,  ***Sensitivity/WES-Sensitivity*** and  ***Sensitivity/WNS-Sensitivity*** contain the codes for generating data and conducting simulations for MES, MNS, WES and WNS setting.
- The R scripts ***Sensitivity/Summary*** contains the code for generating description information and figures of the detection results under all settings.
  

### Workflows
Overall, please first install the *BiRS* package through *BiRS_0.1.0.tar.gz*. 
#### Simulation for comparing detection performance of different methods under M-dependence
1. Please create folders *Simulation_BiRS/M-dependence/Normal_ESBA/* and *Simulation_BiRS/M-dependence/Normal_NSBA/* for saving the simulation results in normal setting with M-dependence covariance structure generated by ***Simulation_BiRS/M-dependence/Normal_ESBA.R*** and ***Simulation_BiRS/M-dependence/Normal_NSBA.R***, respectively, (the code results are saved as .RData files). Then run ***Simulation_BiRS/M-dependence/Normal_ESBA.R*** and ***Simulation_BiRS/M-dependence/Normal_NSBA.R***.
2. Run ***Simulation_BiRS/M-dependence/Describe.R*** for generating description information of the detection results under M-dependence covariance structure.
3. Create folder *Simulation_BiRS/M-dependence/Figures/* for saving the figures of the detection results under M-dependence covariance structure then run ***Simulation_BiRS/M-dependence/Organize.R*** to generate them.

#### Simulation for comparing detection performance of different methods under Weak-dependence
1. Please create folders *Simulation_BiRS/Weak-dependence/Normal_ESBA/* and *Simulation_BiRS/Weak-dependence/Normal_NSBA/* for saving the simulation results in normal setting with Weak-dependence covariance structure generated by ***Simulation_BiRS/Weak-dependence/Normal_ESBA.R*** and ***Simulation_BiRS/Weak-dependence/Normal_NSBA.R***, respectively, (the code results are saved as .RData files). Then run ***Simulation_BiRS/Weak-dependence/Normal_ESBA.R*** and ***Simulation_BiRS/Weak-dependence/Normal_NSBA.R***.
2. Run ***Simulation_BiRS/Weak-dependence/Describe.R*** for generating description information of the detection results under Weak-dependence covariance structure.
3. Create folder *Simulation_BiRS/Weak-dependence/Figures/* for saving the figures of the detection results under Weak-dependence covariance structure then run ***Simulation_BiRS/Weak-dependence/Organize.R*** to generate them.

#### Simulation for comparing detection performance of different methods under genetic settings
1. Similarily, please make sure you have created correct folders to store simulation results.
2. Using cosi2 and the parameter files ***Simulation_BiRS/Genetic_Setting/parama_ukb*** and ***Simulation_BiRS/Genetic_Setting/recParams_ukb*** to generate haplotypes.
3. Create folder *Simulation_BiRS/Genetic_Setting/ukb_simu_txt/* and use the command line to separate the generated haplotype file into 1000 part and run ***Simulation_BiRS/Genetic_Setting/BinaryGenerator*** to get the binary files for sequence data (.bim, .bed, .fam), please create folder *Simulation_BiRS/Genetic_Setting/ukb_simu_binary/*.
4. Using plink2 to merge the binary files in *Simulation_BiRS/Genetic_Setting/ukb_simu_binary/* into ***Simulation_BiRS/Genetic_Setting/genotype.bed (.fam, .bim)***.
5. Run ***Simulation_BiRS/Genetic_Setting/get_corMat*** and ***Simulation_BiRS/Genetic_Setting/get_PDcorMat*** to get the LD matrix, MAF and part numbers (part numbers.RData) for generating data, please create folders *Simulation_BiRS/Genetic_Setting/corMatrix/*, *Simulation_BiRS/Genetic_Setting/corMatrixPD/* and *Simulation_BiRS/Genetic_Setting/MAFS/*.
6. Run ***Simulation_BiRS/Genetic_Setting/Genetic_Setting.R*** to generate the detection results under genetic setting.
7. Run ***Simulation_BiRS/Genetic_Setting/Describe_Genetic.R*** for generating description information of the detection results under genetic setting.
8. Create folder *Simulation_BiRS/Genetic_Setting/Figures/* for saving the figures of the detection results under genetic setting then run ***Simulation_BiRS/Genetic_Setting/Organize_Genetic.R*** to generate them.

#### Simulation for comparing the empirical FWERs of different methods for different settings.
1. Similarily, please make sure you have created correct folders to store simulation results.
2. Run ***Simulation_BiRS/Size-Test/Size-Test.R*** and ***Simulation_BiRS/Size-Test/Genetic_Size.R*** to generate the empirical FWERs of different methods for different settings (saving as .RData file).
3. Run ***Simulation_BiRS/Size-Test/Describe_Size.R*** and ***Simulation_BiRS/Size-Test/Describe_Genetic_Size.R*** for generating tables of the empirical FWERs of different methods for different settings, which are the Table 1 and 2 in the main paper.

#### Simulation for comparing the computational time for different methods.
1. Similarily, please make sure you have created correct folders to store simulation results.
2. Run ***Simulation_BiRS/Genetic_Setting/Speed-test.R*** and ***Simulation_BiRS/Genetic_Setting/Summary_Speed.R*** to generate the computational time for different methods.


#### Data application: pre-treatment the ukbiobank data
1. Please make sure that you can access the ukbiobank main dataset contains the fields list in ***Application_Code/Pretreatment/CTtraits.R***, put them in the suitable directory (one can create the directories in my codes or use its own directory then change the directories in my codes to its own). Then run ***Application_Code/Pretreatment/CTtraits.R*** to generate the .RData file containing phenotypes and covariates.
2. Run ***Application_Code/Pretreatment/covariate.R*** to generate covariates used for Q-SCAN, REGENIE and KnockoffScreen.
3. Please make sure that you can access the ukbiobank genomic field 22418 then generating .vcf file for each chromosome using PLINK 2.0 (the quality control on HWE and MAF performed at the same time as the .vcf file be generated), saving in your own directory or the directory I give in the code.
3. Run ***Application_Code/Pretreatment/SampleQC.R*** and ***Application_Code/Pretreatment/VariantQC.R*** perform further sample and variant quality control on the .vcf file generated prerviously.
4. Run ***Application_Code/Pretreatment/GenerateFinalVariant.R*** to generate the information for the selected variants after sample and variant QC for performing GWAS.
4. Run ***Application_Code/Pretreatment/SampleGLM.R*** to generate data for Q-SCAN to perform GWAS.
5. Run ***Application_Code/Pretreatment/SizeSampleGenerator.R*** to generate the samples for performing permutation test in each chromosome.

#### Data application: performing permutation tests
1. Run ***Application_Code/BiRS/BiCLGenomeSize.R*** and ***Application_Code/BiRS/BiRSGenomeSize.R*** to perform permutation test for calculating empirical FWERs in each chromosome for BiRS-DCF and BiRS-CL.
2. Run ***Application_Code/Q-SCAN/QScanGenomeSize.R*** to perform permutation test for calculating empirical FWER in each chromosome for Q-SCAN.

#### Data application: perform GWASs using different methods
1. Run ***Application_Code/REGENIE/sh_file_generator.R*** to generate .sh files for perform indvidual GWAS in each chromosome.
2. Run ***Application_Code/REGENIE/fit.sh*** to fit the regression model for generating test statistics.
3. Run ***Application_Code/REGENIE/test.sh*** to calculate the p-values of each snp then select the significant snps contained in *GWAS-REGENIE.csv*.
4. Run ***Application_Code/KnockoffScreen/KSGWAS.R*** to perform region-wise GWAS on the genomic data using KnockoffScreen and generate the GWAS results contained in *GWAS-KS.csv*.
5. Run ***Application_Code/QSCAN/QScanGenome.R*** to perform region-wise GWAS on the genomic data using Q-SCAN and generate the GWAS results contained in *GWAS-Scan.csv*.
6. Run ***Application_Code/BiRS/BiRSGenome.R*** to perform region-wise GWAS on the genomic data using BiRS-DCF and generate the GWAS results contained in *GWAS-BiRS.csv*.
7. Run ***Application_Code/4S/4S_GWAS.R*** to perform region-wise GWAS on the genomic data using 4S and generate the GWAS results contained in *GWAS-4S.csv*.

#### Data application: visualize the detected signal regions of GWAS of different methods
1. Please make sure that the .csv files *C50Breast_BiRS.csv*, *C50Breast_Scan*, *C50Breast_KSGWAS*, *C50Breast_Regenie*, *C50Breast_4S*, *C50Breast_chr11_gene.csv*, *C50Breast_chr12_gene.csv*, *C50Breast_chr19_gene.csv* and Rscript ***plot_gene*** in the same directory.
2. Run ***plot_gene*** to visualize the detected signal regions and related genes.

#### Data application: perform WGS using BiRS-DCF
1. Make sure you can access the 200k WGS data in ukbiobank and conduct analysis in RAP.
2. Copy the .RData files *Binary_Sample_field40006_breast_cancer_covariate_phenotype.RData* and *KSID.RData* generated during GWAS to RAP.
3. Run ***Application_Code/WGS/Get_WGS_Samples.R*** to get the samples after quality control.
4. Run ***Application_Code/WGS/Generate_Block_Index.R*** to get the blocks for detection.
5. Run ***Application_Code/WGS/BiRS_WGS_Chr1.R*** to get the WGS results in chromosome 1 (for chromosomes 2--22, change the parameter CHR).
6. Run ***Application_Code/WGS/Summary_WGS.R*** to get the final results.


#### Supplement Simulation: Simulation for comparing detection performance of WBB for all normal settings in Supplementary Material.
1. Please create correct folders to store the simulation results.
2. Run ***Supplementary_Code/BiRS-Cpts/MES-BiRS-Cpts.R***, ***Supplementary_Code/BiRS-Cpts/MNS-BiRS-Cpts.R***, ***Supplementary_Code/BiRS-Cpts/WES-BiRS-Cpts.R*** and ***Supplementary_Code/BiRS-Cpts/WNS-BiRS-Cpts.R*** to get the simulation results of WBB.
3. Run ***Supplementary_Code/BiRS-Cpts/Summary.R*** and ***Supplementary_Code/BiRS-Cpts/plot.R*** to summarize and plot the simulation results.
 
#### Supplement Simulation: Simulation for the sensitivity of truncation parameter s Supplementary Material.
1. Please create correct folders to store the simulation results.
2. Run ***Supplementary_Code/Sensitivity/MES-Sensitive.R***, ***Supplementary_Code/Sensitivity/MNS-Sensitive.R***, ***Supplementary_Code/Sensitivity/WES-Sensitive.R*** and ***Supplementary_Code/Sensitivity/WNS-Sensitive.R*** to get the simulation results of sensitivity.
3. Run ***Supplementary_Code/Sensitivity/Summary.R*** to summarize and plot the simulation results.
   
