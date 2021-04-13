# thesis_chapter_2

# Characterising Second-Generation Immunological Markers of Vaccine Induced Protection 
## Modelling the Roles of Antibody Titre and Avidity in Protection from Malaria Infection Following RTS,S/AS01 Vaccination 

## Overview 

Characterising the immune response after vaccination enables us to try and identify correlates or surrogates of vaccine induced protection. These correlates can facilitate downstream vaccine selection, help to identify ways to increase the immunogenicity of vaccines and aid in licensure decision making. Vaccine induced correlates of protection are often identified following large scale trials or routine implementation whose numbers can power statistical association tests. However during early clinical trials of malaria vaccine candidates human challenge studies are used and the sample sizes are relatively small which can often make statistical tests difficult. Despite the small sample sizes these trials offer an opportunity to characterise immune responses in malaria naive individuals with timed vaccine and infection exposure. In this work we use a Bayesian method to fit a biologically motivated mathematical model of individual level malaria infection to characterise the relationships between Circumsporozoite repeat region specific antibody titre and antibody avidity and protection from infection following vaccination in order to test the utility of both these measures in being able to predict vaccine efficacy.

This repository contains code used to generate the results and figures for the work in my thesis chapter 2:

Many thanks to Michael White, Institut Pasteur, who's original work (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0061395) formed the basis of this study and for his help in designing the MCMC algorithm that was adapted and used here. 

## Repo Content 
The R code presented is grouped according to the results subsections of the paper. 
1. [Data and packages](data)
2. [Immunogenicity Data](R1_Immunogenicity_Data)
3. [Model Fitting](R2_Model_Fitting)
4. [Model Predicted Vaccine Efficacy as a Function of the Immune Response](R3_Efficacy_Function_IR)

## System Requirements  
Performing this analysis requires only a standard computer with enough RAM to support running the code through R. This code has only been tested on a Windows system with the specifications as below:

* Processor: Intel(R) Core(TM) i7-7560U CPU @ 2.40GHz 
* Installed RAM: 16.0GB
* System Type: 64-bit operating system x4 based processor
* Windows specifications: Windows 10 Pro, Version 1803

### Software Requirements  
Users should have R version 3.4.0 or higher, and several packages set up from CRAN.

R can be downloaded for free, the latest version R-3.6.0 is available to download here: https://cran.r-project.org/bin/windows/base/ 

Users should install the following R packages prior to running any analysis, packages can be installed from the R terminal as follows: 

```install.packages(c('ggplot2', 'wesanderson', 'ggpubr', 'compiler', 'MASS', 'binom', 'survival', 'survminer', 'ResourceSelection', 'RColorBrewer', 'fields' , 'NumDeriv', 'patchwork', 'LaCroixColoR' )) ```

Install time of R should take no longer than a few minutes along with package installation. 

## Installation Guide and Instructions for Use 
To perform this analysis, first load the required packages into your R console using the code or RStudio Packages search tool.  

To replicate and reproduce the analyses presented in this paper, do the following:

Download MAL71_data.csv from the [Data](Data) folder of this repository, and use the [R code in the folder](Data/data_processing) to load into your environment. Ensure you have ```setwd("....")``` to the same file location as where you have downloaded the data and as such a suitable location to store processed results. 

Then follow through the remaining sections in the repo folders to reproduce the results as they appear in order in the Results section of the paper. The output of running the code will be a number of MCMC objects, as well as a series of plots representing the output from both running the models and using the model parameter distributions to produce predictive results. These plots form the basis of the plots presented in the paper. 
