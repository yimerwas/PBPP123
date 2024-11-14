# Overview

This repository contains the data and scripts used for the manuscript titled "Formulating three classes of partial borrowing power priors to leverage historical data in process validation". The materials provided here aim to support the reproducibility of our findings and facilitate further research.

## Table of Contents

- [Introduction](#introduction)
- [Files](#files-overview)
- [Data](#data-description)
- [Scripts](#scripts)
- [Contact](#contact)

## Introduction

In the field of pharmaceutical manufacturing, ensuring product quality and safety is paramount, particularly during the process validation stage. The power prior approach, known for its effectiveness in decision-making where historical data is available, holds significant promise for enhancing validation methodologies in this context. However, despite its potential, the application of power priors in process validation has been limited, particularly in the area of partial borrowing of historical information from early development batches.

This manuscript seeks to address this notable gap by exploring the formulation and evaluation of three distinct classes of partial borrowing power priors (PBPP): PBPP1, a partial borrowing power prior with a fixed discounting parameter; PBPP2, an unnormalized partial borrowing power prior with a random discounting parameter; and PBPP3, a normalized partial borrowing power prior with a random discounting parameter. Each of these methodologies utilizes a Bayesian linear mixed model that effectively accounts for both intra-batch and inter-batch variability in drug product potency, thus providing a robust framework for integrating historical data into the process validation framework.

To illustrate the practical applicability of these approaches, we present a simulation study based on real drug product data. This simulation aims to elucidate how our proposed methods can leverage historical data to determine the optimal number of samples per batch required to accurately characterize intra-batch variability for process performance qualification (PPQ). In particular, we underscore the capabilities of PBPP3 in quantifying the discrepancies between historical and current data through its discounting parameter, offering valuable insights for decision-making processes within pharmaceutical manufacturing.

## Files

This repository contains the following files:

- `Data/`: Directory containing the datasets used in the analysis.
- `Scripts/`: Directory containing the scripts for both modeling the motivating data sets and performing simulation.
- `README.md`: This README file providing an overview of the project.

## Data
This repository includes two key datasets located in the `Data/` directory that are critical for the main analysis conducted in this manuscript:

**CurrentData.rds**

 - **Description**: This file contains the current product potency data simulated from the late-stage development process, it includes assay measurements from 10 batches (B1-B10), each containing 5 observations.
 - **Format**: This file is in R's serialized format (RDS), which can be easily loaded into R for further analysis.

**HistoricalData.rds**

 - **Description**: This file contains historical potency data simulated from earlier development batches of the same product, it includes sassay measurements from 20 batches (B1-B20), each with 20 observations.
 - **Format**: This file is also in R's serialized format (RDS).

## Scripts

This repository contains several scripts organized into different directories to facilitate the main analysis and simulations. Below is a summary of the contents:

### Main
- **Main.html**: The output file containing the results and code from the main analysis.
- **Main.Rmd**: The R Markdown file for the main analysis, which includes modeling with the motivating datasets and performs simulations.

### RCodes
The `RCodes` directory includes the following R scripts:
- **Approximate_Ca0.R**: A script for approximating the normalizing factor.
- **CreateStanData_Univ.R**: A script used for creating data in the format required by Stan.
- **OCurve4Assay_Univ.R**: A script for generating operating characteristic (OC) curves and determining the minimum sample size.

### StanCodes
The `StanCodes` directory contains Stan model files for different analyses:
- **BLMM_Current.stan**: Bayesian Linear Mixed Model (BLMM) for current data only.
- **BLMM_Historical.stan**: BLMM for historical data only.
- **BLMM_PartialBorrowing_Fixed_a0.stan**: BLMM modeling using PBPP1 methodology.
- **BLMM_PartialBorrowing_Normalized_Random_a0_2.stan**: BLMM modeling employing PBPP3 methodology.
- **BLMM_PartialBorrowing_Normalized_Random_a0_Prior_2.stan**: BLMM modeling using historical data with a discounting parameter \(a_0\).
- **BLMM_PartialBorrowing_Unnormalized_Random_a0.stan**: BLMM modeling utilizing PBPP2.

### Simulation
The `Simulation` folder comprises:
- **CodeToRunSimulations.R**: A script that contains the code required to run the simulations.
- **results/**: A folder containing simulation results that are imported and plotted in `Main.Rmd`.
- **stan/**: A folder holding the Stan files utilized for the simulations.

## Contact

If you have any questions, comments, or suggestions regarding this project, please feel free to reach out to me:

- **Name**: Yimer Wasihun Kifle 
- **Email**: yimerwas@gmail.com
- **GitHub**: [GitHub Profile](https://github.com/yimerwas)
- **LinkedIn**: [LinkedIn Profile](https://www.linkedin.com/in/yimerwas)


