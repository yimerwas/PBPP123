# Overview

This repository contains the data and scripts used for the manuscript titled "Formulating three classes of partial borrowing power priors to leverage historical data in process validation". The materials provided here aim to support the reproducibility of our findings and facilitate further research.

## Table of Contents

- [Introduction](#introduction)
- [Files Overview](#files-overview)
- [Data Description](#data-description)
- [Scripts](#scripts)
- [Installation Instructions](#installation-instructions)
- [Usage](#usage)
- [License](#license)
- [Contact](#contact)

## Introduction

In the field of pharmaceutical manufacturing, ensuring product quality and safety is paramount, particularly during the process validation stage. The power prior approach, known for its effectiveness in decision-making where historical data is available, holds significant promise for enhancing validation methodologies in this context. However, despite its potential, the application of power priors in process validation has been limited, particularly in the area of partial borrowing of historical information from early development batches.

This manuscript seeks to address this notable gap by exploring the formulation and evaluation of three distinct classes of partial borrowing power priors (PBPP): PBPP1, a partial borrowing power prior with a fixed discounting parameter; PBPP2, an unnormalized partial borrowing power prior with a random discounting parameter; and PBPP3, a normalized partial borrowing power prior with a random discounting parameter. Each of these methodologies utilizes a Bayesian linear mixed model that effectively accounts for both intra-batch and inter-batch variability in drug product potency, thus providing a robust framework for integrating historical data into the process validation framework.

To illustrate the practical applicability of these approaches, we present a simulation study based on real drug product data. This simulation aims to elucidate how our proposed methods can leverage historical data to determine the optimal number of samples per batch required to accurately characterize intra-batch variability for process performance qualification (PPQ). In particular, we underscore the capabilities of PBPP3 in quantifying the discrepancies between historical and current data through its discounting parameter, offering valuable insights for decision-making processes within pharmaceutical manufacturing.

## Files Overview

This repository contains the following files:

- `Data/`: Directory containing the datasets used in the analysis.
- `Scripts/`: Directory containing the scripts for data analysis and visualization.
- `README.md`: This README file providing an overview of the project.

## Data
This repository includes two key datasets located in the Data/ directory that are critical for the main analysis conducted in this manuscript:

**CurrentData.rds**

 - **Description**: This file contains the current product potency data simulated from the late-stage development process, it includes assay measurements from 10 batches (B1-B10), each containing 5 observations.
 - **Format**: This file is in R's serialized format (RDS), which can be easily loaded into R for further analysis.

**HistoricalData.rds**

 - **Description**: This file contains historical potency data simulated from earlier development batches of the same product, it includes sassay measurements from 20 batches (B1-B20), each with 20 observations.
 - **Format**: This file is also in R's serialized format (RDS).


## Scripts

Provide a summary of the scripts available in the `scripts/` directory. Include information on what each script does and how it relates to the data:

- `analysis_script.py`: Description of what this script does.
- `visualization_script.R`: Description of what this script does.
- ...

## Installation Instructions

Outline any prerequisites needed to run your scripts, such as programming languages or libraries. For example:

