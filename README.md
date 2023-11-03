## Bayesian Models for CE Analysis of a Prognostic Model for Sudden Cardiac Death (SCD) in Hypertrophic Cardiomyopathy (HCM)

### Background

Hypertrophic cardiomyopathy(HCM) is a leading cause of sudden cardiac death (SCD) in young adults. Current risk algorithms provide only a crude estimate of risk and fail to account for the different effect size of individual risk factors.

### Aim

The aim of this study was to perform a cost-effectiveness analysis for the new SCD risk prediction model. See for more information:
[C O'Mahony (2014) A novel clinical risk prediction model for sudden cardiac death in hypertrophic cardiomyopathy (HCM risk-SCD). Eur Heart;35(30):2010-20.](https://pubmed.ncbi.nlm.nih.gov/24126876/)

### Overview

We develop a predictive model for assessing risk of SCD among patients with HCM.

-   Using individual level data generate posterior distribution samples of probabilities, costs and health values using WinBUGS
-   Simulate with a Markov model using the posterior samples as inputs
-   Compare outcomes for risk groups in terms of cost-effectiveness

The model structure consists of a decision tree component and a Markov model component. In effect this first model is represented by different cohorts of the data in the fitting process.

<div class="figure" style="text-align: center">
<img src="images/model_diagram.png" alt="Markov model diagram. Bold circles represent starting states and dashed circles represent sink states." width="234" />
<p class="caption">Markov model diagram. Bold circles represent starting states and dashed circles represent sink states.</p>
</div>

### Data

Individual HCM patient level data used in this project:

| Variable     | Description                                                                  |
---------------|-------------------------------------------------------------------------------
| set          | "Imputed dataset"                                                            |
| centre       | "Centre"                                                                     |
| id           | "ID"                                                                         |
| d            | "Censoring Indicator (1=SCD/ICD shock)"                                      |
| time         | "Time (years)"                                                               |
| age          | "Age (years)"                                                                |
| mwt          | "Maximal wall thickness (mm)"                                                |
| mwt2         | "mwt^2"                                                                      |
| la           | "Left atrial size (mm)"                                                      |
| maxlvotg     | "Maximal LVOT gradient (mmHg)"                                               |
| fhxscd       | "Family history of SCD"                                                      |
| nsvt         | "Non Sustained Ventricular Tachycardia"                                      |
| syncope      | "Syncope (blackout)"                                                         

## Workflow

The main files to perform the analysis are:

-   `prep_study_data.R` munges the raw data and save the input data in `data/`
-   `BUGS/script.R` runs the BUGS code in `BUGS/model.txt`
-   `main-ce-analysis.R` performs the cost-effectiveness analysis
-   `pop_counts_plot.R` creates output plots

