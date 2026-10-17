Here are all codes for figures 1-2, B1-B2 and calculations in the manuscript named "Nonlinear-parasite-interactions-alter-the-scaling-up-of-coinfection-patterns".

Analyses were run in R version 4.5.2 and with TMB 1.9.17. Bayesian analyses are run with JAGS 4.3.2

Data files

Data/Coinfection_data_2017.csv: 2017 survey data reported in Hasik et al. (2024)

Data/Parasite_Data_2019.csv: 2019 survey data reported in Hasik et al. (2024)

SummarizedData_2017.csv: Environmental variables collected for 2017 survey data.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
SummarizedData_2019.csv Environmental variables collected for 2019 survey data.

Code files: All code files require the data files to run.

Analysis_Tables3_5_Fig2_TableB2_B3_TableC1_TableC3_C5.R: Analysis of Hasik et al. (2024) dataset using the full mixed-effect generalized linear model described in Jiao et al (2026) which uses the numeric responses for both mites and gregarines. The output from this analysis was used to generate Figure 2, Tables 3-5, Tables B2-B3, Tables C1, and Tables C3-C5.

PartialReduction_analysis_Table6.R:  Analysis of Hasik et al. (2024) dataset using the partially reduced data resolution by treating the response variable (gregarine abundance) as presence-absence. The output from this analysis was used to generate Table 6.

FullReduction_analysis_Table_7_TableC2.R: Analysis of Hasik et al. (2024) dataset using the fully reduced data resolution by treating both the response variable (gregarine abundance) and the biological predictor (mite abundance) as presence-absence. The output from this analysis was used to generate Table7 and Table C2.

cscore_Fig3.R: Simulation of data and application of the nonparametric analysis to the simulated data using the nonparametric C-score. Used to generate figure 3.

Full_Bayesiananalysis_TableB1.R: The Bayesian analysis of the full dataset reproducing the analysis in Analysis_Tables3_5_Fig2_TableB2_B3_TableC1_TableC3_C5.R using a Bayesian model. Used to generate Table B1.

Reduced_Bayesian_TableB4_B5.R: The Bayesian analysis of the reduced datasets treating both gregarines and/or mites as presence-absence data. The output of these analyses was used to generate Tables B4 and B5.

Literature cited:

Hasik, Adam Z., Jason T. Bried, Daniel I. Bolnick, and Adam M. Siepielski. “Is the Local Environment More Important than Within-Host Interactions in Determining Coinfection?” Journal of Animal Ecology 93, no. 10 (2024): 1541–55. https://doi.org/10.1111/1365-2656.14167.

Jiao, J. Haski, A., Siepielski, Ferguson, J. M. Detecting Nonlinear Parasite Interactions: An Integrated Approach in a Damselfly Coinfection System. (2026) Journal of Animal Ecology.


