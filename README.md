# A Two-Stage Modeling Approach to Risk Prediction with High-Dimensional Two-Phase Data
This repository provides the code required to reproduce the simulation results presented in the Two-Stage High-Dimension paper.
## Files:
1. **Format_Data.R**<br>
Contains the main function and helper functions used to generate and format the simulated data.

2. **Model_Estimation.R**<br>
Includes functions for estimating log-odds ratios and their asymptotic variances under the proposed framework (with stage I models based on adaptive lasso or Super Learner), as well as functions for estimating log-odds ratios using the two-stage benchmark model and complete-case (CC) analysis.

4. **Model_Evaluation.R**<br>
Provides functions for calculating AUC estimates and their asymptotic variances under the proposed framework, along with functions for computing AUC estimates for the two-stage benchmark model and CC analysis.

5. **Simulation_Run**<br>
This folder contains all simulation scenarios considered in the study. The scenarios vary according to the following design factors:<br>
   (1) **Missingness mechanism:** Linear, Single-index (Single), Simple Random Sampling (SRS)<br>
   (2) **Availability of Z:** 0.05 or 0.2<br>
   (3) **Correlation used in data generation:** 0, 0.3 or 0.5<br>
Each file corresponds to a specific simulation scenario. The scripts demonstrate how to generate simulation datasets, estimate log-odds ratios and AUC based on various models, and obtain the asymptotic variances of the log-odds ratios and AUC from the proposed method (with stage I models based on adaptive lasso).


**Note**: The simulation may take a long time to complete.
