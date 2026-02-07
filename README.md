# A Two-Stage Modeling Approach to Risk Prediction with High-Dimensional Two-Phase Data
This repository provides the code required to reproduce the simulation results presented in the Two-Stage High-Dimension paper.
## Files:
1. **Format_Data.R**<br>
Contains the main function and helper functions used to generate and format the simulated data.

2. **Model_Estimation.R**<br>
Includes functions for estimating the log-odds ratios and their asymptotic variances using the proposed method, as well as functions for estimating log-odds ratios using the two-stage model (benchmark).

3. **Model_Evaluation.R**<br>
Provides functions for calculating AUC estimates and their asymptotic variances using the proposed method, along with functions for computing AUC estimates for the two-stage model (benchmark).

4. **Simulation_Run**<br>
Inculdes all scenarios in the simulation studies:<br>
(1) Missingness mechanisms: linear, single-index(Single), simple random sampling(SRS)<br>
(2) Availability of Z (P(R=1)): 0.05 or 0.2<br>
(3) correlation when generating data: 0, 0.3 or 0.5<br>
Each file corresponds to one scenario, which demonstrates how to generate a simulation dataset, estimate log-odds ratios and AUC using both the proposed method and the two-stage benchmark model, and obtain the asymptotic variances of the log-odds ratios and AUC from the proposed method.

**Note**: The simulation may take a long time to complete.
