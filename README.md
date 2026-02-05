# A Two-Stage Modeling Approach to Risk Prediction with High-Dimensional Two-Phase Data
This repository provides the code required to reproduce the simulation results presented in the Two-Stage High-Dimension paper.
## Files:
1. **Format_Data.R**<br>
Contains the main function and helper functions used to generate and format the simulated data.

2. **Model_Estimation.R**<br>
Includes functions for estimating the log-odds ratios and their asymptotic variances using the proposed method, as well as functions for estimating log-odds ratios using the two-stage model (benchmark).

3. **Model_Evaluation.R**<br>
Provides functions for calculating AUC estimates and their asymptotic variances using the proposed method, along with functions for computing AUC estimates for the two-stage model (benchmark).

4. **Example_Run.R**
Demonstrates how to generate a simulation dataset, estimate log-odds ratios and AUC using both the proposed method and the two-stage benchmark model, and obtain the asymptotic variances of the log-odds ratios and AUC from the proposed method.

**Note**: The simulation may take a long time to complete.
