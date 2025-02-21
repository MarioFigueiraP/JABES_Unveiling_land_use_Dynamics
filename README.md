# Unveiling Land Use Dynamics: Insights from a Hierarchical Bayesian Spatio-Temporal Modelling of Compositional Data
Code for reproducing the main results presented in the article 'Unveiling Land Use Dynamics: Insights from a Hierarchical Bayesian Spatio-Temporal Modelling of Compositional Data', published in JABES.

## 1- Dealing with 0's in land use data:

The scripts associated with the section on examples related to the presence of zeros in compositional data correspond to [**runme_Zeros_and_Ones_LerouxExample_Beta.R**](./runme_Zeros_and_Ones_LerouxExample_Beta.R) (for univariate data following a Beta distribution) or [**runme_Zeros_and_Ones_LerouxExample_CLR.R**](./runme_Zeros_and_Ones_LerouxExample_CLR.R) (for multivariate compositional data evaluated using a log-ratio transformation).

## 2- Downscaling models

The scripts associated with the downscaling models are [**runme_Beta_Downscaling_Example.R**](./runme_Beta_Downscaling_Example.R), for data following a Beta distribution, or [**runme_Beta_Downscaling_Example.R**](./runme_Beta_Downscaling_Example.R) for the multivariate case using the CLR.

## 3- Big data

The script for Big Data is [**runme_BigData_Example.R**](./runme_BigData_Example.R). This script simulates a large database and demonstrates how to perform an inferential process using a **sequential Consensus** approach.

## 4- Extra

Finally, script [**runme_Leroux_Example_ARL.R**](./runme_Leroux_Example_ARL.R) is an additional code example for simulating compositional data in small areas, with a spatial effect distributed according to a Leroux distribution. The compositions are evaluated using additive log-ratios (ALR), such that in the multivariate normal associated with the ALRs, the correlation matrix between the different compositions is analyzed.
