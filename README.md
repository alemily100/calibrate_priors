# calibrate_priors
This repository contains the code used to simulate results for the paper **Calibration of dose-agnostic priors for Bayesian dose-finding trial designs with joint outcomes**. 

## Background 
This repository presents the code for the simulation studies undertaken within this paper.

## Description of R files
* **functions.R** - code for Section 4: functions required to generate trials and identify the OBD for 5,000 simulations running in parallel. 
  
* **sim.R** - code for Section 4: code required to run PRO-ADD simulation studies. Functions required to run this code are defined in `functions.R`.

* **shape_param_inc.csv** - code for Section 4: .csv file containing matrix of shape parameters to define PRO-nAE burden score simualation scenarios.

* **rate_param_inc.csv** - code for Section 4: .csv file containing rate parameter to define PRO-nAE burden score simualation scenarios.

* **fig3_data.csv** - 1,000 posterior samples of probability of response and PRO-nAE burden score for each dose used to create Figure 3 in the main manuscript. 
