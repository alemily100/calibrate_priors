# calibrate_priors
This repository contains the code used to simulate results for the paper **Calibration of dose-agnostic priors for Bayesian dose-finding trial designs with joint outcomes**. 

## Background 
This repository presents the code for the simulation studies undertaken within this paper.

## Description of R files
* **procrm_functions.R** - code for Section 4.1: functions required to generate PRO-CRM trials and identify optimal doses.

* **procrm_sim.R** - code for Section 4.1: code required to run PRO-CRM simulation studies. Functions required to run this code are defined in `procrm_functions.R`.

* **efftox_functions.R** - code for Section 4.2: functions required to generate Eff+Tox-CRM trials and identify optimal doses.
  
* **efftox_sim.R** - code for Section 4.2: code required to run Eff+Tox-CRM simulation studies. Functions required to run this code are defined in `efftox_functions.R`.
