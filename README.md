# Insight
This repository accompanies the following paper, currently under review:
 [Mihali, A, Broeker, M, Ragalmuto, F, Horga, G. Introspective inference counteracts perceptual distortion, 2022, bioRxiv](https://www.biorxiv.org/content/10.1101/2021.11.13.468497v4)

## Software requirements

This package has been tested on macOS Monterey (12.5.1) with the following:

-  Matlab 9.6 (2019a, MathWorks, Massachusetts, USA) 

- Alongside our Matlab scripts, please also download the [Bayesian adaptive direct search (BADS)  optimization algorithm](https://github.com/lacerbi/bads) by Luigi Acerbi.


## Demo and instructions for use: Data

We have datasets from Exp. 1 and Exp. 2 as `alldata_E1.mat` and respectively `alldata_E2.mat`, which are Nsubj * Ncond. Nsubj is 22 for each, but for demo we will show the data of 1 participant from Exp. 2, stored in `alldata_demo_E2.mat`. There are 4 conditions, combined such that they are in the following order: 

- No-Adapt-See
- No-Adapt-Believe
- Adapt-See
- Adapt-Believe


All analyses were performed on the `alldata` structs. For each subject and condition, there are 5 fields, each with variables of length Ntrials (121):


- stims: uniformly distributed on [-0.3, 0.3]
- resp: 0 (left/CCW) or 1 (right/CW)
- conf: 0 for high confidence, 1 for low confidence 
- resp_times (sec)
- conf_times (sec)


## Scripts

- To fit the psychometric curves with the lapse parameter shared across the 4 conditions, we used the function `psych_curve_fitting_offline_m2.m`, which uses grid search for optimization and calls  `loglike.m`, which calls `function_psi.m`.

- We fit the 5 Bayesian model variants via `fitting_exp_alll.m`. This script is ready to run in the demo version. It calls  `Fitting_model.m`, which finds the parameters that maximize the probability of the data (likelihood) given each model via `Loglike_bm_alll.m`. 

- `analysis_all.m` outputs the main data figures from our manuscript. To overlay the model predictions on top of the data, it calls `Predict_bm_alll.m`. This script will only be able to run once we make the full datasets available. 
