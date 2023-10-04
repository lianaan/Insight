# Insight
This repository accompanies the following paper, currently under review:
 [Mihali, A, Broeker, M, Ragalmuto, F, Horga, G. Introspective inference counteracts perceptual distortion, 2023, bioRxiv]( https://www.biorxiv.org/content/10.1101/2021.11.13.468497v5). In addition to the co-authors and several colleagues listed in the acknowledgements, also thanks to [@astroahad](https://github.com/astroahad)!
 
 This repository contains demo data and analysis and modeling code. We do not provide the experimental code here for now, but mention that to implement the spiral stimuli we adapted the implementation of [Peter Scarfe's Contrast Modulated Spiral](https://peterscarfe.com/contrastModulatedSpiral.html). 


[Watch the video of our MAE illusion implementation here](https://drive.google.com/file/d/1PPtd57ln27324gM14BaQqXUFCyQGcVSq/view?usp=sharing). Fixate on the center white dot (here for 15 seconds) and report the direction of motion of the spiral which appears when the central dot changes color from white to yellow.

## Software requirements
 
This package has been tested on macOS Monterey (12.5.1) with the following:

-  Matlab 9.6 (2019a, MathWorks, Massachusetts, USA) 

- We provide the output results of the Bayesian model fitting procedure.  To run our model fitting scripts, please also download the [Bayesian adaptive direct search (BADS)  optimization algorithm](https://github.com/lacerbi/bads) by Luigi Acerbi.


## Demo and instructions for use: Data

We have datasets from Exp. 1 and Exp. 2 as `alldata_E1.mat` and respectively `alldata_E2.mat`, which are Nsubj * Ncond, with Nsubj being 22 for each. There are 4 conditions, combined such that they are in the following order: 

- 1: No-Adapt-See 
- 2: No-Adapt-Believe
- 3: Adapt-See
- 4: Adapt-Believe

 All the participants performed Adapt-See first, followed by Adapt-Believe (except participant 2 in Experiment 2 who accidentally performed Adapt-Believe before Adapt-See). In Experiment 1, participants then performend No-Adapt-See followed by No-Adapt-Believe; in Experiment 2 the order of these last two conditions was randomized across participants.  
 
All analyses were performed on the `alldata` structs. For each subject and condition, there are 5 fields, each with variables of length Ntrials (121):


- stims: uniformly distributed on [-0.3, 0.3]
- resp: 0 (left/CCW) or 1 (right/CW)
- conf: 0 for high confidence, 1 for low confidence 
- resp_times (sec)
- conf_times (sec)



## Scripts
- For individual model fits, full simulated datasets and fits, complement the scripts from here with the data in the corresponding folders:
[More on data and model fits](https://drive.google.com/drive/folders/1OW1x80jKBBn9jowLEeM6Y8xWxc_w7NO4?usp=drive_link)

- `analysis_all.m` outputs a large number of data figures from our manuscript: Figures 3, 4, 5 and 7 from main, Figures S2, S3 from Supplementary, and, if we change the model index from mi = 1 to mi = 4, we could also get Figure S6.  To overlay the model predictions on top of the data, it calls `Predict_bm_alll.m`. It loads the entire dataset, psychometric curve fits (`psych_curves_fitting_m2_201_E2.mat` for exp.2) and, for exp_i = 2, Bayesian model fits (`params_all_models6_E2_Nsubj_22.mat`). 

## Psychometric curve fitting
- To fit the psychometric curves with the lapse parameter shared across the 4 conditions, we used the function `psych_curve_fitting_offline_m2.m`, which uses grid search for optimization and calls  `loglike.m`, which calls `function_psi.m`.

## Bayesian models fitting
- We fit the 6 Bayesian model variants via `fitting_exp_alll.m`. This script is ready to run but it takes a long time; the results will be saved as `params_all_models_E2_Nsubj_22.mat`, which we provide.`fitting_exp_alll.m`  calls  `Fitting_model.m`, which finds the parameters that maximize the probability of the data given each model (maximum likelihood estimation) via `Loglike_bm_alll.m`.  `Fitting_model.m` loads the data, here `alldata_E2.mat`, as well as the psychometric curve fits, `params_psych_curves_exp2_m1_m2_all.mat` and outputs the best fitting parameters `params_all_models_E2_Nsubj_22.mat`.

## Bayesian models simulations
`sims_mu_like_vs_priorE.m` generates Figure S5 and `parameter_recovery_Bayesian_model_1_plot.m` generates Figure S7.

## Bayesian models simulations full
`fitting_all_model_sim_and_output_plotting.m` generates Figure S8.

## Pupillometry analyses 

- `regressions_pupilAlll.m` generates Figure 6 from our paper (once the full dataset will be made available). Note that it loads `dv_val_mat.mat` - generated by the `analysis_all.m` script described above, as well as `GLMEs_sliding_window_stim_locked_.mat` and `GLMEs_sliding_window_resp_locked_.mat`, generated by the script `regressions_pupil_sliding_window_analysis_FIN.m`. Note that due to their large size, the majority of the .mat files will need to be downloaded from [More on data and model fits](https://drive.google.com/drive/folders/1OW1x80jKBBn9jowLEeM6Y8xWxc_w7NO4?usp=drive_link). The output from the Eyelink eye tracker was parsed with variants of previous scripts: https://github.com/lianaan/Microsacc/tree/master.

## DDM
- While `analysis_all.m` loads `ddm_params_and_model_comp.mat` and `output_list_RT_allEE.mat` and `output_list_choices_allEE.mat` to generate Figure 7, here we show how we implemented the drift diffusion model fitting with [PyDDM](https://github.com/mwshinn/PyDDM) from Shinn et al, 2020.  
  
## Pilot control experiment data

- We also share the data from Exp. 2 control, which we used to generate Figure S9. This is `alldata_E2ctr.mat`, which is also Nsubj * Ncond, with Nsubj being 7. There are 4 conditions: 

- 1: No-Adapt-See 
- 2: No-Adapt-Bias
- 3: Adapt-See
- 4: Adapt-Bias


