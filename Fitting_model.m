function [params_fit_si_ci, nll_si_ci, params_fit, nll, params_fit_best, nll_both] = Fitting_model(sbji,modi)
curr_dir = pwd;
addpath([curr_dir, '/bads-master/'])

exp_i = 2;

load(['params_psych_curves_demo_exp',num2str(exp_i),'_m1_m2_all.mat']);
load(['alldata_demo_E',num2str(exp_i), '.mat']);

s_i = sbji;

mu_NA = params_psych_m1_all(s_i,1,1); 
mu_A = params_psych_m1_all(s_i,3,1); 
mus = [ mu_NA mu_A]; 


N_pts = 121;
N_samp = 800;

options = [];                       % Reset the OPTIONS struct
options.UncertaintyHandling = 1;    % Tell BADS that the objective is noisy
options.NoiseSize           = 1;    % Estimate of noise % confirmed
%options.NoiseFinalSamples   = 100;

if ismember(modi, [1 2 4])
    nvars = 4;
    if modi == 1
        LB = [NaN log(0.001) 0.50000 -10]; % mu is taken from psych curve fits
        UB = [NaN  log(0.15) 0.99999 10];
        PLB = [LB(1) LB(2)   0.5000 -3];
        PUB = [UB(1) UB(2)   0.9999 6.5];
    elseif modi == 2
        LB = [NaN log(0.001)   0.00001   0.50000 ]; % mu is taken from psych curve fits
        UB = [NaN  log(0.15)   0.99999 0.99999 ];
        PLB = [LB(1) LB(2)  0.4 0.50000   ];
        PUB = [UB(1) UB(2)  0.6 0.99999  ];
    elseif modi == 4
        LB = [NaN log(0.001) -20  0.50000  ]; % mu is taken from psych curve fits
        UB = [NaN  log(0.15)  20  0.99999  ];
        PLB = [LB(1) LB(2)  -5  0.5000   ];
        PUB = [UB(1) UB(2)  5  0.9999 ];
    end
elseif ismember(modi, [3 5])
    nvars = 5;
    if modi == 3
        LB = [NaN log(0.001)   0.00001  0.50000 -10 ]; % mu is taken from psych curve fits
        UB = [NaN  log(0.15)   0.99999  0.99999 10 ];
        PLB = [LB(1) LB(2)  0.4  0.5000 -3  ];
        PUB = [UB(1) UB(2)  0.6  0.9999 6.5];
    elseif modi == 5
        LB = [NaN log(0.001) -20  0.50000 -10 ]; % mu is taken from psych curve fits
        UB = [NaN  log(0.15)  20  0.99999 10 ];
        PLB = [LB(1) LB(2)  -5  0.5000 -3  ];
        PUB = [UB(1) UB(2)  5  0.9999 6.5];
    end
end

N_sp = 14;
params_fit = NaN(4,2,max(nvars));
nll = NaN(4,2);
n_grid = 101;
mi = 1; % without decision noise
for cii = 1:2 
   
    mu = mus(cii); 
    
    for ci = (2*cii-1): (2*cii)
        
        s = alldata(s_i,ci).stims;
        resp = alldata(s_i,ci).resp; 
        conf = alldata(s_i,ci).conf; 
        
        fun_nll = @(pars) -Loglike_bm_alll([mu pars],s, resp,conf,modi, mi,N_samp);
        
        params_fit_si = NaN(N_sp,nvars(mi));
        nll_si = NaN(N_sp,1);
        for si = 1:N_sp 
            start_pars = (PUB(2:nvars(mi))-PLB(2:nvars(mi))).*rand(1,nvars(mi)-1) + PLB(2:nvars(mi));
            [params_fit_si(si,1:length(start_pars)), nll_si(si), exitflag, output] = bads(fun_nll, start_pars, LB(2:nvars(mi)), UB(2:nvars(mi)), PLB(2:nvars(mi)), PUB(2:nvars(mi)), options);
        end
        nll_si_ci(ci,1:N_sp) = nll_si;
        params_fit_si_ci(ci,1:N_sp,1:nvars) = params_fit_si;
        indi_min = find(nll_si == min(nll_si));
        nll(ci,mi) = nll_si(indi_min(1));
        params_fit(ci,mi,1:length(start_pars)) = params_fit_si(indi_min(1),1:length(start_pars));
    end
    nll_both(cii) = nll(2*cii-1,mi)+ nll(2*cii,mi);
    
    params_fit_best(2*cii-1, 1:(1+length(start_pars))) = [mus(cii); squeeze(params_fit(2*cii-1,mi,1:length(start_pars)))];
    params_fit_best(2*cii, 1:(1+length(start_pars))) = [mus(cii); squeeze(params_fit(2*cii,mi,1:length(start_pars)))];
end
save(['fits_model_',num2str(modi),'_exp_',num2str(exp_i),'_sbji_',num2str(sbji) ,'.mat'], 'params_fit_si_ci', 'nll_si_ci','params_fit','nll', 'params_fit_best', 'nll_both',    '-mat')
end

