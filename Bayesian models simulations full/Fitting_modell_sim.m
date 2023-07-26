function [params_fit_si_ci, nll_si_ci, params_fit, nll, params_fit_best] = Fitting_modell_sim(mpr_indices,modi)
curr_dir = pwd;
addpath([curr_dir, '/bads-master/'])

mu_like_set = [0 -0.05 -0.1 -0.15];
sigma_enc_set = [0.01 0.05 0.1];
sigma_A_set = [0.01 0.05 0.1];

mu_like = mu_like_set(mpr_indices(1));
sigma_enc = sigma_enc_set(mpr_indices(2));
sigma_A = sigma_A_set(mpr_indices(3));
runi = mpr_indices(4);

%rng(runi);

exp_i = 2;

load(['Asim_data_full_model', '_mu_like_', num2str(mu_like), '_sigma_enc_', num2str(sigma_enc), '_sigma_A_', num2str(sigma_A), '_rep_no_', num2str(runi) '.mat'])

mu_A = -0.055;
mus = mu_A;


N_pts = 121;
N_samp = 800;

options = [];                       % Reset the OPTIONS struct
options.UncertaintyHandling = 1;    % Tell BADS that the objective is noisy
options.NoiseSize           = 1;    % Estimate of noise % confirmed
%options.NoiseFinalSamples   = 100;
options.Display = 'off';

if ismember(modi, [1 2 4 6])
    nvars = 4;
    if modi == 1
        LB = [NaN log(0.0001) 0.50000 -10]; % mu is taken from psych curve fits
        UB = [NaN  log(0.5) 0.99999 10];
        PLB = [LB(1) log(0.001)   0.5000 -3];
        PUB = [UB(1) log(0.15)   0.9999 6.5];
    elseif modi == 2
        LB = [NaN log(0.001)   0.00001   0.50000 ]; % mu is taken from psych curve fits
        UB = [NaN  log(0.15)   0.99999 0.99999 ];
        PLB = [LB(1) LB(2)  0.4 0.50000   ];
        PUB = [UB(1) UB(2)  0.6 0.99999  ];
    elseif modi == 4 | modi == 6
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
params_fit = NaN(2,max(nvars));

nll = NaN;
n_grid = 101;
mi = 1; % without decision noise

cii = 1;
ci = 1;
mu = mus(cii);

s = stims;

fun_nll = @(pars) -Loglike_bm_alll([mu pars],s, resp,conf,modi, mi,N_samp);

params_fit_si = NaN(N_sp,nvars(mi));
nll_si = NaN(N_sp,1);
for si = 1:N_sp
    rng(si);
    start_pars = (PUB(2:nvars(mi))-PLB(2:nvars(mi))).*rand(1,nvars(mi)-1) + PLB(2:nvars(mi));
    [params_fit_si(si,1:length(start_pars)), nll_si(si), exitflag, output] = bads(fun_nll, start_pars, LB(2:nvars(mi)), UB(2:nvars(mi)), PLB(2:nvars(mi)), PUB(2:nvars(mi)), options);
end

nll_si_ci(1:N_sp) = nll_si;
params_fit_si_ci(1:N_sp,1:nvars) = params_fit_si;
indi_min = find(nll_si == min(nll_si));
nll(mi) = nll_si(indi_min(1));
params_fit(mi,1:length(start_pars)) = params_fit_si(indi_min(1),1:length(start_pars));

params_fit_best(2*cii-1, 1:(1+length(start_pars))) = [mus(cii) squeeze(params_fit(2*cii-1,1:length(start_pars)))];

save(['Afits_sim_full_model_', num2str(modi),'_mu_like_', num2str(mu_like), '_sigma_enc_', num2str(sigma_enc), '_sigma_A_', num2str(sigma_A), '_rep_no_', num2str(runi), '.mat'],'params_fit_si_ci', 'nll_si_ci','params_fit','nll', 'params_fit_best', '-mat' );
end



