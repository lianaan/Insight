clear all; close all;

indi_sel = 1:7%1:22;%[1];
exp_i = 2;
Ntrials = 121;  

Nmodels = 6;

Ncond = 4;
N_sp = 14; % 14 tries to find the global minima for nll out of several local minima found by the MLE optimization
params_fit_si_ci_all = NaN(length(indi_sel), Nmodels, Ncond, N_sp, 5);
nll_si_ci_all = NaN(length(indi_sel), Nmodels,Ncond, N_sp);
nll_all = NaN(length(indi_sel), Nmodels,2);
params_fit_best_all =  NaN(length(indi_sel), Nmodels, Ncond,5);
sbj_index_all = indi_sel;
AIC_all = NaN(length(indi_sel), Nmodels);
BIC_all = NaN(length(indi_sel), Nmodels);

%%
for si = 1:length(indi_sel)
    sbji = indi_sel(si)
    for modi = 1:Nmodels
        [params_fit_si_ci, nll_si_ci, params_fit, nll, params_fit_best, nll_both] = Fitting_model_ctr(sbji,modi);
        if ismember(modi, [3 5])
            %params_fit_si_ci_all(si,modi,1:Ncond,1:N_sp,1:5) = params_fit_si_ci;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     1:2,1:N_sp,1:5) = params_fit_si_ci;
            params_fit_best_all(si,modi,1:Ncond,1:5) = params_fit_best;
        else
            params_fit_si_ci_all(si,modi,1:Ncond,1:N_sp,1:4) = params_fit_si_ci;
            params_fit_best_all(si,modi,1:Ncond,1:4) = params_fit_best;
        end
        nll_si_ci_all(si,modi,1:Ncond,1:N_sp) = nll_si_ci; % nll across the 14 iterations
        nll_all(si,modi, 1:Ncond) = nll(1:Ncond,1)';
        nll_sum_all(si,modi) = sum(nll_both);
        
    end
end


%% load the values saved from simulations, and save them

curr_dir = pwd;
addpath([curr_dir, '/0fitsSctr/'])
cd([curr_dir, '/0fitsSctr/'])

for si = 1:length(indi_sel)
    sbji = indi_sel(si);
    
    for modi = 1:Nmodels 
        
        if isfile(['fitsSctr_model_', num2str(modi), '_exp_',num2str(exp_i),'_sbji_', num2str(sbji), '.mat'])
            load(['fitsSctr_model_', num2str(modi), '_exp_',num2str(exp_i),'_sbji_', num2str(sbji), '.mat'])
            
            if ismember(modi, [3 5])
                params_fit_si_ci_all(si,modi,1:Ncond,1:N_sp,1:5) = params_fit_si_ci;
                params_fit_best_all(si,modi,1:Ncond,1:5) = params_fit_best;
            elseif ismember(modi, [1 2 4 6 ])
                params_fit_si_ci_all(si,modi,1:Ncond,1:N_sp,1:4) = params_fit_si_ci(1:Ncond, 1:N_sp, 1:4);
                params_fit_best_all(si,modi,1:Ncond,1:4) = params_fit_best;
            end
            nll_si_ci_all(si,modi,1:Ncond,1:size(nll_si_ci,2)) = nll_si_ci; % nll across the 8 iterations
            nll_all(si,modi, 1:Ncond) = nll(1:Ncond,1)';
            
            nll_sum_all(si,modi) = sum(nll_both);
            if ismember(modi, [3 5])
                AIC_all(si,modi) = 2* nll_sum_all(si,modi) + 2* (Ncond*(5-1));
                BIC_all(si,modi) = 2* nll_sum_all(si,modi) + (Ncond*(5-1)) * log(Ncond*Ntrials);
            else
                AIC_all(si,modi) = 2* nll_sum_all(si,modi)+ 2* (Ncond*(4-1));
                BIC_all(si,modi) = 2* nll_sum_all(si,modi) + (Ncond*(4-1)) * log(Ncond*Ntrials);
            end
            
        end     
    end
end
%%
save ('nll_models_E2ctr_1_6_Nsubj_', ', num2str(length(indi_sel)), '.mat', 'nll_all', '-mat')
save(['params_all_models6_E2ctr_Nsubj_', num2str(length(indi_sel)), '.mat'],'AIC_all', 'BIC_all','nll_sum_all','nll_all', 'nll_si_ci_all', 'params_fit_best_all','params_fit_si_ci_all', '-mat' )
%}
