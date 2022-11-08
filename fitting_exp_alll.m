clear all; close all;
indi_sel = [1];
exp_i = 2;

Ncond = 4;
N_sp = 14;
params_fit_si_ci_all = NaN(length(indi_sel), 5, Ncond, N_sp, 5);
nll_si_ci_all = NaN(length(indi_sel), 5,Ncond, N_sp);
params_fit_all = NaN(length(indi_sel), 5, Ncond, 2, 4); % no need
nll_all = NaN(length(indi_sel), 5,2);
params_fit_best_all =  NaN(length(indi_sel), 5, Ncond,5);
sbj_index_all = indi_sel;
AIC_all = NaN(length(indi_sel), 5);
BIC_all = NaN(length(indi_sel), 5);
 

for si = 1:length(indi_sel)
    sbji = indi_sel(si);
    for modi = 1:5
        [params_fit_si_ci, nll_si_ci, params_fit, nll, params_fit_best, nll_both] = Fitting_model(sbji,modi);
        if ismember(modi, [3 5])
        %params_fit_si_ci_all(si,modi,1:2,1:N_sp, 1:5) =   params_fit_si_ci;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     1:2,1:N_sp,1:5) = params_fit_si_ci;
        params_fit_best_all(si,modi,1:Ncond,1:5) = params_fit_best;
        else
            params_fit_si_ci_all(si,modi,1:Ncond,1:N_sp,1:4) = params_fit_si_ci;
            params_fit_best_all(si,modi,1:Ncond,1:4) = params_fit_best;
        end
        nll_si_ci_all(si,modi,1:Ncond,1:N_sp) = nll_si_ci; % nll across the 8 iterations
        nll_all(si,modi, 1:Ncond) = nll(1:Ncond,1)';
        nll_sum_all(si,modi) = sum(nll_both);
        
        
    end
end
%}
%% load the values saved from simulations, and save them
for si = 1:length(indi_sel) 
    sbji = indi_sel(si);
    for modi = 1:5
        %[params_fit_si_ci, nll_si_ci, params_fit, nll, params_fit_best, nll_both] = Fitting_models_all_online(sbji,modi);
        if isfile(['fits_model_', num2str(modi), '_exp_',num2str(exp_i),'_sbji_', num2str(sbji), '.mat'])
        load(['fits_model_', num2str(modi), '_exp_',num2str(exp_i),'_sbji_', num2str(sbji), '.mat'])
        
        if ismember(modi, [3 5])
            params_fit_si_ci_all(si,modi,1:Ncond,1:N_sp,1:5) = params_fit_si_ci;
            params_fit_best_all(si,modi,1:Ncond,1:5) = params_fit_best;
        elseif ismember(modi, [1 2 4])
            params_fit_si_ci_all(si,modi,1:Ncond,1:N_sp,1:4) = params_fit_si_ci(1:Ncond, 1:N_sp, 1:4);
            params_fit_best_all(si,modi,1:Ncond,1:4) = params_fit_best;
        end
        nll_si_ci_all(si,modi,1:Ncond,1:size(nll_si_ci,2)) = nll_si_ci; % nll across the 8 iterations
        nll_all(si,modi, 1:Ncond) = nll(1:Ncond,1)';
        nll_sum_all(si,modi) = sum(nll_both);
        if ismember(modi, [3 5])
            AIC_all(si,modi) = 2* nll_sum_all(si,modi) + 2* (Ncond*(5-1));
            BIC_all(si,modi) = 2* nll_sum_all(si,modi) + (Ncond*(5-1)) * log(Ncond*61);
        else
            AIC_all(si,modi) = 2* nll_sum_all(si,modi)+ 2* (Ncond*(4-1));
            BIC_all(si,modi) = 2* nll_sum_all(si,modi) + (Ncond*(4-1)) * log(Ncond*61);
        end
        %else
         %   params_fit_si_ci_all(si,modi,1:Ncond,1:N_sp,1:5) = NaN;
        end
    end
end
%%
save(['params_all_models_E2_Nsubj_', num2str(size(nll_sum_all,1)), '.mat'],'AIC_all', 'BIC_all','nll_sum_all','nll_all', 'nll_si_ci_all', 'params_fit_best_all','params_fit_si_ci_all', '-mat' )


