clear all; close all;

exp_i = 2;

%for sbji =  1
sbji = 1;

Ncond = 4;
ngrid = 201;%101; 
mu_range = linspace(-0.2,0.2,ngrid); 
sigma_range = exp(linspace(log(0.0001),log(0.2), ngrid)); 
lambda_range = linspace(0, 0.3,ngrid);

load(['alldata_demo_E',num2str(exp_i),'.mat'])
alldata = alldata(sbji,:);

loglike_max_L = NaN(ngrid);
[mu_est, sigma_est, lambda_est] = deal(NaN(size(alldata,1),Ncond));
n_trials = Ncond*121;
n_pars = 2*Ncond+1;

    for li = 1:ngrid
        lambda_val = lambda_range(li);
        for ci = 1:Ncond
            
            vals = alldata(sbji,ci).stims;
            resp = alldata(sbji,ci).resp;
            
            for mui = 1: ngrid
                mu_val = mu_range(mui);
                for si = 1:ngrid
                    sigma_val = sigma_range(si);
                    
                    likelihood_L(mui,si) = loglike(vals, resp , mu_val, sigma_val, lambda_val);
                end
            end
            [max_val(sbji,ci,li),ind_max]= max(likelihood_L(:));
            [imax, jmax] = ind2sub(size(likelihood_L),ind_max);
            mu_estl(ci,li) = mu_range(imax);
            sigma_estl(ci,li) = sigma_range(jmax);
            
        end
        loglike_max_L(li) = sum(max_val(sbji,:,li),2);
        
    end
    [loglike_max(sbji),li_max] = max(loglike_max_L(:));
    
    lambda_est(sbji) = lambda_range(li_max);
    mu_est(sbji,:) =  squeeze(mu_estl(:,li_max));
    sigma_est(sbji,:) = squeeze(sigma_estl(:,li_max));
    
    AIC(i) = -2*loglike_max(sbji)+2*n_pars ;
    BIC(i) = -2*loglike_max(sbji) + n_pars* log(n_trials);
%end % end of sbji loop

save(['psych_curves_fitting_m2_',num2str(ngrid),'_sbj_',num2str(sbji) ,'.mat'], 'loglike_max', 'AIC','BIC', 'mu_est','sigma_est','lambda_est', '-mat')
%end
