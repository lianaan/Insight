clear all; close all;

exp_i = 2;

Ncond = 4;
ngrid = 201;
mu_range = linspace(-0.2,0.2,ngrid);
sigma_range = exp(linspace(log(0.0001),log(0.2), ngrid));
lambda_range = linspace(0, 0.3,ngrid);

load(['alldata_E',num2str(exp_i),'.mat'])
%load(['alldata_E2ctr.mat'])

n_trials = Ncond*121;
n_pars = 12;

[max_val, mu_est, sigma_est, lambda_est] = deal(NaN(size(alldata,1),Ncond));
loglike_max = nan(size(alldata,1),1);

for sbji =  1: size(alldata,1)
    sbji
    data = alldata(sbji,:);
    
   
    for ci = 1:4
        ci
        likelihood = nan(ngrid,ngrid,ngrid);
        
        vals1 = data(ci).stims;
        resp1 = data(ci).resp;
        for mui = 1: ngrid
            mu_val = mu_range(mui);
            
            for si = 1:ngrid
                sigma_val = sigma_range(si);
                
                for li = 1:ngrid
                    lambda_val = lambda_range(li);
                    likelihood(mui,si,li) = loglike(vals1, resp1 , mu_val, sigma_val, lambda_val);
                end
            end
        end
        [max_val(sbji,ci),ind_max]= max(likelihood(:));
        [imax, jmax, kmax] = ind2sub(size(likelihood),ind_max);
        
        
        mu_est(sbji,ci) = mu_range(imax);
        sigma_est(sbji,ci) = sigma_range(jmax);
        lambda_est(sbji,ci) = lambda_range(kmax);
    end
    loglike_max(sbji) = sum(max_val(sbji,:),2);
    AIC(sbji) = -2*loglike_max(sbji)+2*n_pars ;
    BIC(sbji) = -2*loglike_max(sbji) + n_pars* log(n_trials);
end
save(['psych_curves_fitting_m1_',num2str(ngrid),'_E_', num2str(exp_i),'_all.mat'], 'loglike_max', 'AIC','BIC', 'mu_est','sigma_est','lambda_est', '-mat')
%save(['psych_curves_fitting_m1_',num2str(ngrid),'_E2ctr.mat'], 'loglike_max', 'AIC','BIC', 'mu_est','sigma_est','lambda_est', '-mat')

