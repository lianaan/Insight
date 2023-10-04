function [effectSize,pValue,CI,bootstrapCI] = BCa_CI_Sr(Y,X,sig_level,numBootstraps)
%%%% This function calculates the bias-corrected and accelerated confidence
%%%% interval for an effect size between variables Y and X.
%%%% Ref: Efron B. Better bootstrap confidence intervals. Journal of the American statistical Association. 1987;82(397):171-185.
%%%% Written by Kenneth Wengler,PhD
%%%% Columbia University, Department of Psychiatry
%%%% May 15, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y: dependent variable (Nx1); N subjects
% X: independent varaible(s) (NxV); N subjects, V variables
% sig_level: alpha significnace level (e.g., 0.05 [95%])
% numBootstraps: number of bootstraps (e.g., 10,000)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set parameters
if nargin == 2
    sig_level = 0.05;
    numBootstraps = 100000;
elseif nargin == 3
    numBootstraps = 100000;
end


%% Exlcude NaNs
NaN_flag = (isnan(Y) + sum(isnan(X),2))>0;
Y = Y(~NaN_flag);
X = X(~NaN_flag,:);


%% Calculate effect size

[p,h, stats] = signrank(Y,X);
effectSize  = stats.zval;
pValue = p;



%% Get jacknife distribution
numJackknife = size(Y,1);
for jackknife = 1:size(Y,1)
    if jackknife == 1
        jackknife_index = [2:numJackknife];
    elseif jackknife == numJackknife
        jackknife_index = [1:numJackknife-1];
    else
        jackknife_index = [1:jackknife-1,jackknife+1:numJackknife];
    end
    Y_jk = Y(jackknife_index);
    X_jk = X(jackknife_index,:);
    
    [p,h, stats] = signrank(Y_jk,X_jk);
    effect_size_jackknife(jackknife) = stats.zval;
    
    
end


%% Get bootstrap distribution
for bootstrap = 1:numBootstraps
    bootstrap_index = randi(size(Y,1),size(Y));
    Y_bs = Y(bootstrap_index);
    X_bs = X(bootstrap_index,:);
    
    [p,h, stats] = signrank(Y_bs,X_bs);
    effect_size_bootstrap(bootstrap) = stats.zval;
    
end


%% Get confidence interval
theta_hat = sort(effect_size_bootstrap)';
theta = effectSize;
B = size(theta_hat,1);

z_hat_0 = norminv(sum(theta_hat < theta)/B);

theta_hat_dot = mean(effect_size_jackknife);

num = sum((theta_hat_dot - effect_size_jackknife).^3);
den = 6.*(sum((theta_hat_dot - effect_size_jackknife).^2).^(3/2));
a_hat = num ./ den;

num = z_hat_0 + norminv(sig_level/2);
den = 1 - a_hat .* (z_hat_0 + norminv(sig_level/2));
alpha_1 = normcdf(z_hat_0 + num./den);

num = z_hat_0 + norminv(1 - sig_level/2);
den = 1 - a_hat .* (z_hat_0 + norminv(1 - sig_level/2));
alpha_2 = normcdf(z_hat_0 + num./den);

lb = theta_hat(max(floor(alpha_1*B),1));
ub = theta_hat(min(ceil(alpha_2*B),B));
[CI] = [lb,ub];

bootstrap_lb = prctile(theta_hat,(sig_level/2)*100);
bootstrap_ub = prctile(theta_hat,(1-sig_level/2)*100);
bootstrapCI = [bootstrap_lb,bootstrap_ub];