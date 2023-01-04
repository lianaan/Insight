function prc = Predict_bm_alll(params,s, resp,conf,modi, mi,N_samp)
if modi == 1
    mu_A_encode = params(1); % bias of the measurement, could be 0
    sigma_s = exp(params(2)); % sensory noise
    prob_right = 0.5;
    bd = params(3); %  k_confidence / confidence boundary
    mu_A_decide = params(4); % mu likelihood
elseif modi == 2
    mu_A_encode = params(1); % bias of the measurement, could be 0
    sigma_s = exp(params(2)); % sensory noise
    prob_right = params(3); %0.5;
    bd = params(4); % k_confidence
    mu_A_decide = 0; % mu likelihood
elseif modi == 3
    mu_A_encode = params(1); % bias of the measurement, could be 0
    sigma_s = exp(params(2)); % sensory noise
    prob_right = params(3); %0.5;
    bd = params(4); % k_confidence
    mu_A_decide = params(5); % mu likelihood
elseif modi == 4
    mu_A_encode = params(1); % bias of the measurement, could be 0
    sigma_s = exp(params(2)); % sensory noise
    prob_right = 0.5;
    mu_R = params(3); %  k_choice / response bias
    bd = params(4);% k_confidence
    mu_A_decide = 0; % mu likelihood
elseif modi == 5
    mu_A_encode = params(1); % bias of the measurement, could be 0
    sigma_s = exp(params(2)); % sensory noise
    prob_right = 0.5;
    mu_R = params(3); %  k_choice / response bias
    bd = params(4); % k_confidence
    mu_A_decide = params(5); % mu likelihood
end


Ntrials = length(s);

%encoding
% simulate how the sensory noise would alter the encoding of stimuli s into the measurement x
% by repeating this encoding several times, N_samp times

stims        =  repmat(s', [1,  N_samp]);
% the center of the measurement is offset to stims-mu_A_encode, mu_A_encode could be 0
x            =  normrnd(bsxfun(@minus, stims, mu_A_encode),sigma_s);

%decision model
%1 stim, 2 decision variables dv1 dv2. dv2 is just 1/dv1
% derivation as in doc and Will Adler thesis . dv1 = p(C= 1(CW)| x)/ p(C=
% -1(CCW)|x)
dv1 = log(prob_right/(1-prob_right)) +  log(  (normcdf(-mu_A_decide+max(s)-x, 0, sigma_s)-normcdf(-mu_A_decide-x, 0, sigma_s))./...
    (normcdf(-mu_A_decide-x, 0, sigma_s)- normcdf(-mu_A_decide+min(s)-x,0, sigma_s)));



post = nan(size(x,1), 2, size(x,2));
%if mi == 1 % wo decision noise
post(:,1,:) = [1./(1+exp(-dv1))]; % 1 for CW, 1
post(:,2,:) = [1./(1+exp(dv1))]; % 2 for CCW, 0
if mi == 2  % w decision noise
    postN = post.^alpha;
elseif mi == 1
    postN = post;
end

post = bsxfun(@rdivide,post, sum(post,2)); % normalize to make sure that the 2 posterior values add up to 1
postN = bsxfun(@rdivide,postN, sum(postN,2));

if mi == 1
    
    [post_max,dd]       = max(post,[],2);
    if ismember(modi, [1 2 3 ])
        r_pred = abs(squeeze(dd)-2);
    elseif ismember(modi, [4 5])
        r_pred = dv1 > mu_R; % only for models 4 and 5, with k_choice
    end
    post_max = squeeze(post_max);
    c_pred = post_max > bd;
    
elseif mi == 2
    
    for ti = 1: Ntrials
        for si = 1:N_samp
            ddN(ti,si) = find(mnrnd(1, postN(ti,:,si)) == 1);
            c_pred(ti,si) = postN(ti,ddN(ti,si),si) > bd;
        end
    end
    r_pred = abs(ddN - 2);
end

for ti = 1:Ntrials
    prc(ti,1) = sum(r_pred(ti,:) ==1 & c_pred(ti,:) == 1)/N_samp;
    prc(ti,2) = sum(r_pred(ti,:) ==1 & c_pred(ti,:) == 0)/N_samp;
    prc(ti,3) = sum(r_pred(ti,:) ==0 & c_pred(ti,:) == 1)/N_samp;
    prc(ti,4) = sum(r_pred(ti,:) ==0 & c_pred(ti,:) == 0)/N_samp;
end

prc(prc==0)= 1/N_samp;
prc(prc==1) = 1-1/N_samp;

end
