function LL =loglike(Vals, Resp, mu, sigma, lambda)
%prc : proportion response clockwise

prc=nan(length(mu),length(sigma),length(lambda),length(Vals));

s_tab(1,1,1,1:length(Vals(~isnan(Vals)))) = Vals(~isnan(Vals));
prc(:,:,:,1:length(Vals(~isnan(Vals))))=function_psi(s_tab,mu,sigma, lambda);

%replace 0 with 1/Ntrials
prc(prc==0) = 1/sum(~isnan(Vals));
prc(prc==1) = 1-1/sum(~isnan(Vals));

Respp = Resp(~isnan(Vals));

%calculate loglike as sum across trials
LL  = nansum(log(prc(:,:,:,Respp==1)),4) + nansum(log(1-prc(:,:,:,Respp==0)),4); 

