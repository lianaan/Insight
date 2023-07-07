clear all; close all;
mu_A = -0.055;
k_confidence = 0.8; % as fit in Adapt-Believe

sigma_enc_set = [0.03 0.06 0.09 0.2];
mu_like_set = [0 -0.055 -0.15];
priors_set = [0.5 0.15 0.05];

Ntrials = 50000;
s = linspace(-0.3, 0.3, Ntrials);
resp_all = NaN(2,3,3, Ntrials);
conf_all = NaN(2,3,3, Ntrials);

N_samp = 1; 
nbinz = 11;
binz = [];
for j = 1:nbinz
    binz(j) = quantile(s, j/nbinz);
end
binz = [min(s)*1.001 binz ];
binz_pos = (binz(2:end)+binz(1:end-1))/2;

prop_cw_all = NaN(2,3,3, nbinz);
prop_cf_high_all = NaN(2,3,3, nbinz);


for modi = 1: 2
    
    for i1 = 1:4
        sigma_enc = sigma_enc_set(i1);
        for i2 = 1:3
            
 
            if modi == 1
                mu_like = mu_like_set(i2);
                
                mu_A_encode = mu_A;  
                sigma_s = sigma_enc;
                prob_right = 0.5;
                bd = k_confidence;
                mu_A_decide = mu_like;  
            elseif modi == 2
                prior_right = priors_set(i2);
                
                mu_A_encode = mu_A; 
                sigma_s = sigma_enc; 
                prob_right = prior_right; 
                bd = k_confidence; 
                mu_A_decide = 0; 
            end
            
            stims        =  repmat(s', [1,  N_samp]);
            % the center of the measurement is offset to stims-mu_A, mu_A could be 0
            x            =  normrnd(bsxfun(@minus, stims, mu_A_encode),sigma_s);
            
            
            dv1 = log(prob_right/(1-prob_right)) +  log(  (normcdf(-mu_A_decide+max(s)-x, 0, sigma_s)-normcdf(-mu_A_decide-x, 0, sigma_s))./...
                (normcdf(-mu_A_decide-x, 0, sigma_s)- normcdf(-mu_A_decide+min(s)-x,0, sigma_s)));
            
            post = nan(size(x,1), 2, size(x,2));
            post(:,1,:) = [1./(1+exp(-dv1))]; % 1 for CW, 1
            post(:,2,:) = [1./(1+exp(dv1))]; % 2 for CCW, 0
            
            post = bsxfun(@rdivide,post, sum(post,2)); % normalize to make sure that the 2 posterior values add up to 1

            [post_max,dd]       = max(post,[],2);
            r_pred = abs(squeeze(dd)-2);
            post_max = squeeze(post_max);
            c_pred = post_max > bd;
            
            resp_all(modi,i1,i2, 1:Ntrials) = r_pred;
            conf_all(modi,i1,i2, 1:Ntrials) = c_pred;
        
            
            for j = 1:(nbinz)
                
                indi = find(s>binz(j) & s<=binz(j+1) );
                if length(indi)>0
                    prop_cw_all(modi,i1,i2,j) = nansum(r_pred(indi))/sum(~isnan(r_pred(indi)));
                    prop_cf_high_all(modi,i1,i2,j) = nansum(c_pred(indi))/sum(~isnan(c_pred(indi)));
                    
                else
                    prop_cw_all(modi,i1,i2,j) = nan;
                    prop_cf_high_all(modi,i1,i2,j) = nan;
                    
                end
            end
        end
 
    
end
end
%%
close all;
marginsa = [0.08 0.03 0.09 0.1]; %MARGINS = [LEFT RIGHT BOTTOM TOP]
guttera = [0.06 0.06];
colorz = [255 35 0; 242 120 30;  254 196 79; 255 227 149]/255;

linewi = 1.5;

figure(15) %Figure S5 in paper
set(gcf, 'Position', [100 100 900 280])
for modi = 1: 2
    for ij = 1:3 % the 3 conditions, either for M1 or for M2
        tight_subplot(2,6,1, 3*(modi-1)+ij, guttera, marginsa)
        for kk = 1:4 % across the 3 noises
            plot(binz_pos, squeeze(prop_cw_all(modi,kk,ij,:)),  'Color', colorz(kk,:), 'Linewidth', linewi); hold on; %'o-', 'MarkerFaceColor', colorz(kk,:),'MarkerEdgeColor', colorz(kk,:),
        end
        plot(zeros(1,10), linspace(0,1,10), '--k'); hold on;
        box off
        set(gca, 'tickdir', 'out')
        xlim([-0.3 0.3])
        ylim([0 1])
        
       
        tight_subplot(2,6,2, 3*(modi-1)+ij, guttera, marginsa)
        for kk = 1:4 % across the 3 noises
            plot(binz_pos, squeeze(prop_cf_high_all(modi,kk,ij,:)),  'Color', colorz(kk,:), 'Linewidth', linewi ); hold on;
        end
        plot(zeros(1,10), linspace(0,1,10), '--k'); hold on;
        box off
        set(gca, 'tickdir', 'out')
        xlim([-0.3 0.3])
        ylim([0 1])
    end
    
    
end

%%
psname = 'simulations_mu_likelihood_vs_prior_model.pdf'
%print_pdf(psname)


