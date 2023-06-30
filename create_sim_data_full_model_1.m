clear all; close all;
mu_A = -0.055;

mu_like_set = [0 -0.05 -0.1 -0.15]; 
sigma_enc_set = [0.01 0.05 0.1];
sigma_A_set = [0.01 0.03 0.05];

Ntrials = 121;%1000;
prob_right = 0.5;
k_confidence = 0.9;
for i1 = 1:4 % 4 posssible levels of mu_like
    for i2 = 1:3 % 3 posssible levels of sigma_enc
        for i3 = 1:3 % 3 posssible levels of sigma_A
            for i4 = 1:5 % % 5 simulated datasets for each unique parameter combo 
               
                mu_like = mu_like_set(i1);
                sigma_enc = sigma_enc_set(i2);
                sigma_A = sigma_A_set(i3);
                repi = i4;
                
                A = normrnd(mu_A, sigma_A );
                stims        =  linspace(-0.3, 0.3, Ntrials);
                x            =  normrnd(bsxfun(@minus, stims, A),sigma_enc);
                
                for xi = 1:length(x)
                    xval = x(xi);
                    
                    funn1 = @(s,Av) exp((-(xval-(s-Av)).^2)/(2*sigma_enc^2)).* exp( (-(Av-mu_like).^2)/(2*sigma_A^2));
                    A = normrnd(mu_like, sigma_A, 1, Ntrials);
                    Amin = min(A);
                    Amax = max(A);
                    integr1_1(xi) = integral2(funn1 ,0,max(stims),Amin,Amax);
                    integr2_1(xi) = integral2(funn1 ,min(stims),0,Amin,Amax);
                    dv_sigma(xi) = log(prob_right/(1-prob_right)) + log((1./(max(stims)-0))./(1./(0-min(stims))) * (integr1_1(xi)./integr2_1(xi)));
                    
                end
                
                post(:,1,:) = [1./(1+exp(-dv_sigma))]; % 1 for CW, 1
                post(:,2,:) = [1./(1+exp(dv_sigma))]; % 2 for CCW, 0
                
                
                post = bsxfun(@rdivide,post, sum(post,2)); % normalize to make sure that the 2 posterior values add up to 1
                [post_max,dd]       = max(post,[],2);
                
                r_pred = abs(squeeze(dd)-2);
                resp = r_pred;
                post_max = squeeze(post_max);
                c_pred = post_max > k_confidence;
                conf = c_pred;
                
                filename = ['Asim_data_full_model', '_mu_like_', num2str(mu_like), '_sigma_enc_', num2str(sigma_enc), '_sigma_A_', num2str(sigma_A), '_rep_no_', num2str(repi) '.mat'];
                save(filename,'stims', 'x', 'resp', 'conf', '-mat')
            end
        end
    end
end

