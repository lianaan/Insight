clear all; close all;

%exp_i = 1;
exp_i = 2;
if exp_i == 1
    sbj_list = {'1','2','3','4','5','6','7','8','9','10', '11', '12','13','14','15','16','17','18','19','20','21','22'};
    adapt_type_all = [ 1 -1 -1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1];
elseif exp_i == 2
    sbj_list = {'1','2','3','4','5','6','7','8'}; %'9','10', '12','13','14','15','16','17','18','19','20','21','22','23'};
    % adapt_type_all = -1 for everyone
end

Nsubj = length(sbj_list);

curr_dir = pwd;
Ncond = 4;
Nfiles_per_cond = 2;
if exp_i == 1
    model_pred = 1;
    model_pred_bm = 0;
elseif exp_i == 2
    %mi = 1;  % select the model
    % mi = 1: Bayesian insight model
    % mi = 2: prior model
    % mi = 3: insight + prior model
    % mi = 4: response bias k_choice model
    mi = 5; % response bias + insight model % winning model for the control expt
    % mi = 6: response bias k_choice_and_confidence model
    model_pred = 1;
    model_pred_bm = 1;
end

load(['alldata_E2ctr.mat']);
Ntrials = length(alldata(1,1).stims);
%accuracy in cond Adapt-Believe
for si = 1: length(sbj_list)
    acc_A_BEL(si) = sum((alldata(si,4).stims>0 & alldata(si,4).resp == 1) |...
        (alldata(si,4).stims<0 & alldata(si,4).resp == 0))/ length(alldata(si,4).stims);
end

if exp_i == 1 & model_pred
    load('psych_curves_fitting_m2_101_E1.mat')
    params_psych_m2_all = NaN(Nsubj, Ncond,3);
    for ci = 1:Ncond
        params_psych_m2_all(:,ci,1) = mu_est_all(:,ci);
        params_psych_m2_all(:,ci,2) = sigma_est_all(:,ci);
        if ci == 1
            params_psych_m2_all(:,ci,3) = lambda_est_all;
        end
    end
    
elseif exp_i == 2 & (model_pred | model_pred_bm)
    
    load('psych_curves_fitting_m1_201_E2ctr.mat')
    params_psych_m2_all = NaN(Nsubj, Ncond,3);
    
    for ci = 1:Ncond
        params_psych_m2_all(:,ci,1) = mu_est_all(:,ci);
        params_psych_m2_all(:,ci,2) = sigma_est_all(:,ci);
        if ci == 1
            params_psych_m2_all(:,ci,3) = lambda_est_all(:,ci);%lambda_est_all;
        end
    end
    
    load(['params_all_models6_E2ctr_Nsubj_8.mat'])
    params_bm_all = squeeze(params_fit_best_all(:,mi,:,:));
    params_bm_allV = params_bm_all;
    params_bm_allV(:,:,2) =  exp(params_bm_all(:,:,2));% params_bm_all(:,:,2);
    
end

if model_pred_bm
    N_samp = 500;
end

if exp_i == 2
    insight = abs((squeeze(params_bm_all(:,4,4))));
end

nbinz = 11;
prop_cw_all = nan(Nsubj,Ncond,nbinz);
prop_cw_w_all = nan(Nsubj,Ncond,nbinz);
prop_cw_pr_all = NaN(Nsubj,Ncond, nbinz);
prop_cw_pred_all = NaN(Nsubj,2, Ncond, nbinz);
prop_cf_pred_all = NaN(Nsubj,2, Ncond ,nbinz);

for si = 1 : Nsubj
    
    datac = alldata(si,:);
    
    for ci = 1:Ncond
        if exp_i == 2
            rt_all_trials(si,ci,1:121) = [datac(ci).resp_times];
            
            %compute the decision variable
            mu_A = params_bm_all(si,ci,1);
            sigma_A = exp(params_bm_all(si,ci,2));
            
            prob_right = 0.5;
            mu_A_loglike =  params_bm_all(si,ci,4);
            
            for sampi = 1: N_samp
                x            =  normrnd(bsxfun(@minus, datac(ci).stims, mu_A), sigma_A);
                dv_val_sampi(si,ci,1:121, sampi) = log(prob_right/(1-prob_right)) +  log(  (normcdf(-mu_A_loglike+0.3-x, 0, sigma_A)-normcdf(-mu_A_loglike-x, 0, sigma_A))./...
                    (normcdf(-mu_A_loglike-x, 0, sigma_A)- normcdf(-mu_A_loglike+(-0.3)-x,0, sigma_A)));
                
            end
            dv_val_sampi(dv_val_sampi>2000) = 2000;
            dv_val_sampi(dv_val_sampi<-2000) = -2000;
            
            dv_val(si,ci,1:121)= nanmean(squeeze(dv_val_sampi(si,ci,1:121,:)),2);
        end
        
        binz = [];
        for j = 1:nbinz
            binz(j) = quantile(datac(ci).stims, j/nbinz);
        end
        
        binz = [min(datac(ci).stims)*1.001 binz ];
        binz_posE(si,ci,:) = (binz(2:end)+binz(1:end-1))/2;
        binzz(si,ci,:) = binz;
        
        for j = 1:(nbinz)
            indo = (datac(ci).stims>binz(j) & datac(ci).stims<=binz(j+1));
            indi = find(datac(ci).stims>binz(j) & datac(ci).stims<=binz(j+1) );
            if length(indi)>0
                prop_cw_all(si,ci,j) = nansum(datac(ci).resp(indi))/sum(~isnan(datac(ci).resp(indi)));
                prop_cf_high_all(si,ci,j) = nansum(datac(ci).conf(indi))/sum(~isnan(datac(ci).conf(indi)));
                rt_all(si,ci,j) = nanmedian(datac(ci).resp_times(indi)); %nanmedian(datac(ci).resp_times(indi)); % maybe nanmean? % before it was median
                if exp_i == 2
                    dv_val_binz(si,ci,j) = median(dv_val(si,ci,indi));
                end
            else
                prop_cw_all(si,ci,j) = nan;
                prop_cf_high_all(si,ci,j) = nan;
                rt_all(si,ci,j) = nan;
                if exp_i == 2
                    dv_val_binz(si,ci,j) = nan;
                end
            end
            indi_right = (datac(ci).resp == 1);
            indi_left = (datac(ci).resp == 0);
            if length(indo & indi_right)>5
                rt_right(si,ci,j) = nanmedian(datac(ci).resp_times(indo & indi_right));
            else
                rt_right(si,ci,j) = NaN;
            end
            if length(indo & indi_left)>5
                rt_left(si,ci,j) = nanmedian(datac(ci).resp_times(indo & indi_left));
            else
                rt_left(si,ci,j) = NaN;
            end
        end
        
        if model_pred
            prop_cw_pr_all(si,ci,:) = function_psi(binz_posE(si,ci,:),...
                params_psych_m2_all(si,ci,1), params_psych_m2_all(si,ci,2), params_psych_m2_all(si,1,3));
            
        end
        
        % adjustments for participants in experiment 1 who saw the adaptor spiral rotating CW
        if exp_i == 1 & adapt_type_all(si) == 1 & ci >= 3
            
            binz = binz - 2*(params_psych_m2_all(si,ci,1));
            binz_posE(si,ci,:) = (binz(2:end)+binz(1:end-1))/2;
            binzz(si,ci,:) = binz;
            
            params_psych_m2_all(si,ci,1) = - params_psych_m2_all(si,ci,1);
            prop_cw_pr_all(si,ci,:) = function_psi(binz_posE(si,ci,:),...
                params_psych_m2_all(si,ci,1), params_psych_m2_all(si,ci,2), params_psych_m2_all(si,1,3));
        end
        
        
        if model_pred_bm
            
            % depends on the model, figure out which predict we use
            %Predict_bm_alll(params,s, resp,conf,modi, mi,N_samp)
            
            prc = Predict_bm_alll(squeeze(params_bm_all(si,ci,1:end)),squeeze(binz_posE(si,ci,:))', datac(ci).resp, datac(ci).conf, mi,1,N_samp);
            
            prop_cw_pred_all(si,mi,ci,:) = prc(:,1)+ prc(:,2);
            prop_cf_pred_all(si,mi,ci,:) = prc(:,1)+ prc(:,3);
            
        end
    end
    
end


%%

colorz = [41 138 8; 152 191 100; 0 0 255; 137 195 255 ]/255;
colorz_shade = (colorz + 3*repmat([228 228 228], 4, 1)/255)/4;
marginsa = [0.08 0.03 0.09 0.1]; %MARGINS = [LEFT RIGHT BOTTOM TOP]
guttera = [0.06 0.06];

sz = [7 6 5 4];
eb_w=0.01;%errorbar width
eb_t=1.2; %errorbar line thickness
tlen1 = 0.02;
tlen2 = 0.02;
fontsz = 12;
binz_pos = binz_posE ;%/0.0167;%settings.ifi;
dashed_linez = 0;

indi_sel = 1:1:Nsubj;
msz_all = [5.5 5.5 5.5 5.5];



%%%%%%%%%%%%%% PSYCH CURVES

figure(3)
set(gcf, 'Position', [100 100 740 400])

tight_subplot(2,9,exp_i,[1 2 3], guttera, marginsa)
for  ci = 1:Ncond
    
    he(ci) = fill([mean(squeeze(binz_pos(indi_sel,ci,:)),1) mean(squeeze(binz_pos(indi_sel,ci,end:-1:1)),1)],   ...
        [mean(squeeze(prop_cw_pr_all(indi_sel,ci,:)),1)-std(squeeze(prop_cw_pr_all(indi_sel,ci,:)),1)./sqrt(length(indi_sel))...
        fliplr(mean(squeeze(prop_cw_pr_all(indi_sel,ci,:)),1) + std(squeeze(prop_cw_pr_all(indi_sel,ci,:)),1)./sqrt(length(indi_sel))) ],colorz_shade(ci,:), 'EdgeColor', 'None'); hold on;
    plot(mean(squeeze(binz_pos(indi_sel,ci,:)),1), nanmean(squeeze(prop_cw_all(indi_sel,ci,:)),1),'o','MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:) ); hold on;
    errorbar(mean(squeeze(binz_pos(indi_sel,ci,:)),1), nanmean(squeeze(prop_cw_all(indi_sel,ci,:)),1), nanstd(squeeze(prop_cw_all(indi_sel,ci,:)),1)/sqrt(length(indi_sel)), 'Color',  colorz(ci,:),'Linestyle', 'none','Linewidth',1.5,'CapSize',0); hold on;
    
end
box off
xlim_min = min(min(min(binz_pos)));
xlim_max = max(max(max(binz_pos)));
ylim_min = 0;
ylim_max = 1.01;
if exp_i == 1
    xlim_min = -0.5; xlim_max = 0.5;
elseif exp_i == 2
    xlim_min = -0.3; xlim_max = 0.3;
end
xlim([xlim_min xlim_max])
ylim([ylim_min ylim_max])
set(gca, 'FontSize', fontsz)
plot(zeros(1,nbinz), linspace(0,1,nbinz), '--k'); hold on;
box off
set(gca, 'tickdir', 'out')
set(gca, 'FontSize', fontsz)
set(gca, 'ticklength',[tlen1 tlen2])
set(gca, 'ytick', [0:0.25:1])
set(gca, 'xtick',[-0.5:0.1:0.5])
if exp_i == 2
    xlabel('test stimulus speed clockwise (a.u.)', 'FontName', 'Helvetica', 'FontSize', fontsz)
end
ylabel('proportion response clockwise','FontName', 'Helvetica', 'FontSize', fontsz)

psych_pars_list = {'\mu', '\sigma', '\lambda'};
%since we have 4 planned comparisons
alpha_sidak = 1-(1-0.05).^(1/4); % 0.0127
mszi = 2;
pars_scatter = 0.15;
ylim_minB = -0.14;


nboot = 500;%5000 for final results  %for 95% bootstrapped confidence intervals
ci_bnd_low = 0.025;
ci_bnd_high = 0.975;
if exp_i == 1 | exp_i == 2
    for pi = 1:2
        
        tight_subplot(2,9,exp_i,[3+2*pi-1 4+2*pi-1], guttera, marginsa)
        
        for ci = 1: 4
            
            paraz = squeeze(params_psych_m2_all(indi_sel,ci,pi));
            
            for kk = 1:nboot
                sample = randsample(paraz,length(indi_sel),1);
                paraz_samples(kk) = nanmedian(sample);
            end
            bci_paraz = [quantile(paraz_samples,ci_bnd_low); quantile(paraz_samples,ci_bnd_high)];
            
            hbp(ci)=bar(5*(pi-1)+ci, median(paraz), 'FaceColor', 'None','EdgeColor',colorz(ci,:)); hold on;
            plot((5*(pi-1)+ci)*ones(1, length(indi_sel))- pars_scatter+ 2*pars_scatter*rand(1,length(indi_sel)),  params_psych_m2_all(indi_sel,ci,pi), 'o','MarkerSize',mszi, 'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:)); hold on;
            
            he(ci) = errorbar(5*(pi-1)+ci, median(paraz), median(paraz)- bci_paraz(1),bci_paraz(2)-median(paraz),'CapSize',0); hold on;
            
            set(he(ci),'Color', [0 0 0], 'Linewidth',2)
            errorbarT(he(ci),eb_w,eb_t)
            set(hbp(ci), 'BarWidth', 0.8, 'Linewidth', 1)
        end
        if pi == 1
            ylim([-0.3 0.3])
        else
            ylim([0 0.3])
        end
        box off
        set(gca, 'tickdir', 'out')
        set(gca, 'xtick', [])
        set(gca,'ticklength', [tlen1 tlen2])
        set(gca, 'xticklabels', {})
        %title(psych_pars_list{pi})
        xlim([5*(pi-1)+0.2 5*(pi-1)+4.8])
        set(gca, 'xtick', [2.5 7.5 10.5])
        if pi == 1
            ylabel('fitted parameter values','FontName', 'Helvetica', 'FontSize', fontsz)
        end
        set(gca, 'FontSize', fontsz)
    end
end


if exp_i == 1 | exp_i == 2 % show median and 95% bootstrapped confidence intervals
    
    
    for pi =  3 % psych curve parameter 3
        
        ci = 1;
        paraz = squeeze(params_psych_m2_all(indi_sel,ci,pi));
        
        for kk = 1:nboot
            sample = randsample(paraz,length(indi_sel),1);
            paraz_samples(kk) = nanmedian(sample);
        end
        bci_paraz = [quantile(paraz_samples,ci_bnd_low); quantile(paraz_samples,ci_bnd_high)];
        
        tight_subplot(2,9,exp_i,[8], guttera, marginsa)
        hbp(ci)=bar(5*(pi-1)+ci, median(paraz), 'FaceColor', 'None','EdgeColor',[0.5 0.5 0.5]); hold on;
        plot((5*(pi-1)+ci)* ones(1, length(indi_sel))- pars_scatter + 2*pars_scatter*rand(1,length(indi_sel)),  params_psych_m2_all(indi_sel,ci,pi), 'o','MarkerSize',mszi, 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.5 0.5 0.5]); hold on;
        he(ci) = errorbar(5*(pi-1)+ci, median(paraz),median(paraz)- bci_paraz(1),bci_paraz(2)-median(paraz) ,'CapSize',0); hold on;
        set(he(ci),'Color', [0 0 0], 'Linewidth',2)
        errorbarT(he(ci),eb_w,eb_t)
        set(hbp(ci), 'BarWidth', 0.8, 'Linewidth', 1)
        set(gca, 'FontSize', fontsz)
        box off
        set(gca, 'tickdir', 'out')
        set(gca,'ticklength', [tlen1 tlen2])
        set(gca, 'xticklabels', {})
        ylim([0 0.21])
        
    end
end

compensation = -(squeeze(params_psych_m2_all(:,3,1))-squeeze(params_psych_m2_all(:,4,1)));

[p,h, stats] = signrank(compensation(indi_sel))
compensation_baseline = (squeeze(params_psych_m2_all(:,1,1))-squeeze(params_psych_m2_all(:,2,1)));
[p,h] = signrank(compensation_baseline(indi_sel));

compensation_norm = compensation./abs(squeeze(params_psych_m2_all(:,3,1)))
[p, h, stats ]= signrank(compensation_norm(indi_sel)-1);

mid_green = (colorz(1,:)+colorz(2,:))/2;
mid_blue = (colorz(3,:)+colorz(4,:))/2;


for kk = 1:nboot
    sample = randsample(compensation_norm,length(indi_sel),1);
    comp_norm_samples(kk) = nanmedian(sample);
end
bci_comp_norm = [quantile(comp_norm_samples,ci_bnd_low); quantile(comp_norm_samples,ci_bnd_high)];

tight_subplot(2,9,exp_i,[9], guttera, marginsa)
plot(linspace(0.5,2.5,20), zeros(1,20), '--k'); hold on;
plot(1*ones(1,length(indi_sel))-1/8+1/4*rand(1,length(indi_sel)),compensation_norm(indi_sel), 'o','MarkerSize', mszi, 'MarkerEdgeColor',mid_blue,'MarkerFaceColor',mid_blue); hold on;
xlim([0.5 1.5])
b22 = bar(1, median(squeeze(compensation_norm(indi_sel))), 'FaceColor', 'None','EdgeColor',mid_blue); hold on;
set(b22, 'BarWidth', 0.4, 'Linewidth', 1)
hee2 = errorbar(1, median(squeeze(compensation_norm(indi_sel))),median(squeeze(compensation_norm(indi_sel)))- bci_comp_norm(1),bci_comp_norm(2)-median(squeeze(compensation_norm(indi_sel))) ,'CapSize',0); hold on;
set(hee2,'Color', [0 0 0], 'Linewidth',2)
errorbarT(hee2,eb_w,eb_t)
set(gca, 'FontSize', fontsz)
set(gca, 'ticklength',[tlen1 tlen2])
box off
set(gca, 'tickdir', 'out')
set(gca, 'xtick',[1 2])
set(gca, 'xticklabel',{  ''},'FontName', 'Helvetica', 'FontSize', fontsz)
ylim([-2 6])
ylabel('  MAE compensation', 'FontName', 'Helvetica', 'FontSize', fontsz)
[h,p,stats] = ttest(compensation(indi_sel), compensation_baseline(indi_sel))
plot(linspace(0.5,2.5,10),1*ones(1,10), '-k', 'LineWidth',1.5); hold on;
%%

psname =['Fig3_manuscript_nboot_', num2str(nboot),'_exp_',num2str(exp_i),'CTR.pdf']
%print_pdf(psname)


%% prep for confidence plots for figure 4 and supp, based on Maldonado Moscoso et al 2020


gauss_fit_just_mu = @(mu) sum((a_val*exp(-((s0-mu)/c_val).^2)-s1).^2);
options = [];                       % Reset the OPTIONS struct
options.UncertaintyHandling = 1;    % Tell BADS that the objective is noisy
options.NoiseSize           = 1;

LB = -2; UB = 2; PLB = -0.2; PUB = 0.2;
LB_all = [0 -0.2 0.01]; UB_all = [10 0.2 10];

for ci = 1: Ncond
    binz_all_condi = [];
    prop_cf_high_all_condi = [];
    rt_all_condi = [];
    
    for si = 1: Nsubj
        binz_all_condi = [binz_all_condi; squeeze(binz_posE(si,ci,:))];
        prop_cf_high_all_condi = [prop_cf_high_all_condi; 1-squeeze(prop_cf_high_all(si,ci,:))];
        rt_si_ci = squeeze(rt_all(si,ci,:));
        rt_si_ci = (rt_si_ci - nanmean(rt_si_ci))./nanstd(rt_si_ci);
        rt_all_condi = [rt_all_condi; rt_si_ci];
        
        indi_cf = find(squeeze(prop_cf_high_all(si,ci,:))== min(squeeze(prop_cf_high_all(si,ci,:))));
        if length(indi_cf) == 1
            mu_estim_cf(si,ci) = binz_posE(si,ci,indi_cf(1));
        else
            mu_estim_cf(si,ci)= mean(binz_posE(si,ci, indi_cf));
        end
        
        indi_rt = find(squeeze(rt_all(si,ci,:))== max(squeeze(rt_all(si,ci,:))));
        mu_estim_rt(si,ci) = binz_posE(si,ci,indi_rt(1));
        
    end
    
end

mu_cf_fit_all = mu_estim_cf;
mu_rt_fit_all = mu_estim_rt;

delta_mu_estim_cf_MAE = mu_estim_cf(:,3)-mu_estim_cf(:,1);
delta_mu_estim_rt_MAE = mu_estim_rt(:,3)-mu_estim_rt(:,1);
MAE_strength = -(squeeze(params_psych_m2_all(:,3,1))-squeeze(params_psych_m2_all(:,1,1)));

[r_cf_MAE, p_cf_MAE]= corr(delta_mu_estim_cf_MAE, MAE_strength, 'type', 'Spearman')
[r_rt_MAE, p_rt_MAE]= corr(delta_mu_estim_rt_MAE, MAE_strength, 'type', 'Spearman')
delta_mu_estim_cf_compensation = mu_estim_cf(:,4)-mu_estim_cf(:,3);
delta_mu_estim_rt_compensation = mu_estim_rt(:,4)-mu_estim_rt(:,3);
[r_cf_compensation, p_cf_compensation]= corr(delta_mu_estim_cf_compensation,compensation, 'type', 'Spearman')
[r_rt_compensation, p_rt_compensation]= corr(delta_mu_estim_rt_compensation, compensation, 'type', 'Spearman')
%% Correlations from mus
%close all;
figure(12)
set(gcf, 'Position', [100 100 500 340])
marginsa3 = [0.1200    0.100    0.1400    0.090]; %MARGINS = [LEFT RIGHT BOTTOM TOP]
guttera3 = [ 0.1 0.16];
if exp_i == 1
    limis = [- 0.5 0.6];
else
    limis = [-0.2 0.3];
end
msz = 3;

for ci = 3: Ncond
    tight_subplot(2,2,1,ci-2, guttera3, marginsa3)
    
    plot(squeeze(params_psych_m2_all(indi_sel,ci, 1)), squeeze(mu_cf_fit_all(indi_sel,ci)),'o', 'MarkerSize', msz, 'Color', colorz(ci,:), 'MarkerFaceColor',colorz(ci,:) ); hold on;
    axis equal
    lsline
    box off
    set(gca, 'tickdir', 'out')
    [r_cf(ci), p_cf(ci)]= corr(squeeze(params_psych_m2_all(indi_sel,ci, 1)), squeeze(mu_cf_fit_all(indi_sel,ci)), 'type', 'Spearman');
    title(['\rho =', num2str(r_cf(ci), '%.2f'), ', p = ', num2str(p_cf(ci), '%.2f')]); hold on;
    xlim(limis)
    ylim(limis)
    set(gca, 'FontSize', fontsz)
    set(gca, 'ticklength',2*[tlen1 tlen2])
    if ci == 3
        ylabel(' \mu from the confidence curves', 'FontName', 'Helvetica','FontSize', fontsz)
    end
    
    tight_subplot(2,2,2,ci-2, guttera3, marginsa3)
    plot(squeeze(params_psych_m2_all(indi_sel,ci, 1)), squeeze(mu_rt_fit_all(indi_sel,ci)),'o','MarkerSize', msz, 'Color',colorz(ci,:), 'MarkerFaceColor',colorz(ci,:) ); hold on;
    axis equal
    lsline
    box off
    set(gca, 'tickdir', 'out')
    [r_rt(ci), p_rt(ci)]= corr(squeeze(params_psych_m2_all(indi_sel,ci, 1)), squeeze(mu_rt_fit_all(indi_sel,ci)), 'type', 'Spearman');
    title(['\rho =', num2str(r_rt(ci), '%.2f'), ', p = ', num2str(p_rt(ci), '%.2f')]); hold on;
    xlabel(' \mu from the psychometric curve', 'FontName', 'Helvetica','FontSize', fontsz)
    xlim(limis)
    ylim(limis)
    set(gca, 'FontSize', fontsz)
    set(gca, 'ticklength',2*[tlen1 tlen2])
    if ci == 3
        ylabel(' \mu from the RT curves', 'FontName', 'Helvetica','FontSize', fontsz)
    end
end
%%
psname = ['Corrs_btwn_mu_based_on_peaks_averaged_exp_', num2str(exp_i), 'CTR.pdf']
%print_pdf(psname)
%%
% mu estim from the 3 curves
mus_psych = squeeze(params_psych_m2_all(:,:, 1));
% cond 3:
[p,h] = signrank(mus_psych(:,3),  mu_cf_fit_all(:,3)); % 0.49
[p,h] = signrank(mus_psych(:,3),  mu_rt_fit_all(:,3)); % 0.98

% cond 4:
[p,h] = signrank(mus_psych(:,4),  mu_cf_fit_all(:,4)); % 0.306
[p,h] = signrank(mus_psych(:,4),  mu_rt_fit_all(:,4)); % 0.92

for ci = 1:4
    [p_mus_ci_psych_cf(ci),h] = signrank(mus_psych(:,ci),  mu_cf_fit_all(:,ci));%signrank(mus_psych(:,ci),  mu_cf_fit_all(:,ci));
    [p_mus_ci_psych_rt(ci),h] = signrank(mus_psych(:,ci),  mu_rt_fit_all(:,ci));%signrank(mus_psych(:,ci),  mu_rt_fit_all(:,ci));
end


%% confidence and RT track with psychometric curves: confidence curves paper
%close all;
%hf2 = figure(2)
clear h;
figure(4)
set(gcf, 'Position', [100 100 600 800])

linewi = 1;%1.1;%0.9;
sz_dot = 2;
dashed_linez = 0;

%indi_sel = 1:1:Nsubj;%[5 6 7 8 9 10];%[1 2 5 6  9 10];%1: 1: Nsubj;% [1 2 4 5 8 9 10 ];%1: 1: Nsubj;
% AGGREGATE DATA FROM EXP 1

nboot = 5000;

if exp_i == 1
    tight_subplot(4,3,1,1, guttera, marginsa)
    for  ci = 1:Ncond
        
        h(ci) = plot(mean(squeeze(binz_pos(indi_sel,ci,:)),1), nanmean(squeeze(prop_cw_all(indi_sel,ci,:)),1),'o','MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:)); hold on;
        errorbar(mean(squeeze(binz_pos(indi_sel,ci,:)),1), nanmean(squeeze(prop_cw_all(indi_sel,ci,:)),1), nanstd(squeeze(prop_cw_all(indi_sel,ci,:)),1)/sqrt(length(indi_sel)), 'Color',  colorz(ci,:), 'Linewidth', linewi, 'CapSize',0); hold on;
        
        if dashed_linez
            
            plot(nanmean(params_psych_m2_all(indi_sel,ci,1)) *ones(1,10), linspace(0,1,10),'--', 'Color', colorz(ci,:),'Linewidth', linewi ); hold on;
            plot(linspace(-0.5,nanmean(params_psych_m2_all(indi_sel,ci,1)),10), 0.5*ones(1,10),'--', 'Color', colorz(ci,:),'Linewidth',linewi ); hold on;
            
        end
    end
    box off
    
    if exp_i == 1
        xlim_min = -0.5; xlim_max = 0.5;
    elseif exp_i == 2
        xlim_min = -0.3; xlim_max = 0.3;
    end
    ylim_min = 0;
    ylim_max = 1.01;
    xlim([xlim_min xlim_max])
    ylim([ylim_min ylim_max])
    plot(zeros(1,nbinz), linspace(0,1,nbinz), '--k'); hold on;
    box off
    set(gca, 'xtick',[-0.5:0.1:0.5])
    set(gca, 'ytick', [0:0.25:1])
    set(gca, 'tickdir', 'out')
    set(gca, 'FontSize', fontsz)
    set(gca, 'ticklength',[tlen1 tlen2])
    %xlabel('test stimulus speed clockwise (a.u.)') % degrees of visual angle (dva)
    ylabel('proportion response clockwise')
    title('Exp 1: group')
    
    
    tight_subplot(4,3,2,1,guttera, marginsa)
    for  ci = 1:Ncond
        
        h(ci) = plot(mean(squeeze(binz_pos(indi_sel,ci,:)),1), nanmean(squeeze(prop_cf_high_all(indi_sel,ci,:)),1),'rs','MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:)); hold on;
        errorbar(mean(squeeze(binz_pos(indi_sel,ci,:)),1), nanmean(squeeze(prop_cf_high_all(indi_sel,ci,:)),1), nanstd(squeeze(prop_cf_high_all(indi_sel,ci,:)),1)/sqrt(length(indi_sel)), 'Color',  colorz(ci,:), 'Linewidth', linewi, 'CapSize',0); hold on;
        if dashed_linez
            plot(nanmean(params_psych_m2_all(indi_sel,ci,1)) *ones(1,10), linspace(0,1,10),'--', 'Color', colorz(ci,:),'Linewidth',linewi ); hold on;
        end
    end
    box off
    xlim([xlim_min xlim_max])
    ylim([ylim_min ylim_max])
    plot(zeros(1,nbinz), linspace(0,1,nbinz), '--k'); hold on;
    box off
    set(gca, 'xtick',[-0.5:0.1:0.5])
    set(gca, 'ytick', [0:0.25:1])
    set(gca, 'tickdir', 'out')
    set(gca, 'FontSize', fontsz)
    set(gca, 'ticklength',[tlen1 tlen2])
    % xlabel('test stimulus speed clockwise (a.u.)') % degrees of visual angle (dva)
    ylabel('proportion high confidence')
    
    
    
    tight_subplot(4,3,3,1,guttera, marginsa)
    for  ci = 1:Ncond
        
        h(ci) = plot(mean(squeeze(binz_pos(indi_sel,ci,:)),1), nanmean(squeeze(rt_all(indi_sel,ci,:)),1),'^','MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:)); hold on;
        errorbar(mean(squeeze(binz_pos(indi_sel,ci,:)),1), nanmean(squeeze(rt_all(indi_sel,ci,:)),1), nanstd(squeeze(rt_all(indi_sel,ci,:)),1)/sqrt(length(indi_sel)), 'Color',  colorz(ci,:), 'Linewidth', linewi, 'CapSize',0); hold on;
        if dashed_linez
            plot(nanmean(params_psych_m2_all(indi_sel,ci,1)) *ones(1,10), linspace(0,2.25*ylim_max,10),'--', 'Color', colorz(ci,:),'Linewidth',linewi ); hold on;
        end
    end
    box off
    
    xlim([xlim_min xlim_max])
    ylim([1 2.25*ylim_max])
    plot(zeros(1,nbinz), linspace(0,3,nbinz), '--k'); hold on;
    box off
    set(gca, 'xtick',[-0.5:0.1:0.5])
    set(gca, 'tickdir', 'out')
    set(gca, 'FontSize', fontsz)
    set(gca, 'ticklength',[tlen1 tlen2])
    xlabel('test stimulus speed clockwise (a.u.)')
    ylabel('reaction times (s)')
    
    
    tight_subplot(4,3,4,1,guttera, marginsa)
    
    for ci = 1: Ncond
        pgr_psych(ci) = plot(nanmean(params_psych_m2_all(indi_sel,ci,1)) , 20-5*(ci-1), 'o', 'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor',colorz(ci,:)); hold on;
        pgr_cf(ci) = plot(nanmean(mu_estim_cf(indi_sel,ci,1)), 19-5*(ci-1),'rs','MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor',colorz(ci,:)); hold on;
        pgr_rt(ci) = plot(nanmean(mu_estim_rt(indi_sel,ci,1)), 18-5*(ci-1), '^','MarkerFaceColor', colorz(ci,:),'MarkerEdgeColor',colorz(ci,:)); hold on;
        
        errorbar(nanmean(params_psych_m2_all(indi_sel,ci,1)) , 20-5*(ci-1), nanstd(params_psych_m2_all(indi_sel,ci,1)), 'horizontal','Color', colorz(ci,:), 'Capsize', 0,'Linewidth',linewi); hold on;
        errorbar(nanmean(mu_estim_cf(indi_sel,ci,1)) , 19-5*(ci-1), nanstd(mu_estim_cf(indi_sel,ci)), 'horizontal','Color', colorz(ci,:),'Capsize', 0,'Linewidth',linewi); hold on;
        errorbar(nanmean(mu_estim_rt(indi_sel,ci,1)) , 18-5*(ci-1), nanstd(mu_estim_rt(indi_sel,ci)), 'horizontal','Color', colorz(ci,:),'Capsize', 0,'Linewidth',linewi); hold on;
    end
    plot(zeros(1,nbinz), linspace(0,22,nbinz), '--k'); hold on;
    ylim([0 22])
    box off
    set(gca, 'tickdir', 'out')
    set(gca, 'FontSize', fontsz)
    set(gca, 'ticklength',[tlen1 tlen2])
    
    
elseif exp_i == 2
    % DECIDE ON A PARTICIPANT
    %sii = 13;
    %in the control experiment, only 8 participants total
    sii = 8;
    tight_subplot(4,3,1,2, guttera, marginsa)
    
    for  ci = 1:Ncond
        h(ci) = plot(squeeze(binz_pos(sii,ci,:)), squeeze(prop_cw_all(sii,ci,:)),'o','MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:)); hold on;
        h(ci) = plot(squeeze(binz_pos(sii,ci,:)), squeeze(prop_cw_pr_all(sii,ci,:)),'Color', colorz(ci,:), 'Linewidth', linewi); hold on;
        if dashed_linez
            if ci>=3
                plot(params_psych_m2_all(sii,ci,1) *ones(1,10), linspace(0,1,10),'--', 'Color', colorz(ci,:),'Linewidth',linewi ); hold on;
                plot(linspace(-0.3,params_psych_m2_all(sii,ci,1),10), 0.5*ones(1,10),'--', 'Color', colorz(ci,:),'Linewidth',linewi ); hold on;
            end
        end
    end
    if exp_i == 1
        xlim_min = -0.5; xlim_max = 0.5;
    elseif exp_i == 2
        xlim_min = -0.3; xlim_max = 0.3;
    end
    ylim_min = 0;
    ylim_max = 1.01;
    xlim([xlim_min xlim_max])
    ylim([ylim_min ylim_max])
    plot(zeros(1,nbinz), linspace(0,1,nbinz), '--k'); hold on;
    box off
    set(gca, 'xtick',[-0.3:0.1:0.3])
    set(gca, 'ytick', [0:0.25:1])
    set(gca, 'tickdir', 'out')
    set(gca, 'FontSize', fontsz)
    set(gca, 'ticklength',[tlen1 tlen2])
    %xlabel('test stimulus speed clockwise (a.u.)') % degrees of visual angle (dva)
    ylabel('prop response clockwise')
    title('Exp 2 ctr: participant 8')
    
    
    tight_subplot(4,3,2,2,guttera, marginsa)
    for  ci = 1:Ncond
        h(ci) = plot(squeeze(binz_pos(sii,ci,:)), squeeze(squeeze(prop_cf_high_all(sii,ci,:))),'rs','MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:)); hold on;
        hp(ci) = plot(squeeze(binz_pos(sii,ci,:)), squeeze(squeeze(prop_cf_high_all(sii,ci,:))),'Color', colorz(ci,:), 'Linewidth', linewi); hold on;
        if dashed_linez
            plot(squeeze(params_psych_m2_all(sii,ci,1)) *ones(1,10), linspace(0,1,10),'--', 'Color', colorz(ci,:),'Linewidth',linewi ); hold on;
        end
    end
    box off
    xlim([xlim_min xlim_max])
    ylim([ylim_min ylim_max])
    plot(zeros(1,nbinz), linspace(0,1,nbinz), '--k'); hold on;
    box off
    set(gca, 'xtick',[-0.3:0.1:0.3])
    set(gca, 'ytick', [0:0.25:1])
    set(gca, 'tickdir', 'out')
    set(gca, 'FontSize', fontsz)
    set(gca, 'ticklength',[tlen1 tlen2])
    %xlabel('test stimulus speed clockwise (a.u.)') % degrees of visual angle (dva)
    ylabel('proportion high confidence')
    
    
    
    tight_subplot(4,3,3,2,guttera, marginsa)
    for  ci = 1:Ncond
        
        h(ci) = plot(squeeze(binz_pos(sii,ci,:)), squeeze(squeeze(rt_all(sii,ci,:))),'^','MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:)); hold on;
        hp(ci) = plot(squeeze(binz_pos(sii,ci,:)), squeeze(squeeze(rt_all(sii,ci,:))),'Color', colorz(ci,:), 'Linewidth', linewi); hold on;
        if dashed_linez
            plot(squeeze(params_psych_m2_all(sii,ci,1))*ones(1,10), linspace(0,2.25*ylim_max,10),'--', 'Color', colorz(ci,:),'Linewidth',linewi ); hold on;
        end
        
    end
    box off
    xlim([xlim_min xlim_max])
    ylim([1 2.25*ylim_max])
    plot(zeros(1,nbinz), linspace(0,3,nbinz), '--k'); hold on;
    box off
    set(gca, 'xtick',[-0.3:0.1:0.3])
    set(gca, 'tickdir', 'out')
    set(gca, 'FontSize', fontsz)
    set(gca, 'ticklength',[tlen1 tlen2])
    xlabel('test stimulus speed clockwise (a.u.)')
    ylabel('reaction times (s)')
    
    
    % AGGREGATE DATA FROM EXP 2
    
    tight_subplot(4,3,1,3, guttera, marginsa)
    for  ci = 1:Ncond
        
        h(ci) = plot(mean(squeeze(binz_pos(indi_sel,ci,:)),1), nanmean(squeeze(prop_cw_all(indi_sel,ci,:)),1),'o','MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:) ); hold on;
        errorbar(mean(squeeze(binz_pos(indi_sel,ci,:)),1), nanmean(squeeze(prop_cw_all(indi_sel,ci,:)),1), nanstd(squeeze(prop_cw_all(indi_sel,ci,:)),1)/sqrt(length(indi_sel)), 'Color',  colorz(ci,:), 'Linewidth', linewi, 'CapSize',0); hold on;
        % nan is needed to accomodate for exp 1
        if dashed_linez
            plot(nanmean(params_psych_m2_all(indi_sel,ci,1)) *ones(1,10), linspace(0,1,10),'--', 'Color', colorz(ci,:),'Linewidth',linewi ); hold on;
            plot(linspace(-0.3,nanmean(params_psych_m2_all(indi_sel,ci,1)),10), 0.5*ones(1,10),'--', 'Color', colorz(ci,:),'Linewidth',linewi ); hold on;
        end
    end
    box off
    
    
    if exp_i == 1
        xlim_min = -0.5; xlim_max = 0.5;
    elseif exp_i == 2
        xlim_min = -0.3; xlim_max = 0.3;
    end
    ylim_min = 0;
    ylim_max = 1.01;
    xlim([xlim_min xlim_max])
    ylim([ylim_min ylim_max])
    plot(zeros(1,nbinz), linspace(0,1,nbinz), '--k'); hold on;
    box off
    set(gca, 'xtick',[-0.3:0.1:0.3])
    set(gca, 'ytick', [0:0.25:1])
    set(gca, 'tickdir', 'out')
    set(gca, 'FontSize', fontsz)
    set(gca, 'ticklength',[tlen1 tlen2])
    %xlabel('test stimulus speed clockwise (a.u.)') % degrees of visual angle (dva)
    %ylabel('proportion response clockwise')
    title('Exp 2 ctr: group (N = 8)')
    
    
    
    tight_subplot(4,3,2,3,guttera, marginsa)
    for  ci = 1:Ncond
        
        h(ci) = plot(mean(squeeze(binz_pos(indi_sel,ci,:)),1), nanmean(squeeze(prop_cf_high_all(indi_sel,ci,:)),1),'rs','MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:) ); hold on;
        errorbar(mean(squeeze(binz_pos(indi_sel,ci,:)),1), nanmean(squeeze(prop_cf_high_all(indi_sel,ci,:)),1), nanstd(squeeze(prop_cf_high_all(indi_sel,ci,:)),1)/sqrt(length(indi_sel)), 'Color',  colorz(ci,:), 'Linewidth', linewi, 'CapSize',0); hold on;
        if dashed_linez
            plot(nanmean(params_psych_m2_all(indi_sel,ci,1)) *ones(1,10), linspace(0,1,10),'--', 'Color', colorz(ci,:),'Linewidth',linewi ); hold on;
        end
    end
    box off
    xlim([xlim_min xlim_max])
    ylim([ylim_min ylim_max])
    plot(zeros(1,nbinz), linspace(0,3,nbinz), '--k'); hold on;
    box off
    set(gca, 'xtick',[-0.3:0.1:0.3])
    set(gca, 'ytick', [0:0.25:1])
    set(gca, 'tickdir', 'out')
    set(gca, 'FontSize', fontsz)
    set(gca, 'ticklength',[tlen1 tlen2])
    %xlabel('test stimulus speed clockwise (a.u.)') % degrees of visual angle (dva)
    %ylabel('proportion high confidence')
    
    
    tight_subplot(4,3,3,3,guttera, marginsa)
    for  ci = 1:Ncond
        h(ci) = plot(mean(squeeze(binz_pos(indi_sel,ci,:)),1), nanmean(squeeze(rt_all(indi_sel,ci,:)),1),'^','MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:) ); hold on;
        errorbar(mean(squeeze(binz_pos(indi_sel,ci,:)),1), nanmean(squeeze(rt_all(indi_sel,ci,:)),1), nanstd(squeeze(rt_all(indi_sel,ci,:)),1)/sqrt(length(indi_sel)), 'Color',  colorz(ci,:), 'Linewidth', linewi, 'CapSize',0); hold on;
        if dashed_linez
            plot(nanmean(params_psych_m2_all(indi_sel,ci,1)) *ones(1,10), linspace(0,2.25*ylim_max,10),'--', 'Color', colorz(ci,:),'Linewidth',linewi ); hold on;
        end
        
    end
    box off
    plot(zeros(1,nbinz), linspace(0,3,nbinz), '--k'); hold on;
    xlim([xlim_min xlim_max])
    ylim([1 2.25*ylim_max])
    box off
    set(gca, 'xtick',[-0.3:0.1:0.3])
    set(gca, 'tickdir', 'out')
    set(gca, 'FontSize', fontsz)
    set(gca, 'ticklength',[tlen1 tlen2])
    xlabel('test stimulus speed clockwise (a.u.)')
    %ylabel('reaction times (s)')
    
    
    tight_subplot(4,3,4,3,guttera, marginsa)
    for ci = 1: Ncond
        
        vala_psych = squeeze(params_psych_m2_all(indi_sel,ci,1));
        vala_cf = squeeze(mu_estim_cf(indi_sel,ci,1));
        vala_rt = squeeze(mu_estim_rt(indi_sel,ci,1));
        clear psych_samples; clear cf_samples; clear rt_samples;
        
        for kk = 1:nboot
            sample_psych = randsample(vala_psych,length(indi_sel),1);
            psych_samples(kk) = nanmedian(sample_psych);
            
            sample_cf = randsample(vala_cf,length(indi_sel),1);
            cf_samples(kk) = nanmedian(sample_cf);
            
            sample_rt = randsample(vala_rt,length(indi_sel),1);
            rt_samples(kk) = nanmedian(sample_rt);
            
        end
        bci_psych = [quantile(psych_samples,ci_bnd_low); quantile(psych_samples,ci_bnd_high)];
        bci_cf = [quantile(cf_samples,ci_bnd_low); quantile(cf_samples,ci_bnd_high)];
        bci_rt = [quantile(rt_samples,ci_bnd_low); quantile(rt_samples,ci_bnd_high)];
        
        pgr_psych(ci) = plot(nanmedian(vala_psych) , 20-5*(ci-1), 'o', 'MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor',colorz(ci,:)); hold on;
        pgr_cf(ci) = plot(nanmedian(vala_cf), 19-5*(ci-1),'rs','MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor',colorz(ci,:)); hold on;
        pgr_rt(ci) = plot(nanmedian(vala_rt), 18-5*(ci-1), '^','MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:),'MarkerEdgeColor',colorz(ci,:)); hold on;
        
        errorbar(nanmedian(vala_psych) , 20-5*(ci-1),bci_psych(1),bci_psych(2), 'horizontal','Color', colorz(ci,:), 'Capsize', 0,'Linewidth',linewi); hold on;
        errorbar(nanmedian(vala_cf) , 19-5*(ci-1), bci_cf(1),bci_cf(2), 'horizontal','Color', colorz(ci,:),'Capsize', 0,'Linewidth',linewi); hold on;
        errorbar(nanmedian(vala_rt) , 18-5*(ci-1), bci_rt(1),bci_rt(2), 'horizontal','Color', colorz(ci,:),'Capsize', 0,'Linewidth',linewi); hold on;
    end
    plot(zeros(1,nbinz), linspace(0,22,nbinz), '--k'); hold on;
    ylim([0 22])
    xlim([-0.3 0.3])
    box off
    set(gca, 'tickdir', 'out')
    set(gca, 'FontSize', fontsz)
    set(gca, 'ticklength',[tlen1 tlen2])
    l_pgr = legend([pgr_psych(1) pgr_cf(1) pgr_rt(1)], {'Psychometric', 'Confidence', 'RT'}, 'FontName', 'Helvetica', 'FontSize', fontsz)
    set(l_pgr, 'Position', [0.4309    0.1779    0.1908    0.0634])
    
    
end

%%
psname = ['Fig4_Psych_confidence_RT_Means_exp_',num2str(exp_i), 'CTR.pdf']
%print_pdf(psname)
%%

hf2 = figure(13)
set(gcf, 'Position', [100 100 600 700])

linewi = 1.1;
dashed_linez = 0;

if exp_i == 2
    si_vec = [ 1 3 6]; % pick 3 participants from exp 2 ctr
    for  sii_ind = 1:3
        sii = si_vec(sii_ind);
        tight_subplot(3,3,1,sii_ind, guttera, marginsa)
        
        for  ci = 1:Ncond
            h(ci) = plot(squeeze(binz_pos(sii,ci,:)), squeeze(prop_cw_all(sii,ci,:)),'o','MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:)); hold on;
            h(ci) = plot(squeeze(binz_pos(sii,ci,:)), squeeze(prop_cw_pr_all(sii,ci,:)),'Color', colorz(ci,:), 'Linewidth', linewi); hold on;
            if dashed_linez
                if ci>=3
                    plot(params_psych_m2_all(sii,ci,1) *ones(1,10), linspace(0,1,10),'--', 'Color', colorz(ci,:),'Linewidth',linewi ); hold on;
                    plot(linspace(-0.3,params_psych_m2_all(sii,ci,1),10), 0.5*ones(1,10),'--', 'Color', colorz(ci,:),'Linewidth',linewi ); hold on;
                end
            end
        end
        
        if exp_i == 1
            xlim_min = -0.5; xlim_max = 0.5;
        elseif exp_i == 2
            xlim_min = -0.3; xlim_max = 0.3;
        end
        ylim_min = 0;
        ylim_max = 1.01;
        xlim([xlim_min xlim_max])
        ylim([ylim_min ylim_max])
        plot(zeros(1,nbinz), linspace(0,1,nbinz), '--k'); hold on;
        box off
        set(gca, 'xtick',[-0.3:0.1:0.3])
        set(gca, 'ytick', [0:0.25:1])
        set(gca, 'tickdir', 'out')
        set(gca, 'FontSize', fontsz)
        set(gca, 'ticklength',[tlen1 tlen2])
        %xlabel('test stimulus speed clockwise (a.u.)') % degrees of visual angle (dva)
        if sii_ind == 1
            ylabel('proportion response clockwise')
        end
        title(['Exp 2 ctr: participant ', num2str(sii)])
        
        
        tight_subplot(3,3,2,sii_ind,guttera, marginsa)
        for  ci = 1:Ncond
            h(ci) = plot(squeeze(binz_pos(sii,ci,:)), squeeze(squeeze(prop_cf_high_all(sii,ci,:))),'o','MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:)); hold on;
            hp(ci) = plot(squeeze(binz_pos(sii,ci,:)), squeeze(squeeze(prop_cf_high_all(sii,ci,:))),'Color', colorz(ci,:), 'Linewidth', linewi); hold on;
            if dashed_linez
                plot(squeeze(params_psych_m2_all(sii,ci,1)) *ones(1,10), linspace(0,1,10),'--', 'Color', colorz(ci,:),'Linewidth',linewi ); hold on;
            end
        end
        box off
        xlim([xlim_min xlim_max])
        ylim([ylim_min ylim_max])
        plot(zeros(1,nbinz), linspace(0,1,nbinz), '--k'); hold on;
        box off
        set(gca, 'xtick',[-0.3:0.1:0.3])
        set(gca, 'ytick', [0:0.25:1])
        set(gca, 'tickdir', 'out')
        set(gca, 'FontSize', fontsz)
        set(gca, 'ticklength',[tlen1 tlen2])
        %xlabel('test stimulus speed clockwise (a.u.)') % degrees of visual angle (dva)
        if sii_ind == 1
            ylabel('proportion high confidence')
        end
        
        
        tight_subplot(3,3,3,sii_ind,guttera, marginsa)
        for  ci = 1:Ncond
            
            h(ci) = plot(squeeze(binz_pos(sii,ci,:)), squeeze(squeeze(rt_all(sii,ci,:))),'o','MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:)); hold on;
            hp(ci) = plot(squeeze(binz_pos(sii,ci,:)), squeeze(squeeze(rt_all(sii,ci,:))),'Color', colorz(ci,:), 'Linewidth', linewi); hold on;
            if dashed_linez
                plot(squeeze(params_psych_m2_all(sii,ci,1)) *ones(1,10), linspace(0,2.25*ylim_max,10),'--', 'Color', colorz(ci,:),'Linewidth',linewi ); hold on;
            end
        end
        box off
        xlim([xlim_min xlim_max])
        ylim([1 2.25*ylim_max])
        plot(zeros(1,nbinz), linspace(0,3,nbinz), '--k'); hold on;
        box off
        set(gca, 'xtick',[-0.3:0.1:0.3])
        set(gca, 'tickdir', 'out')
        set(gca, 'FontSize', fontsz)
        set(gca, 'ticklength',[tlen1 tlen2])
        if sii_ind == 2
            xlabel('test stimulus speed clockwise (a.u.)')
        end
        if sii_ind == 1
            ylabel('reaction times (s)')
        end
    end
    
end
%%
psname = 'Fig3_CTR_Supp_Psych_confidence_RT_paper.pdf'
print_pdf(psname)

%% MODEL FIGURE
%close all;

figure(5)
set(gcf, 'Position', [100 100 740 700])

indi_sel = 1:1:Nsubj;
xlim_min = -0.3;%min(min(min(binz_pos)));
xlim_max = 0.3;%max(max(max(binz_pos)));


tight_subplot(5,4,3,1,guttera, marginsa)
for  ci = 1:Ncond
    
    h(ci) = plot(mean(squeeze(binz_pos(indi_sel,ci,:)),1), nanmean(squeeze(rt_all(indi_sel,ci,:)),1),'o','MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:)); hold on;
    errorbar(mean(squeeze(binz_pos(indi_sel,ci,:)),1), nanmean(squeeze(rt_all(indi_sel,ci,:)),1), nanstd(squeeze(rt_all(indi_sel,ci,:)),1)/sqrt(length(indi_sel)), 'Color',  colorz(ci,:), 'LineStyle', 'none','CapSize',0, 'Linewidth', linewi); hold on;
    
end
box off
xlim([xlim_min xlim_max])
ylim([1 2])
plot(zeros(1,nbinz), linspace(0,4,nbinz), '--k'); hold on;
box off
set(gca, 'xtick',[-0.3:0.1:0.3])
set(gca, 'tickdir', 'out')
set(gca, 'FontSize', fontsz)
set(gca, 'ticklength',[tlen1 tlen2])
ylabel('reaction times (s)')



% decision variable
if exp_i == 2
    tight_subplot(5,4,3,2,guttera, marginsa)
    for  ci = 1:Ncond
        
        % decision variable
        if exp_i == 2
            %{
            tight_subplot(5,4,3,2,guttera, marginsa)
            for  ci = 1:Ncond
                h(ci) = plot(mean(squeeze(binz_pos(indi_sel,ci,:)),1), -nanmean(log(abs(squeeze(dv_val_binz(indi_sel,ci,:)))),1),'o','MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:)); hold on;
                errorbar(mean(squeeze(binz_pos(indi_sel,ci,:)),1), -nanmean(log(abs(squeeze(dv_val_binz(indi_sel,ci,:)))),1), nanstd(squeeze(log(abs(dv_val_binz(indi_sel,ci,:)))),1)/sqrt(length(indi_sel)), 'Color',  colorz(ci,:),'CapSize',0, 'Linewidth', linewi); hold on;
            end
            box off
            xlim([xlim_min xlim_max])
            plot(zeros(1,nbinz), linspace(-5,1,nbinz), '--k'); hold on;
            box off
            set(gca, 'xtick',[-0.3:0.1:0.3])
            set(gca, 'tickdir', 'out')
            set(gca, 'FontSize', fontsz)
            set(gca, 'ticklength',[tlen1 tlen2])
            set(gca, 'tickdir', 'out')
            %xlabel('test stim')
            %ylabel('-decision variable')
            ylabel('-log (|d|)')
            xlabel('test stimulus speed clockwise (a.u.)')
            %}
        end
        
        msz_all = [5.5 5.5 5.5 5.5];
        
        % Bayesian model 1
        tight_subplot(5,4,1,3, guttera, marginsa)
        mi = 1;
        for  ci = 1:Ncond
            if exp_i == 2
                he(ci) = fill([mean(squeeze(binz_pos(indi_sel,ci,:)),1) mean(squeeze(binz_pos(indi_sel,ci,end:-1:1)),1)],   ...
                    [mean(squeeze(prop_cw_pred_all(indi_sel,mi,ci,:)),1)-std(squeeze(prop_cw_pred_all(indi_sel,mi,ci,:)),1)./sqrt(length(indi_sel))...
                    fliplr(mean(squeeze(prop_cw_pred_all(indi_sel,mi,ci,:)),1) + std(squeeze(prop_cw_pred_all(indi_sel,mi,ci,:)),1)./sqrt(length(indi_sel))) ],colorz_shade(ci,:), 'EdgeColor', 'None'); hold on;
            end
            plot(nanmean(squeeze(binz_pos(:,ci,:)),1), nanmean(squeeze(prop_cw_all(indi_sel,ci,:)),1),'o','MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:) ); hold on;
            errorbar(nanmean(squeeze(binz_pos(:,ci,:)),1), nanmean(squeeze(prop_cw_all(indi_sel,ci,:)),1), nanstd(squeeze(prop_cw_all(indi_sel,ci,:)),1)/sqrt(length(indi_sel)), 'Color',  colorz(ci,:), 'LineStyle', 'none','CapSize',0,'Linewidth',linewi); hold on;
        end
        box off
        xlim_min = -0.3;%min(min(min(binz_pos)));
        xlim_max = 0.3;%max(max(max(binz_pos)));
        ylim_min = 0;
        ylim_max = 1.01;
        xlim([xlim_min xlim_max])
        ylim([ylim_min ylim_max])
        plot(zeros(1,nbinz), linspace(0,1,nbinz), '--k'); hold on;
        box off
        set(gca, 'xtick', [-0.3:0.1:0.3])
        set(gca, 'xtick',[-0.3:0.1:0.3])
        set(gca, 'ytick', [0:0.25:1])
        set(gca, 'tickdir', 'out')
        set(gca, 'FontSize', fontsz)
        set(gca, 'ticklength',[tlen1 tlen2])
        set(gca, 'tickdir', 'out')
        %xlabel('test stim')
        ylabel('prop resp CW')
        %h(ci) = plot(mean(squeeze(binz_pos(indi_sel,ci,:)),1), -nanmean(log(abs(squeeze(dv_val_binz(indi_sel,ci,:)))),1),'o','MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:)); hold on;
        %errorbar(mean(squeeze(binz_pos(indi_sel,ci,:)),1), -nanmean(log(abs(squeeze(dv_val_binz(indi_sel,ci,:)))),1), nanstd(squeeze(log(abs(dv_val_binz(indi_sel,ci,:)))),1)/sqrt(length(indi_sel)), 'Color',  colorz(ci,:),'CapSize',0, 'Linewidth', linewi); hold on;
    end
    
end

msz_all = [5.5 5.5 5.5 5.5];

% Bayesian model 1
tight_subplot(5,4,1,3, guttera, marginsa)
%mi = 1;
mi = 5;
for  ci = 1:Ncond
    if exp_i == 2
        he(ci) = fill([mean(squeeze(binz_pos(indi_sel,ci,:)),1) mean(squeeze(binz_pos(indi_sel,ci,end:-1:1)),1)],   ...
            [mean(squeeze(prop_cw_pred_all(indi_sel,mi,ci,:)),1)-std(squeeze(prop_cw_pred_all(indi_sel,mi,ci,:)),1)./sqrt(length(indi_sel))...
            fliplr(mean(squeeze(prop_cw_pred_all(indi_sel,mi,ci,:)),1) + std(squeeze(prop_cw_pred_all(indi_sel,mi,ci,:)),1)./sqrt(length(indi_sel))) ],colorz_shade(ci,:), 'EdgeColor', 'None'); hold on;
    end
    plot(nanmean(squeeze(binz_pos(:,ci,:)),1), nanmean(squeeze(prop_cw_all(indi_sel,ci,:)),1),'o', 'Color',colorz(ci,:),'Linewidth', 0.1,'MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:),'Linewidth',0.01 ); hold on;
    errorbar(nanmean(squeeze(binz_pos(:,ci,:)),1), nanmean(squeeze(prop_cw_all(indi_sel,ci,:)),1), nanstd(squeeze(prop_cw_all(indi_sel,ci,:)),1)/sqrt(length(indi_sel)), 'Color',  colorz(ci,:), 'LineStyle', 'none','CapSize',0, 'Linewidth', linewi); hold on;
    
end
box off
xlim_min = -0.3;%min(min(min(binz_pos)));
xlim_max = 0.3;%max(max(max(binz_pos)));
ylim_min = 0;
ylim_max = 1.01;
xlim([xlim_min xlim_max])
ylim([ylim_min ylim_max])
plot(zeros(1,nbinz), linspace(0,1,nbinz), '--k'); hold on;
box off
set(gca, 'xtick', [-0.3:0.1:0.3])
set(gca, 'xtick',[-0.3:0.1:0.3])
set(gca, 'ytick', [0:0.25:1])
set(gca, 'tickdir', 'out')
set(gca, 'FontSize', fontsz)
set(gca, 'ticklength',[tlen1 tlen2])
set(gca, 'tickdir', 'out')
%xlabel('test stim')
ylabel('prop resp CW')


tight_subplot(5,4,2,3, guttera, marginsa)
for  ci = 1:Ncond
    if exp_i == 2
        he(ci) = fill([mean(squeeze(binz_pos(indi_sel,ci,:)),1) mean(squeeze(binz_pos(indi_sel,ci,end:-1:1)),1)],   ...
            [mean(squeeze(prop_cf_pred_all(indi_sel,mi,ci,:)),1)-std(squeeze(prop_cf_pred_all(indi_sel,mi,ci,:)), 1)./sqrt(length(indi_sel))...
            fliplr(mean(squeeze(prop_cf_pred_all(indi_sel,mi,ci,:)),1) + std(squeeze(prop_cf_pred_all(indi_sel,mi,ci,:)),1)./sqrt(length(indi_sel))) ],colorz_shade(ci,:), 'EdgeColor', 'None'); hold on;
    end
    plot(nanmean(squeeze(binz_pos(:,ci,:)),1), nanmean(squeeze(prop_cf_high_all(indi_sel,ci,:)),1),'o', 'Color',colorz(ci,:),'Linewidth', 1.1,'MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:) ); hold on;
    errorbar(nanmean(squeeze(binz_pos(:,ci,:)),1), nanmean(squeeze(prop_cf_high_all(indi_sel,ci,:)),1), nanstd(squeeze(prop_cf_high_all(indi_sel,ci,:)),1)/sqrt(length(indi_sel)), 'Color',  colorz(ci,:),'LineStyle','none','CapSize',0, 'Linewidth', linewi); hold on;
    
end
box off
xlim_min = -0.3;%min(min(min(binz_pos)));
xlim_max = 0.3;%max(max(max(binz_pos)));
ylim_min = 0;
ylim_max = 1.01;
xlim([xlim_min xlim_max])
ylim([ylim_min ylim_max])
set(gca, 'ytick', [0:0.2:1])
plot(zeros(1,nbinz), linspace(0,1,nbinz), '--k'); hold on;
box off
set(gca, 'tickdir', 'out')
set(gca, 'xtick',[-0.3:0.1:0.3])
set(gca, 'ytick', [0:0.25:1])
set(gca, 'tickdir', 'out')
set(gca, 'FontSize', fontsz)
set(gca, 'ticklength',[tlen1 tlen2])
xlabel('test stimulus speed clockwise (deg / sec)')
ylabel('prop high conf')
pars_sub = 0.15;

if exp_i == 2
    tight_subplot(5,4,4,3, guttera, marginsa)
    for pi = 1%5%1 %4 % Bayesian model params
        
        for ci = 1: 4
            hb(ci)=bar(5*(pi-1)+ci, mean(squeeze(params_bm_allV(indi_sel,ci,pi))), 'FaceColor', 'None','EdgeColor',colorz(ci,:), 'Linewidth',linewi); hold on;
            plot((5*(pi-1)+ci)*ones(1, length(indi_sel))-pars_sub+2*pars_sub* rand(1, length(indi_sel)),  params_bm_allV(indi_sel,ci,pi), 'o','MarkerSize',mszi, 'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:)); hold on;
            he(ci) = errorbar(5*(pi-1)+ci, mean(squeeze(params_bm_allV(indi_sel,ci,pi))), std(squeeze(params_bm_allV(indi_sel,ci,pi)))/sqrt(length(indi_sel)),'CapSize',0); hold on;
            set(he(ci),'Color', [0 0 0], 'Linewidth',2)
            errorbarT(he(ci),eb_w,eb_t)
        end
        box off
        set(gca, 'xtick', [1:1:4])
        set(gca, 'tickdir', 'out')
        set(gca, 'FontSize', fontsz)
        set(gca,'ticklength', [tlen1 tlen2])
        set(gca, 'xticklabels', {})
    end
    ylim([-0.2 0.1])
    title('\mu_{encoding} (fixed) ', 'FontName', 'Helvetica', 'FontSize', 1.3*fontsz)
    
    tight_subplot(5,4,5,3, guttera, marginsa)
    for pi = 2 %4 % Bayesian model params
        
        for ci = 1: 4
            hb(ci)=bar(ci, mean(squeeze(params_bm_allV(indi_sel,ci,pi))), 'FaceColor', 'None','EdgeColor',colorz(ci,:), 'Linewidth', linewi); hold on;
            plot((ci)*ones(1, length(indi_sel))-pars_sub+2*pars_sub* rand(1, length(indi_sel)),  params_bm_allV(indi_sel,ci,pi), 'o','MarkerSize',mszi, 'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:)); hold on;
            he(ci) = errorbar(ci, mean(squeeze(params_bm_allV(indi_sel,ci,pi))), std(squeeze(params_bm_allV(indi_sel,ci,pi)))/sqrt(length(indi_sel)),'CapSize',0); hold on;
            set(he(ci),'Color', [0 0 0], 'Linewidth',2)
            errorbarT(he(ci),eb_w,eb_t)
        end
        box off
        set(gca, 'xtick', [1:1:4])
        set(gca, 'tickdir', 'out')
        set(gca, 'FontSize', fontsz)
        set(gca,'ticklength', [tlen1 tlen2])
        set(gca, 'xticklabels', {})
        
    end
    ylim([0 0.3])
    xlabel('\sigma_{encoding}', 'FontName', 'Helvetica', 'FontSize', 1.3*fontsz)
    
    
    tight_subplot(5,4,4,4, guttera, marginsa)
    for pi = 5  % Bayesian model params
        
        for ci = 1: 4
            
            hb(ci)=bar(ci, mean(squeeze(params_bm_all(indi_sel,ci,pi))), 'FaceColor', 'None','EdgeColor',colorz(ci,:), 'Linewidth',linewi); hold on;
            plot((ci)*ones(1, length(indi_sel))-pars_sub+2*pars_sub* rand(1, length(indi_sel)),  params_bm_all(indi_sel,ci,pi), 'o','MarkerSize',mszi, 'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:)); hold on;
            he(ci) = errorbar(ci, mean(squeeze(params_bm_all(indi_sel,ci,pi))), std(squeeze(params_bm_all(indi_sel,ci,pi)))/sqrt(length(indi_sel)),'CapSize',0); hold on;
            set(he(ci),'Color', [0 0 0], 'Linewidth',2)
            errorbarT(he(ci),eb_w,eb_t)
        end
        box off
        set(gca, 'xtick', [1:1:4])
        set(gca, 'tickdir', 'out')
        set(gca, 'FontSize', fontsz)
        set(gca,'ticklength', [tlen1 tlen2])
        set(gca, 'xticklabels', {})
    end
    xlabel('\mu_{likelihood}', 'FontName', 'Helvetica', 'FontSize', 1.3*fontsz)
    
    
    tight_subplot(5,4,5,4, guttera, marginsa)
    for pi = 4  % Bayesian model params
        
        for ci = 1: 4
            hb(ci)=bar(ci, mean(squeeze(params_bm_all(indi_sel,ci,pi))), 'FaceColor', 'None','EdgeColor',colorz(ci,:), 'Linewidth',linewi); hold on;
            plot((ci)*ones(1, length(indi_sel))-pars_sub+2*pars_sub* rand(1, length(indi_sel)),  params_bm_all(indi_sel,ci,pi), 'o','MarkerSize',mszi, 'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:)); hold on;
            he(ci) = errorbar(ci, mean(squeeze(params_bm_all(indi_sel,ci,pi))), std(squeeze(params_bm_all(indi_sel,ci,pi)))/sqrt(length(indi_sel)),'CapSize',0); hold on;
            set(he(ci),'Color', [0 0 0], 'Linewidth',2)
            errorbarT(he(ci),eb_w,eb_t)
        end
        box off
        set(gca, 'tickdir', 'out')
        set(gca, 'xtick', [1:1:4])
        set(gca, 'FontSize', fontsz)
        set(gca,'ticklength', [tlen1 tlen2])
        set(gca, 'xticklabels', {})
    end
    xlabel('k_{confidence}', 'FontName', 'Helvetica', 'FontSize', 1.3*fontsz)
    
    tight_subplot(5,4,5,2, guttera, marginsa)
    for pi = 3  % Bayesian model params
        
        for ci = 1: 4
            hb(ci)=bar(ci, mean(squeeze(params_bm_all(indi_sel,ci,pi))), 'FaceColor', 'None','EdgeColor',colorz(ci,:), 'Linewidth',linewi); hold on;
            plot((ci)*ones(1, length(indi_sel))-pars_sub+2*pars_sub* rand(1, length(indi_sel)),  params_bm_all(indi_sel,ci,pi), 'o','MarkerSize',mszi, 'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:)); hold on;
            he(ci) = errorbar(ci, mean(squeeze(params_bm_all(indi_sel,ci,pi))), std(squeeze(params_bm_all(indi_sel,ci,pi)))/sqrt(length(indi_sel)),'CapSize',0); hold on;
            set(he(ci),'Color', [0 0 0], 'Linewidth',2)
            errorbarT(he(ci),eb_w,eb_t)
        end
        box off
        set(gca, 'tickdir', 'out')
        set(gca, 'xtick', [1:1:4])
        set(gca, 'FontSize', fontsz)
        set(gca,'ticklength', [tlen1 tlen2])
        set(gca, 'xticklabels', {})
    end
    xlabel('k_{choice}', 'FontName', 'Helvetica', 'FontSize', 1.3*fontsz)
    
    
end


%%

%% Bayesian model comparison
nll_fit_all_models = [sum(nll_all,3)]';
AIC = AIC_all';
BIC = BIC_all';
Nmodels = 6;
diff_models_nll = nll_fit_all_models - nll_fit_all_models(1,:);
diff_models_AIC = AIC - AIC(1,:); %bsxfun(@minus,AIC(1,:), AIC );
diff_models_BIC = BIC - BIC(1,:); %bsxfun(@minus,BIC(1,:), BIC );

sample0 = [];
sample=[];
sample2=[];
nboot = 50000;  % for quick run of the script, %1000000 for plotting


for rmi = 1:Nmodels
    d_nll_sum(rmi) = sum(squeeze(diff_models_nll(rmi,:)));
    d_AIC_sum(rmi) = sum(squeeze(diff_models_AIC(rmi,:)));
    d_BIC_sum(rmi) = sum(squeeze(diff_models_BIC(rmi,:)));
    for kk=1:nboot
        
        sample0=randsample(diff_models_nll(rmi,:),Nsubj,1);
        d_nll_sums(kk,rmi) = sum(sample0);
        
        sample=randsample(diff_models_AIC(rmi,:),Nsubj,1);
        d_AIC_sums(kk,rmi) = sum(sample);
        
        sample2=randsample(diff_models_BIC(rmi,:),Nsubj,1);
        d_BIC_sums(kk, rmi) = sum(sample2);
    end
    
end

ci_bnd_low = 0.025;
ci_bnd_high = 0.975;

for rmi = 1:Nmodels
    bci_nll(rmi,1:2) = [quantile(squeeze(d_nll_sums(:,rmi)),ci_bnd_low); quantile(squeeze(d_nll_sums(:,rmi)),ci_bnd_high)];
    bci_aic(rmi,1:2) = [quantile(squeeze(d_AIC_sums(:,rmi)),ci_bnd_low); quantile(squeeze(d_AIC_sums(:,rmi)),ci_bnd_high)];
    bci_bic(rmi,1:2) = [quantile(squeeze(d_BIC_sums(:,rmi)),ci_bnd_low); quantile(squeeze(d_BIC_sums(:,rmi)),ci_bnd_high)];
end

%% add model comparison subplots on figure 5

figure(5)

colorz_modelcomp = repmat([178 178 178],Nmodels, 1)/255;

%FIND COLORZ FOR EACH MODEL

linewi = 1;

tlen1 = 0.014;
tlen2 = 0.014;
fontsz = 12;

ylim_min_1 = -800;
ylim_min_2 = -50;


tight_subplot(5,4,1,1, guttera, marginsa)
for mi = 1:Nmodels
    bar(mi, d_AIC_sum(mi), 'FaceColor', colorz_modelcomp(mi,:), 'EdgeColor', colorz_modelcomp(mi,:)); hold on;
    errorbar(mi, d_AIC_sum(mi),d_AIC_sum(mi)-bci_aic(mi,1),bci_aic(mi,2)-d_AIC_sum(mi), 'Color', 'k','Capsize', 0, 'LineWidth',linewi); hold on;
end
xlim([0.5 Nmodels+0.5])
ylim([ylim_min_1 300])
%ylim([480 580])
box off
set(gca, 'tickdir', 'out')
set(gca, 'FontSize', fontsz)
set(gca, 'ticklength',[tlen1 tlen2])
set(gca, 'xtick', 1:1:Nmodels)
set(gca, 'xticklabels', {'M1','M2', 'M3','M4', 'M5', 'M6' })
%ylabel('AIC(model) - AIC(winning model)')
ylabel('\Delta AIC')
title('\Delta AIC = AIC(Mi) - AIC(M1)')


tight_subplot(5,4,1,2, guttera, marginsa)
for mi = 1:Nmodels
    bb(mi)=bar(mi, d_BIC_sum(mi),  'FaceColor', colorz_modelcomp(mi,:), 'EdgeColor', colorz_modelcomp(mi,:)); hold on;
    ee(mi) = errorbar(mi, d_BIC_sum(mi),d_BIC_sum(mi)-bci_bic(mi,1),bci_bic(mi,2)-d_BIC_sum(mi) ,'Color', 'k', 'Capsize', 0, 'LineWidth',linewi); hold on;
end

ll5 = legend([bb(1), bb(2), bb(3) bb(4) bb(5) bb(6)], {'fit \mu_{likelihood}', 'fit prior', 'fit \mu_{likelihood} and prior', ...
    'fit k_{choice}', 'fit \mu_{likelihood} and k_{choice}',  'fit k_{choice and confidence}' }, 'FontName', 'Helvetica', 'FontSize', 9)
set(ll5, 'Position', [0.2108    0.5929    0.1811    0.1357]);
%ylim([480 580])
box off
set(gca, 'tickdir', 'out')
set(gca, 'FontSize', fontsz)
set(gca, 'ticklength',[tlen1 tlen2])
set(gca, 'xtick', 1:1:Nmodels)
set(gca, 'xticklabels', {'M1','M2', 'M3', 'M4', 'M5', 'M6'})
%ylabel('BIC(model) - BIC(winning model)')
ylabel('\Delta BIC')
xlim([0.5 Nmodels+0.5])
ylim([ylim_min_1 300])
%%
psname = ['Fig5_parts_June2023_exp_',num2str(exp_i) ,'CTR_Nsubj_',num2str(Nsubj),'.pdf'];
print_pdf(psname)



%% correlations of model (1) parameters with MAE compensation
delta_mu_lik = squeeze(params_bm_all(:,4,4))- squeeze(params_bm_all(:,3,4))
[r,p] = corr( compensation, delta_mu_lik, 'type', 'Spearman','rows','complete')
%r = -0.95, 10^-11

delta_mu_noise = exp(squeeze(params_bm_all(:,4,2)))- exp(squeeze(params_bm_all(:,3,2)))
[r,p] = corr( compensation, delta_mu_noise, 'type', 'Spearman','rows','complete')
%r = 0.1073, p = 0.6345
%[r,p] = corr( compensation_norm, delta_prior_norm, 'type', 'Spearman','rows','complete')

delta_mu_k_conf = squeeze(params_bm_all(:,4,3))- squeeze(params_bm_all(:,3,3));
[r,p] = corr(compensation, delta_mu_k_conf, 'type', 'Spearman','rows','complete')
%r = 0.03, p =  0.88

%{
%% GLMEs for RT with stim strength and RT with the perceptual certainty |d| - results associated with Figure 5D

load('stims_set_load.mat')
Ntrials = 121;
cond_for_lme = [];
stim_strength_for_lme = [];
stim_strength_rel_for_lme = [];
dec_var_for_lme = [];
rt_for_lme = [];
subject_for_lme = [];
rt_all_tr = rt_all_trials;
for i = 1: Nsubj
    subject_for_lme = [subject_for_lme; i*ones(1, Ntrials*Ncond)'];
    cond_for_lme = [cond_for_lme; [ones(1,121) 2*ones(1,121) 3*ones(1,121) 4*ones(1,121)]'];
    stim_strength_for_lme = [ stim_strength_for_lme; [repmat(abs(stims_set_load),1,4)]'];
    rt_for_lme = [rt_for_lme; squeeze(rt_all_tr(i,1,:)); squeeze(rt_all_tr(i,2,:)); squeeze(rt_all_tr(i,3,:)); squeeze(rt_all_tr(i,4,:))];
    dec_var_for_lme = [dec_var_for_lme; squeeze(abs(dv_val(i,1,:))); squeeze(abs(dv_val(i,2,:))); squeeze(abs(dv_val(i,3,:))); squeeze(abs(dv_val(i,4,:)))];
    %rt_for_lme = [rt_for_lme; squeeze(rt_all_tr(i,3,:)); squeeze(rt_all_tr(i,4,:)); squeeze(rt_all_tr(i,1,:)); squeeze(rt_all_tr(i,2,:))];
    %dec_var_for_lme = [dec_var_for_lme; squeeze(abs(dv_val(i,3,:))); squeeze(abs(dv_val(i,4,:))); squeeze(abs(dv_val(i,1,:))); squeeze(abs(dv_val(i,2,:)))];
end
%%

[~, ~, stim_strength_for_lme_ranking] = unique(stim_strength_for_lme);
[~, ~, dec_var_for_lme_ranking] = unique(dec_var_for_lme);
[~, ~, rt_for_lme_ranking] = unique(rt_for_lme);

tbl_abs_stim = table(stim_strength_for_lme_ranking,rt_for_lme_ranking,cond_for_lme, subject_for_lme,'VariableNames',{'StimStrength','RT','Condition', 'Subject'});
tbl_abs_stim.Condition = nominal(tbl_abs_stim.Condition);
tbl_abs_stim.Subject = nominal(tbl_abs_stim.Subject);

tbl_abs_dv = table(dec_var_for_lme_ranking,rt_for_lme_ranking,cond_for_lme, subject_for_lme,'VariableNames',{'DecVar','RT','Condition', 'Subject'});
tbl_abs_dv.Condition = nominal(tbl_abs_dv.Condition);
tbl_abs_dv.Subject = nominal(tbl_abs_dv.Subject);
%%

lme_abs_stim = fitlme(tbl_abs_stim,'RT ~ 1 + StimStrength+Condition + (1+StimStrength+Condition|Subject)');
lme_abs_dv = fitlme(tbl_abs_dv,'RT ~ 1 + DecVar+Condition + (1+DecVar+Condition|Subject)');
%%
% AIC and BIC
[lme_abs_stim.ModelCriterion.AIC lme_abs_stim.ModelCriterion.BIC;...
    lme_abs_dv.ModelCriterion.AIC lme_abs_dv.ModelCriterion.BIC]

% coefficients and tStat and p values
[lme_abs_stim.Coefficients(2:end,1) lme_abs_stim.Coefficients(2:end,4)  lme_abs_stim.Coefficients(2:end,6)]
[lme_abs_dv.Coefficients(2:end,1) lme_abs_dv.Coefficients(2:end,4)  lme_abs_dv.Coefficients(2:end,6)]


%% same GLME analysis only for Adapt-See and Adapt-Believe

cond_for_lme_R = [];
stim_strength_for_lme_R = [];
stim_strength_rel_for_lme_R = [];
dec_var_for_lme_R= [];
rt_for_lme_R = [];
subject_for_lme_R = [];
%rt_all_tr = rt_all_trials;
for i = 1: Nsubj
    subject_for_lme_R = [subject_for_lme_R; i*ones(1, Ntrials*2)'];
    cond_for_lme_R = [cond_for_lme_R; [3*ones(1,121) 4*ones(1,121)]'];
    stim_strength_for_lme_R = [ stim_strength_for_lme_R; [repmat(abs(stims_set_load),1,2)]'];
    rt_for_lme_R = [rt_for_lme_R; squeeze(rt_all_tr(i,3,:)); squeeze(rt_all_tr(i,4,:))];
    dec_var_for_lme_R = [dec_var_for_lme_R; squeeze(abs(dv_val(i,3,:))); squeeze(abs(dv_val(i,4,:)))];
    %rt_for_lme = [rt_for_lme; squeeze(rt_all_tr(i,3,:)); squeeze(rt_all_tr(i,4,:)); squeeze(rt_all_tr(i,1,:)); squeeze(rt_all_tr(i,2,:))];
    %dec_var_for_lme = [dec_var_for_lme; squeeze(abs(dv_val(i,3,:))); squeeze(abs(dv_val(i,4,:))); squeeze(abs(dv_val(i,1,:))); squeeze(abs(dv_val(i,2,:)))];
end
%%

[~, ~, stim_strength_for_lme_ranking_R] = unique(stim_strength_for_lme_R);
[~, ~, dec_var_for_lme_ranking_R] = unique(dec_var_for_lme_R);
[~, ~, rt_for_lme_ranking_R] = unique(rt_for_lme_R);

tbl_abs_stim_R = table(stim_strength_for_lme_ranking_R,rt_for_lme_ranking_R,cond_for_lme_R, subject_for_lme_R,'VariableNames',{'StimStrength','RT','Condition', 'Subject'});
tbl_abs_stim_R.Condition = nominal(tbl_abs_stim_R.Condition);
tbl_abs_stim_R.Subject = nominal(tbl_abs_stim_R.Subject);

tbl_abs_dv_R = table(dec_var_for_lme_ranking_R,rt_for_lme_ranking_R,cond_for_lme_R, subject_for_lme_R,'VariableNames',{'DecVar','RT','Condition', 'Subject'});
tbl_abs_dv_R.Condition = nominal(tbl_abs_dv_R.Condition);
tbl_abs_dv_R.Subject = nominal(tbl_abs_dv_R.Subject);
%%

lme_abs_stim_R = fitlme(tbl_abs_stim_R,'RT ~ 1 + StimStrength+Condition + (1+StimStrength+Condition|Subject)')
lme_abs_dv_R = fitlme(tbl_abs_dv_R,'RT ~ 1 + DecVar+Condition + (1+DecVar+Condition|Subject)')


%% confidence-consistency relationship, as in Boundy-Singer et al, 2022
figure
set(gcf, 'Position', [100 100 560 280])
gutteraa = [ 0.0900    0.0600];
marginsaa = [0.0800    0.0800    0.1400    0.1000]; %MARGINS = [LEFT RIGHT BOTTOM TOP]
%for si = 1: Nsubj
tight_subplot(1,2,1,1, gutteraa, marginsaa)
Nsubj = 22;
for bi = 1:11
    plot(nanmean(prop_cw_all(:,3,bi)), nanmean(prop_cf_high_all(:,3,bi)), 'o-', 'MarkerFaceColor', colorz(3,:),  'MarkerEdgeColor', colorz(3,:)); hold on;
    errorbar(nanmean(prop_cw_all(:,3,bi)), nanmean(prop_cf_high_all(:,3,bi)), nanstd(prop_cw_all(:,3,bi))/sqrt(Nsubj),nanstd(prop_cw_all(:,3,bi))/sqrt(Nsubj) , 'linestyle', 'none', 'color', colorz(3,:), 'linewidth', .8);
    errorbar(nanmean(prop_cw_all(:,3,bi)), nanmean(prop_cf_high_all(:,3,bi)), nanstd(prop_cf_high_all(:,3,bi))/sqrt(Nsubj),nanstd(prop_cf_high_all(:,3,bi))/sqrt(Nsubj) , 'horizontal', 'linestyle', 'none', 'color', colorz(3,:), 'linewidth', .8);
    
    plot(nanmean(prop_cw_all(:,4,bi)), nanmean(prop_cf_high_all(:,4,bi)), 'o-', 'MarkerFaceColor', colorz(4,:),  'MarkerEdgeColor', colorz(4,:)); hold on;
    errorbar(nanmean(prop_cw_all(:,4,bi)), nanmean(prop_cf_high_all(:,4,bi)), nanstd(prop_cw_all(:,4,bi))/sqrt(Nsubj),nanstd(prop_cw_all(:,4,bi))/sqrt(Nsubj) , 'linestyle', 'none', 'color', colorz(4,:), 'linewidth', .8);
    errorbar(nanmean(prop_cw_all(:,4,bi)), nanmean(prop_cf_high_all(:,4,bi)), nanstd(prop_cf_high_all(:,4,bi))/sqrt(Nsubj),nanstd(prop_cf_high_all(:,4,bi))/sqrt(Nsubj) , 'horizontal', 'linestyle', 'none', 'color', colorz(4,:), 'linewidth', .8);
end
%errorbar(x, y, error_y, error_y, 'horizontal', 'linestyle', 'none', 'color', 'k', 'linewidth', .8);
%end

set(gca, 'tickdir', 'out')
box off
xlim([0 1.05])
ylim([0 1.05])
vala1 = squeeze((nanmean(prop_cw_all(:,3,:),1)));
vala2 = squeeze((nanmean(prop_cf_high_all(:,3,:),1)));
vala3 = squeeze((nanmean(prop_cw_all(:,4,:),1)));
vala4 = squeeze((nanmean(prop_cf_high_all(:,4,:),1)));
plot(vala1, vala2, 'Color', colorz(3,:)); hold on;
plot(vala3, vala4, 'Color', colorz(4,:)); hold on;
xlabel('proportion response clockwise', 'FontName', 'Helvetica', 'FontSize', fontsz)
ylabel('proportion high confidence','FontName', 'Helvetica', 'FontSize', fontsz)
title(['Nsubj = ', num2str(Nsubj)], 'FontName', 'Helvetica', 'FontSize', fontsz)

tight_subplot(1,2,1,2, gutteraa, marginsaa)
Nsubj = 22;
for bi = 1:11
    plot(nanmean(prop_cw_all(:,3,bi)), nanmean(rt_all(:,3,bi)), 'o-', 'MarkerFaceColor', colorz(3,:),  'MarkerEdgeColor', colorz(3,:)); hold on;
    errorbar(nanmean(prop_cw_all(:,3,bi)), nanmean(rt_all(:,3,bi)), nanstd(prop_cw_all(:,3,bi))/sqrt(Nsubj),nanstd(prop_cw_all(:,3,bi))/sqrt(Nsubj) , 'linestyle', 'none', 'color', colorz(3,:), 'linewidth', .8);
    errorbar(nanmean(prop_cw_all(:,3,bi)), nanmean(rt_all(:,3,bi)), nanstd(rt_all(:,3,bi))/sqrt(Nsubj),nanstd(rt_all(:,3,bi))/sqrt(Nsubj) , 'horizontal', 'linestyle', 'none', 'color', colorz(3,:), 'linewidth', .8);
    
    plot(nanmean(prop_cw_all(:,4,bi)), nanmean(rt_all(:,4,bi)), 'o-', 'MarkerFaceColor', colorz(4,:),  'MarkerEdgeColor', colorz(4,:)); hold on;
    errorbar(nanmean(prop_cw_all(:,4,bi)), nanmean(rt_all(:,4,bi)), nanstd(prop_cw_all(:,4,bi))/sqrt(Nsubj),nanstd(prop_cw_all(:,4,bi))/sqrt(Nsubj) , 'linestyle', 'none', 'color', colorz(4,:), 'linewidth', .8);
    errorbar(nanmean(prop_cw_all(:,4,bi)), nanmean(rt_all(:,4,bi)), nanstd(rt_all(:,4,bi))/sqrt(Nsubj),nanstd(rt_all(:,4,bi))/sqrt(Nsubj) , 'horizontal', 'linestyle', 'none', 'color', colorz(4,:), 'linewidth', .8);
end

set(gca, 'tickdir', 'out')
box off
xlim([0 1.05])
vala1 = squeeze((nanmean(prop_cw_all(:,3,:),1)));
vala2 = squeeze((nanmean(rt_all(:,3,:),1)));
vala3 = squeeze((nanmean(prop_cw_all(:,4,:),1)));
vala4 = squeeze((nanmean(rt_all(:,4,:),1)));
plot(vala1, vala2, 'Color', colorz(3,:)); hold on;
plot(vala3, vala4, 'Color', colorz(4,:)); hold on;
xlabel('proportion response clockwise', 'FontName', 'Helvetica', 'FontSize', fontsz)
ylabel('reaction times (s)','FontName', 'Helvetica', 'FontSize', fontsz)

psname =['confidence_consistency_relationship_exp', num2str(exp_i),'.pdf']
%print_pdf(psname)

% conditional response functions - generate one figure per subject
% rt_all_tr
%divide into 5 quantiles --find proportion resp yes for each quantile
% this is across several stimulus strengths


for si = 1 : Nsubj
    
    datac = alldata(si,:);
    
    for ci = 1:Ncond
        % define the RT quantiles
        rtz = rt_all_trials(si,ci,1:121);
        
        binz_rt = [];
        for j = 1:4
            binz_rt(j) = quantile(rt_all_trials(si,ci,1:121), j/5);
        end
        binz_rt(5) = quantile(rt_all_trials(si,ci,1:121), 0.95);
        
        binz_rt = [0 binz_rt ];
        binz_posE_rt(si,ci,:) = (binz_rt(2:end)+binz_rt(1:end-1))/2;
        binzz_rt(si,ci,:) = binz_rt;
        
        for j = 1:5
            ind_rt = squeeze((rtz>binz_rt(j) & rtz<=binz_rt(j+1)));
            %ind_cw =  datac(ci).stims >= 0;
            %ind_ccw = datac(ci).stims <= 0;
            ind_cw =  datac(ci).conf == 1;
            ind_ccw = datac(ci).conf == 0;
            if sum(ind_rt)>5
                prop_cw_rt(si,ci,j) = nansum(datac(ci).resp(ind_rt))/sum(~isnan(datac(ci).resp(ind_rt)));
            else
                prop_cw_rt(si,ci,j)  = NaN;
            end
            
            ind_rt_cw = ind_rt & ind_cw';
            if sum(ind_rt_cw)>5
                prop_cw_rt_cw(si,ci,j) = nansum(datac(ci).resp(ind_rt_cw))/sum(~isnan(datac(ci).resp(ind_rt_cw)));
            else
                prop_cw_rt_cw(si,ci,j)  = NaN;
            end
            
            ind_rt_ccw = ind_rt & ind_ccw';
            if sum(ind_rt_ccw)>5
                prop_cw_rt_ccw(si,ci,j) = nansum(datac(ci).resp(ind_rt_ccw))/sum(~isnan(datac(ci).resp(ind_rt_ccw)));
            else
                prop_cw_rt_ccw(si,ci,j)  = NaN;
            end
        end
        %{
        figure(si)
        set(gcf, 'Position', [100 100 600 240])
        tight_subplot(1,3,1,1, gutteraa, marginsaa)
        plot( squeeze(binz_posE_rt(si,ci,:)), squeeze(prop_cw_rt(si,ci,:)),'-', 'Color', colorz(ci,:), 'Linewidth',2); hold on; %'MarkerEdgeColor', colorz(ci,:)
        ylim([0 1])
        box off
        if ci == 1
            xlabel('RT quantile')
            ylabel('prop resp CW')
        end
        %}
    end
end


% conditional response functions - ALL subjects, RT plots everything
figure(7)
set(gcf, 'Position', [100 100 700 540])
guttera_rt = [0.0800    0.0800];
marginsa_rt =  [0.100    0.0700    0.100    0.1900];
tight_subplot(3,4,1,1, guttera_rt, marginsa_rt)
for ci = 1: Ncond
    plot( mean(squeeze(binz_posE_rt(:,ci,:)),1), mean(squeeze(prop_cw_rt(:,ci,:)),1),'-', 'Color', colorz(ci,:), 'Linewidth',2); hold on; %'MarkerEdgeColor', colorz(ci,:)
    errorbar( mean(squeeze(binz_posE_rt(:,ci,:)),1), mean(squeeze(prop_cw_rt(:,ci,:)),1),std(squeeze(prop_cw_rt(:,ci,:)),1)/sqrt(Nsubj),'-', 'Color', colorz(ci,:), 'Linewidth',2); hold on; %'MarkerEdgeColor', colorz(ci,:)
    ylim([0 1])
end
plot(linspace(0,3,20), 0.5*ones(1,20), '--k', 'Linewidth', 1.2); hold on;
box off
%xlabel('RT quantile')
xlabel('Reaction time  (sec)')
ylabel('prop resp CW')
set(gca, 'tickdir', 'out')

tight_subplot(3,4,1,2, guttera_rt, marginsa_rt)
% schematic for the starting point and drift bias
lm1 = plot( mean(squeeze(binz_posE_rt(:,ci,:)),1), [0.53 0.51 0.56 0.56 0.53], '-', 'Color', [0.5 0.5 0.5],'Linewidth',2); hold on;
lm2 = plot( mean(squeeze(binz_posE_rt(:,ci,:)),1), [0.8 0.61 0.52 0.53 0.53], '-', 'Color', [0.9 0.0 0.9],'Linewidth',2); hold on;
lm3 = plot( mean(squeeze(binz_posE_rt(:,ci,:)),1), [0.8 0.76 0.72 0.68 0.64], '-', 'Color', [0.4 0.0 0.9],'Linewidth',2); hold on;

plot( mean(squeeze(binz_posE_rt(:,ci,:)),1), [0.2 0.43 0.45 0.49 0.49], '--', 'Color', [0.9 0.0 0.9],'Linewidth',2); hold on;
plot( mean(squeeze(binz_posE_rt(:,ci,:)),1), [0.22 0.26 0.30 0.34 0.38], '--', 'Color', [0.4 0.0 0.9],'Linewidth',2); hold on;
ylim([0 1])
box off
plot(linspace(0,3,20), 0.5*ones(1,20), '--k', 'Linewidth', 1.2); hold on;
set(gca, 'tickdir', 'out')
set(gca, 'FontSize', fontsz)
set(gca, 'ticklength',[tlen1 tlen2])
llm = legend([lm1 lm2 lm3], {'Unbiased', 'Starting point', 'Drift bias'});
set(llm, 'Position',  [0.33    0.45    0.15    0.08])
% model comparison subplot
load('ddm_params_and_model_comp.mat')
% delta AIC and BIC from base model

lw = 1.5;
fontsz = 12;
pars_scatter = 0.15;
mszi = 2;

AIC_full = sum(AIC_all_ddm,3);
BIC_full = sum(BIC_all_ddm,3);

diff_models_AIC = AIC_full - AIC_full(2,:);
diff_models_BIC = BIC_full - BIC_full(2,:);

sample0 = [];
sample=[];
sample2=[];
nboot = 100000;%1000000 when generating figure

for rmi = 1:4
    
    d_AIC_sum(rmi) = sum(squeeze(diff_models_AIC(rmi,:)));
    d_BIC_sum(rmi) = sum(squeeze(diff_models_BIC(rmi,:)));
    for kk=1:nboot
        
        sample=randsample(diff_models_AIC(rmi,:),Nsubj,1);
        d_AIC_sums(kk,rmi) = sum(sample);
        
        sample2=randsample(diff_models_BIC(rmi,:),Nsubj,1);
        d_BIC_sums(kk, rmi) = sum(sample2);
    end
    
end

ci_bnd_low = 0.025;
ci_bnd_high = 0.975;

for rmi = 1:4
    bci_aic(rmi,1:2) = [quantile(squeeze(d_AIC_sums(:,rmi)),ci_bnd_low); quantile(squeeze(d_AIC_sums(:,rmi)),ci_bnd_high)];
    bci_bic(rmi,1:2) = [quantile(squeeze(d_BIC_sums(:,rmi)),ci_bnd_low); quantile(squeeze(d_BIC_sums(:,rmi)),ci_bnd_high)];
end


colorz_mod_comp = repmat([178 178 178],5, 1)/255;
linewi = 1;
tlen1 = 0.014;
tlen2 = 0.014;
%linewi = 1.1;
fontsz = 12;
ylim_min_1 = -200;%-1400;
ylim_min_2 = -200;

tight_subplot(3,4,1,3, guttera_rt, marginsa_rt)
for mi = 1:4
    bar(mi, d_AIC_sum(mi), 'FaceColor', colorz_mod_comp(mi,:), 'EdgeColor', colorz_mod_comp(mi,:)); hold on;
    errorbar(mi, d_AIC_sum(mi),d_AIC_sum(mi)-bci_aic(mi,1),bci_aic(mi,2)-d_AIC_sum(mi), 'Color', 'k','Capsize', 0, 'LineWidth',linewi); hold on;
end
xlim([0.5 4.5])
box off
set(gca, 'tickdir', 'out')
set(gca, 'FontSize', fontsz)
set(gca, 'ticklength',[tlen1 tlen2])
set(gca, 'xtick', 1:1:4)
set(gca, 'xticklabels', [])
title('\Delta AIC')
ylim([ylim_min_1 1400])

tight_subplot(3,4,1,4, guttera_rt, marginsa_rt)
for mi = 1:4
    bb(mi) = bar(mi, d_BIC_sum(mi),  'FaceColor', colorz_mod_comp(mi,:), 'EdgeColor', colorz_mod_comp(mi,:)); hold on;
    ee(mi) = errorbar(mi, d_BIC_sum(mi),d_BIC_sum(mi)-bci_bic(mi,1),bci_bic(mi,2)-d_BIC_sum(mi) ,'Color', 'k', 'Capsize', 0, 'LineWidth',linewi); hold on;
end
ll = legend([bb(1), bb(2), bb(3) bb(4)], {'Base', 'Starting point', 'Drift bias', ...
    'Starting point + drift bias'});
set(ll, 'Position',  [0.27    0.3    0.23    0.11])
box off
set(gca, 'tickdir', 'out')
set(gca, 'FontSize', fontsz)
set(gca, 'ticklength',[tlen1 tlen2])
set(gca, 'xtick', 1:1:4)
set(gca, 'xticklabels', [])
title('\Delta BIC')
xlim([0.5 4.5])
ylim([ylim_min_1 1400])


ddm_plot = 1;
%if ddm_plot == 1
load('output_list_RT_allEE.mat');
load('output_list_choices_allEE.mat');
prop_cw_all_ddm = NaN(Nsubj, Ncond, nbinz);
rt_all_ddm = NaN(Nsubj, Ncond, nbinz);
RT_pred_binz = NaN(Nsubj, Ncond, nbinz);
prop_cw_pred_binz = NaN(Nsubj, Ncond, nbinz);

rt_all_left_ddm = NaN(Nsubj, Ncond, nbinz);
rt_all_right_ddm = NaN(Nsubj, Ncond, nbinz);
RT_pred_left_binz = NaN(Nsubj, Ncond, nbinz);
RT_pred_right_binz = NaN(Nsubj, Ncond, nbinz);

binz_poss = (binz(2:end)+binz(1:end-1))/2;
for si = 1: 22
    datac = alldata(si,1:4);
    for ci = 1: Ncond
        
        RT_pred = output_list_RT(4*(si-1)+ci,1:121);
        prop_corr_pred = output_list_choices(4*(si-1)+ci,1:121);
        prop_CW_pred= NaN(1,Ntrials);
        prop_CW_pred(datac(ci).stims<0) = 1- prop_corr_pred(datac(ci).stims<0);
        %prop_CW_pred(datac(ci).stims<=0) = 1- prop_corr_pred(datac(ci).stims<=0);
        prop_CW_pred(datac(ci).stims>=0) = prop_corr_pred(datac(ci).stims>=0);
        
        prop_CW_predE = prop_CW_pred;
        prop_CW_pred(prop_CW_pred>0.5) = 1;
        prop_CW_pred(prop_CW_pred<0.5) = 0;
        
        for j = 1:(nbinz)
            indo = (datac(ci).stims>binz(j) & datac(ci).stims<=binz(j+1));
            indi = find(datac(ci).stims>binz(j) & datac(ci).stims<=binz(j+1) );
            indo2 = ~isnan(RT_pred);
            ind_sel = indo & indo2;
            indi_right = (datac(ci).resp == 1);
            indi_left = (datac(ci).resp == 0);
            if length(ind_sel)>10
                prop_cw_all_ddm(si,ci,j) = nansum(datac(ci).resp(ind_sel))/sum(~isnan(datac(ci).resp(ind_sel)));
                %{
            if binz_pos(j)<0  % correct would mean counterclockwise
                prop_cw_pred_binz(si,ci,j) = 1-nanmean(prop_corr_pred(indi));
            elseif binz_pos(j)>0 %correct means clockwise
                prop_cw_pred_binz(si,ci,j) = nanmean(prop_corr_pred(indi));
            end
                %}
                prop_cw_pred_binz(si,ci,j) = nanmean(prop_CW_pred(ind_sel));
                
                rt_all_ddm(si,ci,j) = nanmedian(datac(ci).resp_times(ind_sel));
                if length(ind_sel & indi_left) >10%0
                    rt_all_left_ddm(si,ci,j) = nanmedian(datac(ci).resp_times(ind_sel & indi_left));
                    RT_pred_left_binz(si,ci,j) = nanmedian(RT_pred(ind_sel & indi_left));
                end
                if length(ind_sel & indi_right) >10%0
                    rt_all_right_ddm(si,ci,j) = nanmedian(datac(ci).resp_times(ind_sel & indi_right));
                    RT_pred_right_binz(si,ci,j) = nanmedian(RT_pred(ind_sel & indi_right));
                end
                RT_pred_binz(si,ci,j) = nanmedian(RT_pred(ind_sel)); %nanmedian(datac(ci).resp_times(indi)); % maybe nanmean? % before it was median
            else
                prop_cw_all_ddm(si,ci,j) = nan;
                prop_cw_pred_binz(si,ci,j) = nan;
                rt_all_ddm(si,ci,j) = nan;
                RT_pred_binz(si,ci,j)= nan;
            end
        end
    end
end

tight_subplot(3,4,2,1, guttera_rt, marginsa_rt)
for  ci = 1:Ncond
    
    he(ci) = fill([nanmean(squeeze(binz_pos(indi_sel,ci,:)),1) mean(squeeze(binz_pos(indi_sel,ci,end:-1:1)),1)],   ...
        [nanmean(squeeze(prop_cw_pred_binz(indi_sel,ci,:)),1)-nanstd(squeeze(prop_cw_pred_binz(indi_sel,ci,:)),1)./sqrt(length(indi_sel))...
        fliplr(nanmean(squeeze(prop_cw_pred_binz(indi_sel,ci,:)),1) + nanstd(squeeze(prop_cw_pred_binz(indi_sel,ci,:)),1)./sqrt(length(indi_sel))) ],colorz_shade(ci,:), 'EdgeColor', 'None'); hold on;
    
    plot(nanmean(squeeze(binz_pos(:,ci,:)),1), nanmean(squeeze(prop_cw_all_ddm(indi_sel,ci,:)),1),'o','MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:) ); hold on;
    errorbar(nanmean(squeeze(binz_pos(:,ci,:)),1), nanmean(squeeze(prop_cw_all_ddm(indi_sel,ci,:)),1), nanstd(squeeze(prop_cw_all_ddm(indi_sel,ci,:)),1)/sqrt(length(indi_sel)), 'Color',  colorz(ci,:), 'LineStyle', 'none','CapSize',0,'Linewidth',linewi); hold on;
    
end
box off
xlim_min = -0.3;%min(min(min(binz_pos)));
xlim_max = 0.3;%max(max(max(binz_pos)));
ylim_min = 0;
ylim_max = 1.01;
xlim([xlim_min xlim_max])
ylim([ylim_min ylim_max])
plot(zeros(1,nbinz), linspace(0,1,nbinz), '--k'); hold on;
set(gca, 'ticklength',[tlen1 tlen2])
box off
set(gca, 'xtick', [-0.3:0.1:0.3])
set(gca, 'ytick', [0:0.25:1])
set(gca, 'tickdir', 'out')
set(gca, 'FontSize', fontsz)
set(gca, 'ticklength',[tlen1 tlen2])
set(gca, 'tickdir', 'out')
%xlabel('test stim')
ylabel('prop resp CW')


tight_subplot(3,4,3,1,guttera_rt, marginsa_rt)
for  ci = [1 2 4 3]%1:Ncond
    
    he(ci) = fill([binz_poss binz_poss(end:-1:1)],   ...
        [nanmean(squeeze(RT_pred_binz(:,ci,:)),1)-nanstd(squeeze(RT_pred_binz(:,ci,:)),1)./sqrt(Nsubj)...
        fliplr(nanmean(squeeze(RT_pred_binz(:,ci,:)),1) + nanstd(squeeze(RT_pred_binz(:,ci,:)),1)./sqrt(Nsubj)) ],colorz_shade(ci,:), 'EdgeColor', 'None'); hold on;
    h(ci) = plot(mean(squeeze(binz_pos(indi_sel,ci,:)),1), nanmean(squeeze(rt_all_ddm(indi_sel,ci,:)),1),'o','MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:)); hold on;
    errorbar(mean(squeeze(binz_pos(indi_sel,ci,:)),1), nanmean(squeeze(rt_all_ddm(indi_sel,ci,:)),1), nanstd(squeeze(rt_all_ddm(indi_sel,ci,:)),1)/sqrt(length(indi_sel)), 'Color',  colorz(ci,:),'LineStyle', 'none','CapSize',0, 'Linewidth', linewi); hold on;
    
end
box off
xlim([xlim_min xlim_max])
ylim([1 2])
plot(zeros(1,nbinz), linspace(0,4,nbinz), '--k'); hold on;
box off
set(gca, 'xtick',[-0.3:0.1:0.3])
set(gca, 'tickdir', 'out')
set(gca, 'FontSize', fontsz)
set(gca, 'ticklength',[tlen1 tlen2])
ylabel('reaction times (s)')
xlabel('test stimulus speed clockwise (a.u.)', 'FontName', 'Helvetica', 'FontSize', fontsz)
%end

mi = 2;
for pi = 1:4
    
    tight_subplot(3,4, 2+mod(1+pi,2), 3+floor((-0.5+pi)/2),guttera_rt, marginsa_rt)
    for ci = 1:4
        bar(ci, mean(params_ddm(mi,:,ci,pi)),'FaceColor', 'none', 'EdgeColor',colorz(ci,:),'Linewidth',lw); hold on;
        errorbar(ci, mean(params_ddm(mi,:,ci,pi)), std(params_ddm(mi,:,ci,pi))/sqrt(Nsubj), 'Color','k','Linewidth',2, 'CapSize', 0); hold on;
        plot(ci*ones(1, Nsubj)- pars_scatter+ 2*pars_scatter*rand(1,Nsubj),  squeeze(params_ddm(mi,:,ci,pi)), 'o','MarkerSize',mszi, 'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:)); hold on;
    end
    box off
    %ylim([-25 5])
    xlim([0.2 4.8])
    set(gca, 'FontSize', fontsz)
    set(gca, 'xtick', [1:1:4])
    set(gca, 'xticklabels', [])
    set(gca, 'tickdir', 'out')
    set(gca, 'ticklength',[tlen1 tlen2])
    if pi == 1
        title(' Mean drift rate', 'FontName', 'Helvetica', 'FontSize', fontsz)
    elseif pi == 2
        title(' Decision bound', 'FontName', 'Helvetica', 'FontSize', fontsz)
    elseif pi == 3
        title(' Nondecision time', 'FontName', 'Helvetica', 'FontSize', fontsz)
    elseif pi == 4
        title(' Starting point bias', 'FontName', 'Helvetica', 'FontSize', fontsz)
    end
end


psname = 'RT_plot_all.pdf'
%print_pdf(psname)

% psychometric curves divided into high conf and low conf
prop_cw_high_conf = NaN(Nsubj, Ncond, nbinz);
% binz for stim is always the same
mi = 1;
for si = 1 : Nsubj
    
    datac = alldata(si,:);
    
    for ci = 1:Ncond
        %datac = alldata(si,ci);
        prc_si = Predict_bm_alll(squeeze(params_bm_all(si,ci,1:4)),datac(ci).stims, datac(ci).resp, datac(ci).conf, mi,1,N_samp);
        
        prop_cw_predd_all(si,ci,:) = prc_si(:,1)+ prc_si(:,2);
        prop_cf_predd_all(si,ci,:) = prc_si(:,1)+ prc_si(:,3);
        
        for j = 1:(nbinz)
            indi = datac(ci).stims>binz(j) & datac(ci).stims<=binz(j+1) ;
            ind_cw =  datac(ci).conf == 1;
            ind_ccw = datac(ci).conf == 0;
            
            ind_high_conf = indi & ind_cw;
            if sum(ind_high_conf)>3
                prop_cw_high_conf(si,ci,j) = nansum(datac(ci).resp(ind_high_conf))/sum(~isnan(datac(ci).resp(ind_high_conf)));
                prop_cw_high_conf_pred(si,ci,j) = nanmean(prop_cw_predd_all(si,ci,ind_high_conf));
            else
                prop_cw_high_conf(si,ci,j)  = NaN;
                prop_cw_high_conf_pred(si,ci,j) = NaN;
            end
            
            ind_low_conf = indi & ind_ccw;
            if sum(ind_low_conf)>3
                prop_cw_low_conf(si,ci,j) = nansum(datac(ci).resp(ind_low_conf))/sum(~isnan(datac(ci).resp(ind_low_conf)));
                prop_cw_low_conf_pred(si,ci,j) = nanmean(prop_cw_predd_all(si,ci,ind_low_conf));
            else
                prop_cw_low_conf(si,ci,j)  = NaN;
                prop_cw_low_conf_pred(si,ci,j) = NaN;
            end
        end
        
    end
    %{
    figure(si)
    set(gcf, 'Position', [100 100 600 240])
    tight_subplot(1,3,1,1, gutteraa, marginsaa)
    plot( squeeze(binz_posE_rt(si,ci,:)), squeeze(prop_cw_rt(si,ci,:)),'-', 'Color', colorz(ci,:), 'Linewidth',2); hold on; %'MarkerEdgeColor', colorz(ci,:)
    ylim([0 1])
    box off
    if ci == 1
        xlabel('stim strength')
        ylabel('prop resp CW')
    end
    %}
end
%end


figure
set(gcf, 'Position', [100 100 700 240])
tight_subplot(1,3,1,1, gutteraa, marginsaa)

for  ci = 1:Ncond
    
    he(ci) = fill([nanmean(squeeze(binz_pos(indi_sel,ci,:)),1) nanmean(squeeze(binz_pos(indi_sel,ci,end:-1:1)),1)],   ...
        [nanmean(squeeze(prop_cw_pr_all(indi_sel,ci,:)),1)-nanstd(squeeze(prop_cw_pr_all(indi_sel,ci,:)),1)./sqrt(length(indi_sel))...
        fliplr(nanmean(squeeze(prop_cw_pr_all(indi_sel,ci,:)),1) + nanstd(squeeze(prop_cw_pr_all(indi_sel,ci,:)),1)./sqrt(length(indi_sel))) ],colorz_shade(ci,:), 'EdgeColor', 'None'); hold on;
    plot(nanmean(squeeze(binz_pos(indi_sel,ci,:)),1), nanmean(squeeze(prop_cw_all(indi_sel,ci,:)),1),'o','MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:) ); hold on;
    errorbar(nanmean(squeeze(binz_pos(indi_sel,ci,:)),1), nanmean(squeeze(prop_cw_all(indi_sel,ci,:)),1), nanstd(squeeze(prop_cw_all(indi_sel,ci,:)),1)/sqrt(length(indi_sel)), 'Color',  colorz(ci,:),'Linestyle', 'none','Linewidth',1.5,'CapSize',0); hold on;
    
end
box off
xlabel('stim strength')
ylabel('prop resp CW')

%prop_cw_high_conf

tight_subplot(1,3,1,2, gutteraa, marginsaa)
for  ci = 1:Ncond
    
    he(ci) = fill([nanmean(squeeze(binz_pos(indi_sel,ci,:)),1) nanmean(squeeze(binz_pos(indi_sel,ci,end:-1:1)),1)],   ...
        [nanmean(squeeze(prop_cw_high_conf_pred(indi_sel,ci,:)),1)-nanstd(squeeze(prop_cw_high_conf_pred(indi_sel,ci,:)),1)./sqrt(length(indi_sel))...
        fliplr(nanmean(squeeze(prop_cw_high_conf_pred(indi_sel,ci,:)),1) + nanstd(squeeze(prop_cw_high_conf_pred(indi_sel,ci,:)),1)./sqrt(length(indi_sel))) ],colorz_shade(ci,:), 'EdgeColor', 'None'); hold on;
    plot(nanmean(squeeze(binz_pos(indi_sel,ci,:)),1), nanmean(squeeze(prop_cw_high_conf(indi_sel,ci,:)),1),'o','MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:) ); hold on;
    errorbar(nanmean(squeeze(binz_pos(indi_sel,ci,:)),1), nanmean(squeeze(prop_cw_high_conf(indi_sel,ci,:)),1), nanstd(squeeze(prop_cw_high_conf(indi_sel,ci,:)),1)/sqrt(length(indi_sel)), 'Color',  colorz(ci,:),'Linestyle', 'none','Linewidth',1.5,'CapSize',0); hold on;
    
end
box off
xlabel('stim strength')
ylabel('prop resp CW')
set(gca, 'tickdir', 'out')
%title('stims CW')
title('High conf')

tight_subplot(1,3,1,3, gutteraa, marginsaa)
for  ci = 1:Ncond
    
    he(ci) = fill([nanmean(squeeze(binz_pos(indi_sel,ci,:)),1) nanmean(squeeze(binz_pos(indi_sel,ci,end:-1:1)),1)],   ...
        [nanmean(squeeze(prop_cw_low_conf_pred(indi_sel,ci,:)),1)-nanstd(squeeze(prop_cw_low_conf_pred(indi_sel,ci,:)),1)./sqrt(length(indi_sel))...
        fliplr(nanmean(squeeze(prop_cw_low_conf_pred(indi_sel,ci,:)),1) + nanstd(squeeze(prop_cw_low_conf_pred(indi_sel,ci,:)),1)./sqrt(length(indi_sel))) ],colorz_shade(ci,:), 'EdgeColor', 'None'); hold on;
    plot(nanmean(squeeze(binz_pos(indi_sel,ci,:)),1), nanmean(squeeze(prop_cw_low_conf(indi_sel,ci,:)),1),'o','MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:) ); hold on;
    errorbar(nanmean(squeeze(binz_pos(indi_sel,ci,:)),1), nanmean(squeeze(prop_cw_low_conf(indi_sel,ci,:)),1), nanstd(squeeze(prop_cw_low_conf(indi_sel,ci,:)),1)/sqrt(length(indi_sel)), 'Color',  colorz(ci,:),'Linestyle', 'none','Linewidth',1.5,'CapSize',0); hold on;
    
end
box off
xlabel('stim strength')
ylabel('prop resp CW')
set(gca, 'tickdir', 'out')
%title('stims CCW')
title('Low conf')
set(gca, 'tickdir', 'out')
psname = 'psych_curves_broken_by_high_vs_low_conf.pdf'
%print_pdf(psname)


diff_in_mu_likelihood = params_bm_allV(:,4,4)- params_bm_allV(:,3,4)
diff_in_starting_point = squeeze(params_ddm(2,:,4,4))-squeeze(params_ddm(2,:,3,4))

diff_in_sigma_encoding = params_bm_allV(:,4,2)- params_bm_allV(:,3,2)
diff_in_mean_drift_rate = squeeze(params_ddm(2,:,4,1))-squeeze(params_ddm(2,:,3,1))

[r,p]= corr(diff_in_starting_point', diff_in_mu_likelihood, 'type', 'Spearman')

[r,p]= corr(diff_in_mean_drift_rate', diff_in_sigma_encoding, 'type', 'Spearman')


[r,p]= corr(diff_in_starting_point', compensation, 'type', 'Spearman')

%}