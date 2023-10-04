clear all; close all;

%exp_i = 1;
exp_i = 2;
if exp_i == 1
    sbj_list = {'1','2','3','4','5','6','7','8','9','10', '11', '12','13','14','15','16','17','18','19','20','21','22'};
    adapt_type_all = [ 1 -1 -1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1];
elseif exp_i == 2
    sbj_list = {'1','2','3','4','5', '6','7'};
    
end


Nsubj = length(sbj_list);

curr_dir = pwd;
Ncond = 4;
Nfiles_per_cond = 2;
if exp_i == 1
    model_pred = 0;
    model_pred_bm = 1;
elseif exp_i == 2
    %mi = 1;  % select the model
    % mi = 1: Bayesian insight model
    % mi = 2: prior model
    % mi = 3: insight + prior model
    % mi = 4: response bias k_choice model
    mi = 5; % response bias + insight model % winning model for the control expt
    % mi = 6: response bias k_choice_and_confidence model
    model_pred = 0;
    model_pred_bm = 1;
end

load(['alldata_E2ctr.mat']);
%exclude subjects 2 and 3... 
%alldata = alldata([1 4:8],:); 
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
   
    load('psych_curves_fitting_m1_201_E2ctr_Nsubj_7.mat')
    params_psych_m2_all = NaN(Nsubj, Ncond,3);
    
    for ci = 1:Ncond
        params_psych_m2_all(:,ci,1) = mu_est_all(:,ci);
        params_psych_m2_all(:,ci,2) = sigma_est_all(:,ci);
        if ci == 1
            params_psych_m2_all(:,ci,3) = lambda_est_all(:,ci);%lambda_est_all;
        end
    end
    
    
    load(['params_all_models6_E2ctr_Nsubj_7.mat'])
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
    %in the control experiment, only 5 participants total
    sii = 5;
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
    title('Exp 2 ctr: group (N = 7)')
    
    
    
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
    si_vec = [ 2 3 5]; % pick 3 participants from exp 2 ctr
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
%print_pdf(psname)

%% MODEL FIGURE
%close all;

figure(5)
set(gcf, 'Position', [100 100 740 700])

indi_sel = 1:1:Nsubj;
xlim_min = -0.3;%min(min(min(binz_pos)));
xlim_max = 0.3;%max(max(max(binz_pos)));


tight_subplot(5,4,3,1,guttera, marginsa)
for  ci = 1:Ncond
    
    h(ci) = plot(mean(squeeze(binz_pos(indi_sel,ci,:)),1), nanmean(squeeze(rt_all(indi_sel,ci,:)),1),'o-','Color',colorz(ci,:),'MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:),'Linewidth', linewi); hold on;
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
        %mi = 1;
        for  ci = 1:Ncond
            if exp_i == 2 & model_pred_bm == 1
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
    if model_pred_bm
    if exp_i == 2 
        he(ci) = fill([mean(squeeze(binz_pos(indi_sel,ci,:)),1) mean(squeeze(binz_pos(indi_sel,ci,end:-1:1)),1)],   ...
            [mean(squeeze(prop_cw_pred_all(indi_sel,mi,ci,:)),1)-std(squeeze(prop_cw_pred_all(indi_sel,mi,ci,:)),1)./sqrt(length(indi_sel))...
            fliplr(mean(squeeze(prop_cw_pred_all(indi_sel,mi,ci,:)),1) + std(squeeze(prop_cw_pred_all(indi_sel,mi,ci,:)),1)./sqrt(length(indi_sel))) ],colorz_shade(ci,:), 'EdgeColor', 'None'); hold on;
    end
    plot(nanmean(squeeze(binz_pos(:,ci,:)),1), nanmean(squeeze(prop_cw_all(indi_sel,ci,:)),1),'o', 'Color',colorz(ci,:),'Linewidth', 0.1,'MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:),'Linewidth',0.01 ); hold on;
    errorbar(nanmean(squeeze(binz_pos(:,ci,:)),1), nanmean(squeeze(prop_cw_all(indi_sel,ci,:)),1), nanstd(squeeze(prop_cw_all(indi_sel,ci,:)),1)/sqrt(length(indi_sel)), 'Color',  colorz(ci,:), 'LineStyle', 'none','CapSize',0, 'Linewidth', linewi); hold on;
    else
       plot(nanmean(squeeze(binz_pos(:,ci,:)),1), nanmean(squeeze(prop_cw_all(indi_sel,ci,:)),1),'o-', 'Color',colorz(ci,:),'Linewidth', 0.1,'MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:),'Linewidth',linewi); hold on;
    errorbar(nanmean(squeeze(binz_pos(:,ci,:)),1), nanmean(squeeze(prop_cw_all(indi_sel,ci,:)),1), nanstd(squeeze(prop_cw_all(indi_sel,ci,:)),1)/sqrt(length(indi_sel)), 'Color',  colorz(ci,:), 'LineStyle', 'none','CapSize',0, 'Linewidth', linewi); hold on;
    end 
    
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
    if model_pred_bm
    if exp_i == 2
        he(ci) = fill([mean(squeeze(binz_pos(indi_sel,ci,:)),1) mean(squeeze(binz_pos(indi_sel,ci,end:-1:1)),1)],   ...
            [mean(squeeze(prop_cf_pred_all(indi_sel,mi,ci,:)),1)-std(squeeze(prop_cf_pred_all(indi_sel,mi,ci,:)), 1)./sqrt(length(indi_sel))...
            fliplr(mean(squeeze(prop_cf_pred_all(indi_sel,mi,ci,:)),1) + std(squeeze(prop_cf_pred_all(indi_sel,mi,ci,:)),1)./sqrt(length(indi_sel))) ],colorz_shade(ci,:), 'EdgeColor', 'None'); hold on;
    end
    plot(nanmean(squeeze(binz_pos(:,ci,:)),1), nanmean(squeeze(prop_cf_high_all(indi_sel,ci,:)),1),'o', 'Color',colorz(ci,:),'Linewidth', 1.1,'MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:) ); hold on;
    errorbar(nanmean(squeeze(binz_pos(:,ci,:)),1), nanmean(squeeze(prop_cf_high_all(indi_sel,ci,:)),1), nanstd(squeeze(prop_cf_high_all(indi_sel,ci,:)),1)/sqrt(length(indi_sel)), 'Color',  colorz(ci,:),'LineStyle','none','CapSize',0, 'Linewidth', linewi); hold on;
    else
        plot(nanmean(squeeze(binz_pos(:,ci,:)),1), nanmean(squeeze(prop_cf_high_all(indi_sel,ci,:)),1),'o-', 'Color',colorz(ci,:),'Linewidth', 1.1,'MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:) ); hold on;
    errorbar(nanmean(squeeze(binz_pos(:,ci,:)),1), nanmean(squeeze(prop_cf_high_all(indi_sel,ci,:)),1), nanstd(squeeze(prop_cf_high_all(indi_sel,ci,:)),1)/sqrt(length(indi_sel)), 'Color',  colorz(ci,:),'LineStyle','none','CapSize',0, 'Linewidth', linewi); hold on;
    end
    
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
psname = 'Model_fig_expt_controlEE.pdf'
%print_pdf(psname)

%%
%mu_likelihood
[p,h,stats]=signrank(squeeze(params_bm_all(indi_sel,3,5)),squeeze(params_bm_all(indi_sel,4,5)) )
% k_choice
[p,h,stats]=signrank(squeeze(params_bm_all(indi_sel,3,3)),squeeze(params_bm_all(indi_sel,4,3)) )


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

%% AIC and BIC medians and bootstrapped 95% CI
%rmi = 1;
rmi = 4;
nboot = 500000; 
for kk=1:nboot
        
        sampleeA=randsample(AIC_all(:,rmi),Nsubj,1);
        medians_AIC(kk) = median(sampleeA);
        
        sampleeB=randsample(BIC_all(:,rmi),Nsubj,1);
        medians_BIC(kk) = median(sampleeB);
        
end
   
ci_bnd_low = 0.025;
ci_bnd_high = 0.975;


medians_AIC_95CI(1:2) = [quantile(medians_AIC,ci_bnd_low); quantile(medians_AIC,ci_bnd_high)]
medians_BIC_95CI(1:2) = [quantile(medians_BIC,ci_bnd_low); quantile(medians_BIC,ci_bnd_high)]    

%%
[p,h,stats]=signrank((AIC_all(:,4)), (AIC_all(:,1)))
[p,h,stats]=signrank((AIC_all(:,5)), (AIC_all(:,1)), 'tail', 'left')

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
psname = ['Fig5_parts_Sept2023_exp_',num2str(exp_i) ,'CTR_Nsubj_',num2str(Nsubj),'.pdf'];
%print_pdf(psname)


