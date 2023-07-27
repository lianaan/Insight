clear all; close all;

colorz = [41 138 8; 152 191 100; 0 0 255; 137 195 255 ]/255;
colorz_shade_mid = (colorz + 4*repmat([228 228 228], 4, 1)/255)/5;
colorz_shade = (colorz + 3*repmat([228 228 228], 4, 1)/255)/4;

marginsa = [0.08 0.03 0.06 0.1]; %MARGINS = [LEFT RIGHT BOTTOM TOP]
guttera = [0.03 0.06];


load('dv_val_mat.mat')
load('stims_set_load.mat')
load('alldata_E2.mat')
Nsubj = length(alldata);

stims = stims_set_load; % always same order

Ncond = 4;
Ntrials = 121;
nbinz = 11;

binz = [];
for j = 1:nbinz
    binz(j) = quantile(stims, j/nbinz);
end

binz = [min(stims)*1.001 binz ];
binz_posE = (binz(2:end)+binz(1:end-1))/2;

stim_off_ts = 1035+3036+1035; 
stim_on_ts = 1035+3036; 
adaptor_on_ts = 1035;


window_len = 500;
window_start_post = 2000;
window_len_sliding = 200;
time_step = 20; 

load('pupil_all_tr_cond_subj_ALL.mat')
pupil_all_tr_condE = NaN(Ncond, Ntrials, 12000);
for si = 1: Nsubj
    
    data = alldata(si,:);
    pupil_all_tr_cond = squeeze(pupil_all_tr_cond_subj_ALL(si,:,:,:));
    
    valaa = [squeeze(pupil_all_tr_cond(1,:,1:end)); squeeze(pupil_all_tr_cond(2,:,1:end));...
        squeeze(pupil_all_tr_cond(3,:,1:end)); squeeze(pupil_all_tr_cond(4,:,1:end))];
    nanstd_valaa = nanstd(valaa);
    nanstd_valaa(nanstd_valaa<0.000001) = NaN;
    
    maxx_all = max(max([squeeze(pupil_all_tr_cond(1,:,:)); squeeze(pupil_all_tr_cond(2,:,:));...
        squeeze(pupil_all_tr_cond(3,:,:)); squeeze(pupil_all_tr_cond(4,:,:))]));
    pupil_all_tr_cond_div_by_max = pupil_all_tr_cond./maxx_all;
    
    
    for ci = 1:Ncond
        for ti = 1:Ntrials
            pupil_all_tr_condE(ci,ti,1:1:end) = squeeze(pupil_all_tr_cond_div_by_max(ci,ti,:));
        end
        
        rt_all_tr(si,ci,:) = data(1,ci).resp_times;
        resp_all_tr(si,ci,:) = data(1,ci).resp;
        cf_all_tr(si,ci,:) = data(1,ci).conf;
    end
    pupil_all_tr_condE_all_sbj(si, 1:Ncond, 1:Ntrials, 1:12000)= pupil_all_tr_condE;
   
end




ts_length0 = 7200;
ts_length = 7200;
ts_bef = 4500;
ts_aft = 2000; 
ts_sell = ts_bef+ts_aft+1;

for si = 1:Nsubj
    
    data = alldata(si,:);
 
    for ci = 1:Ncond
        
        clear va_all; clear va_all_r_locked;
        va_all = squeeze(pupil_all_tr_condE_all_sbj(si,ci,:,1:ts_length));
        
        for tri = 1:Ntrials
            ts_beg = floor(1035+3036+rt_all_tr(si,ci,tri)*1000-ts_bef);
            ts_end = floor(1035+3036+rt_all_tr(si,ci,tri)*1000+ts_aft);
            if ts_end>12000
                va_all_r_locked(tri,1:1:(12000-ts_beg+1)) = squeeze(pupil_all_tr_condE_all_sbj(si,ci,tri,ts_beg:1:12000));
                va_all_r_locked(tri,(12000-ts_beg+1)+1:1:ts_sell) = NaN(1, ts_sell-(12000-ts_beg+1));
            else
                va_all_r_locked(tri, 1:(ts_bef+ts_aft+1)) = ...
                    squeeze(pupil_all_tr_condE_all_sbj(si,ci,tri,ts_beg:ts_end))';
            end
        end
        pupil_all_tr_condE_all_sbj_r_locked(si,ci,1:Ntrials,:) = va_all_r_locked;
        
    end
end
%%

%% actually do the regressions across every window ---takes time --- see below for the option to load
% the saved values
%{
for wi = 1: number_of_wi
    
    if rem(wi,10) == 0
        wi
    end
    window_interval = (adaptor_on_ts + 500+ time_step*(wi-1)):1:  (adaptor_on_ts + 500+ time_step*(wi-1)+window_len);
    pupil_valz = [];
    
    for ci = 1:Ncond
        pupil_valz(1:Nsubj,ci,1:Ntrials) = nanmedian(squeeze(pupil_all_tr_condE_all_sbj(1:Nsubj, ci, 1:Ntrials, window_interval)),3);
        if r_locked
            pupil_valz(1:Nsubj,ci,1:Ntrials)=  nanmedian(squeeze(pupil_all_tr_condE_all_sbj_r_locked(1:Nsubj, ci, 1:Ntrials, window_interval)),3);
        end
    end
    
    pupils_for_lme = [];
    cond_for_lme = [];
    stim_strength_for_lme = [];
    stim_strength_rel_for_lme = [];
    rt_for_lme = [];
    subject_for_lme = [];
    dec_var_for_lme = [];
    conf_for_lme = [];
    load('dv_val_mat.mat')
    for i = 1: Nsubj
        
        subject_for_lme = [subject_for_lme; i*ones(1, Ntrials*Ncond)'];
        pupils_for_lme = [pupils_for_lme; squeeze(pupil_valz(i,1,:)); squeeze(pupil_valz(i,2,:));...
            squeeze(pupil_valz(i,3,:)); squeeze(pupil_valz(i,4,:))];
        cond_for_lme = [cond_for_lme; [ones(1,121) 2*ones(1,121) 3*ones(1,121) 4*ones(1,121)]'];
        stim_strength_for_lme = [ stim_strength_for_lme; [repmat(abs(stims_set_load),1,4)]'];
        
        rt_for_lme = [rt_for_lme; squeeze(rt_all_tr(i,1,:)); squeeze(rt_all_tr(i,2,:)); squeeze(rt_all_tr(i,3,:)); squeeze(rt_all_tr(i,4,:))];
        dec_var_for_lme = [dec_var_for_lme; squeeze(abs(dv_val(i,1,:))); squeeze(abs(dv_val(i,2,:))); squeeze(abs(dv_val(i,3,:))); squeeze(abs(dv_val(i,4,:)))];
        conf_for_lme = [conf_for_lme; squeeze(cf_all_tr(i,1,:)); squeeze(cf_all_tr(i,2,:));...
            squeeze(cf_all_tr(i,3,:)); squeeze(cf_all_tr(i,4,:))];
        
    end
    
    [~, ~, stim_strength_for_lme_ranking] = unique(stim_strength_for_lme);
    [~, ~, dec_var_for_lme_ranking] = unique(dec_var_for_lme);
    [~, ~, pupils_for_lme_ranking] = unique(pupils_for_lme);
    [~, ~, rt_for_lme_ranking] = unique(rt_for_lme);
    
    tbl_abs_stim = table(stim_strength_for_lme_ranking,pupils_for_lme_ranking,cond_for_lme, subject_for_lme,'VariableNames',{'StimStrength','PupilArea','Condition', 'Subject'});
    tbl_abs_stim.Condition = nominal(tbl_abs_stim.Condition);
    tbl_abs_stim.Subject = nominal(tbl_abs_stim.Subject);
    
    tbl_abs_dv = table(dec_var_for_lme_ranking, pupils_for_lme_ranking,cond_for_lme, subject_for_lme,'VariableNames',{'DecVar','PupilArea','Condition', 'Subject'});
    tbl_abs_dv.Condition = nominal(tbl_abs_dv.Condition);
    tbl_abs_dv.Subject = nominal(tbl_abs_dv.Subject);
    
    lme_abs_stim = fitlme(tbl_abs_stim,'PupilArea ~ 1 + StimStrength+Condition + (1+StimStrength+Condition|Subject)');
    lme_abs_dv = fitlme(tbl_abs_dv,'PupilArea ~ 1 + DecVar+Condition + (1+DecVar+Condition|Subject)');
    
    lme_abs_stim_all{wi} = lme_abs_stim;
    lme_abs_dv_all{wi} = lme_abs_dv;
    
    
    R2_stim(wi) = lme_abs_stim.Rsquared.adjusted;
    AIC_stim(wi) = lme_abs_stim.ModelCriterion.AIC;
    BIC_stim(wi) = lme_abs_stim.ModelCriterion.BIC;
    
    t_stim(wi) = lme_abs_stim.Coefficients.tStat(2);
    p_stim(wi) = lme_abs_stim.Coefficients.pValue(2);
    
    t_stim_cond2(wi) = lme_abs_stim.Coefficients.tStat(3);
    p_stim_cond2(wi) = lme_abs_stim.Coefficients.pValue(3);
    
    t_stim_cond3(wi) = lme_abs_stim.Coefficients.tStat(4);
    p_stim_cond3(wi) = lme_abs_stim.Coefficients.pValue(4);
    
    t_stim_cond4(wi) = lme_abs_stim.Coefficients.tStat(5);
    p_stim_cond4(wi) = lme_abs_stim.Coefficients.pValue(5);
    
    % .....
    R2_dv(wi) = lme_abs_dv.Rsquared.adjusted;
    AIC_dv(wi) = lme_abs_dv.ModelCriterion.AIC;
    BIC_dv(wi) = lme_abs_dv.ModelCriterion.BIC;
    
    t_dv(wi) = lme_abs_dv.Coefficients.tStat(2);
    p_dv(wi) = lme_abs_dv.Coefficients.pValue(2);
    
    t_dv_cond2(wi) = lme_abs_dv.Coefficients.tStat(3);
    p_dv_cond2(wi) = lme_abs_dv.Coefficients.pValue(3);
    
    t_dv_cond3(wi) = lme_abs_dv.Coefficients.tStat(4);
    p_dv_cond3(wi) = lme_abs_dv.Coefficients.pValue(4);
    
    t_dv_cond4(wi) = lme_abs_dv.Coefficients.tStat(5);
    p_dv_cond4(wi) = lme_abs_dv.Coefficients.pValue(5);
    
    
    
end


if r_locked == 0
    %save ('GLMEs_sliding_window_stim_locked_.mat', 'lme_abs_stim_all','lme_abs_dv_all', '-mat')
elseif r_locked == 1
    %save ('GLMEs_sliding_window_resp_locked_.mat', 'lme_abs_stim_all','lme_abs_dv_all', '-mat')
end
%}

%%
r_locked = 0;


stim_on_ts = 1035+3036; 
adaptor_on_ts = 1035;
ts_end = 8003;


time_step = 20;%100;
window_len = 500;
window_start_post = 2000;
number_of_wi = (7200- adaptor_on_ts - 500)/time_step;
number_of_wi = floor(number_of_wi);
window_len_sliding = 200;


if r_locked == 0
    load('GLMEs_sliding_window_stim_locked_.mat')
elseif r_locked == 1
    load('GLMEs_sliding_window_resp_locked_.mat')
end
for wi = 1: number_of_wi
    lme_abs_stim = lme_abs_stim_all{wi};
    lme_abs_dv = lme_abs_dv_all{wi};
    
    
    R2_stim(wi) = lme_abs_stim.Rsquared.adjusted;
    AIC_stim(wi) = lme_abs_stim.ModelCriterion.AIC;
    BIC_stim(wi) = lme_abs_stim.ModelCriterion.BIC;
    
    t_stim(wi) = lme_abs_stim.Coefficients.tStat(2);
    p_stim(wi) = lme_abs_stim.Coefficients.pValue(2);
    
    t_stim_cond2(wi) = lme_abs_stim.Coefficients.tStat(3);
    p_stim_cond2(wi) = lme_abs_stim.Coefficients.pValue(3);
    
    t_stim_cond3(wi) = lme_abs_stim.Coefficients.tStat(4);
    p_stim_cond3(wi) = lme_abs_stim.Coefficients.pValue(4);
    
    t_stim_cond4(wi) = lme_abs_stim.Coefficients.tStat(5);
    p_stim_cond4(wi) = lme_abs_stim.Coefficients.pValue(5);
    
    % .....
    R2_dv(wi) = lme_abs_dv.Rsquared.adjusted;
    AIC_dv(wi) = lme_abs_dv.ModelCriterion.AIC;
    BIC_dv(wi) = lme_abs_dv.ModelCriterion.BIC;
    
    t_dv(wi) = lme_abs_dv.Coefficients.tStat(2);
    p_dv(wi) = lme_abs_dv.Coefficients.pValue(2);
    
    t_dv_cond2(wi) = lme_abs_dv.Coefficients.tStat(3);
    p_dv_cond2(wi) = lme_abs_dv.Coefficients.pValue(3);
    
    t_dv_cond3(wi) = lme_abs_dv.Coefficients.tStat(4);
    p_dv_cond3(wi) = lme_abs_dv.Coefficients.pValue(4);
    
    t_dv_cond4(wi) = lme_abs_dv.Coefficients.tStat(5);
    p_dv_cond4(wi) = lme_abs_dv.Coefficients.pValue(5);
    
end


figure

colorz_pink = [232 144 156]./255;
colorz_pink_shade  = (colorz_pink+ 3*[1 1 1])./4; 
linewi = 1.1;
window_len = 500;
window_len_sliding = 200;
time_step = 20;
fontsz = 12;

stim_locked = 1;
resp_locked = 0;

tlen1 = 0.014;
tlen2 = 0.014;

colorz = [41 138 8; 152 191 100; 0 0 255; 137 195 255 ]/255;
colorz_shade_mid = (colorz + 4*repmat([228 228 228], 4, 1)/255)/5;
colorz_shade = (colorz + 3*repmat([228 228 228], 4, 1)/255)/4;


color_red = [230 0 38]./255;
color_red_shade = (color_red + 7*repmat([254 254 254], 1, 1)/255)/8;
color_black = colorz(3,:);
color_black_shade = (color_black + 2*repmat([254 254 254], 1, 1)/255)/3;
color_purple = colorz(4,:);
color_purple_shade = (color_purple + 2*repmat([254 254 254], 1, 1)/255)/3;
p_thresh = 0.01;

%tight_subplot(4,7, 3,[1 2 3 4], guttera, marginsa)
subplot(2,1,1)
Npts_fill = 100;
x_fill1 = linspace(stim_on_ts+window_start_post, stim_on_ts+window_start_post+window_len, Npts_fill);
fill([x_fill1,x_fill1(end:-1:1)],  [-11*ones(1,Npts_fill) fliplr(42*ones(1,Npts_fill))], colorz_pink_shade, 'EdgeColor', 'None'); hold on;
plot(1:1:ts_end, zeros(1,ts_end), '-', 'Color', 'k'); hold on;
if stim_locked
    plot((1000)*ones(1,100), linspace(-200,1200,100), '--','Color', [0.5 0.5 0.5],'Linewidth',linewi ); hold on;
    plot((1000+3000)*ones(1,100), linspace(-200,1200,100), '--k','Linewidth',linewi ); hold on;
    plot((1000+3000+1000)*ones(1,100), linspace(-200,1200,100), '--k','Linewidth',linewi ); hold on;
elseif resp_locked
    plot((ts_bef+1)*ones(1,100), linspace(-0.1,40,100), '--k','Linewidth',linewi ); hold on;
end
box off


ms_size = 1;
for wi = 1 : number_of_wi-1
    window_interval = (adaptor_on_ts+ 500+ time_step*(wi-1)):1:  (adaptor_on_ts+ 500+ time_step*(wi-1)+window_len_sliding);
    window_interval_center = (window_interval(1)+ window_interval(end))/2;
    if p_stim(wi) < p_thresh
        lF1 = plot(window_interval_center, t_stim(wi), 'o', 'MarkerSize', ms_size, 'MarkerFaceColor', color_red, 'MarkerEdgeColor', color_red); hold on;
    else
        lF1_nonsig = plot(window_interval_center, t_stim(wi), 'o', 'MarkerSize', ms_size, 'MarkerFaceColor', color_red_shade, 'MarkerEdgeColor', color_red_shade); hold on;
    end
    
    if p_stim_cond3(wi) < p_thresh
        lF2 = plot(window_interval_center, t_stim_cond3(wi), 'o', 'MarkerSize', ms_size, 'MarkerFaceColor', color_black, 'MarkerEdgeColor', color_black); hold on;
    else
        lF2_nonsig=plot(window_interval_center, t_stim_cond3(wi), 'o', 'MarkerSize', ms_size, 'MarkerFaceColor', color_black_shade, 'MarkerEdgeColor', color_black_shade); hold on;
    end
    
    if p_stim_cond4(wi) < p_thresh
        lF3 = plot(window_interval_center, t_stim_cond4(wi), 'o', 'MarkerSize', ms_size, 'MarkerFaceColor', color_purple, 'MarkerEdgeColor', color_purple); hold on;
    else
        lF3_nonsig =  plot(window_interval_center, t_stim_cond4(wi), 'o', 'MarkerSize', ms_size, 'MarkerFaceColor', color_purple_shade, 'MarkerEdgeColor', color_purple_shade); hold on;
    end
    
end
xlim([0 7001])
ylim([-10 10])
set(gca, 'ytick', [-10:5:10])
set(gca, 'tickdir', 'out')
set(gca, 'FontSize', fontsz)
set(gca, 'ticklength',[tlen1 tlen2])
ylabel('T-statistic', 'FontName', 'Helvetica', 'FontSize', fontsz)
xlabel('time (ms)')
ll1 = legend([ lF1 lF2 lF3], {'|s|', 'Adapt-See', 'Adapt-Believe'})
legend boxoff
set(ll1, 'Location', 'NorthWest')

subplot(2,1,2)

x_fill1 = linspace(stim_on_ts+window_start_post, stim_on_ts+window_start_post+window_len, Npts_fill);
fill([x_fill1,x_fill1(end:-1:1)],  [-11*ones(1,Npts_fill) fliplr(42*ones(1,Npts_fill))], colorz_pink_shade, 'EdgeColor', 'None'); hold on;
plot(1:1:ts_end, zeros(1,ts_end), '-', 'Color', 'k'); hold on;
if stim_locked
    plot((1000)*ones(1,100), linspace(-200,1200,100), '--','Color', [0.5 0.5 0.5],'Linewidth',linewi ); hold on;
    plot((1000+3000)*ones(1,100), linspace(-200,1200,100), '--k','Linewidth',linewi  ); hold on;
    plot((1000+3000+1000)*ones(1,100), linspace(-200,1200,100), '--k','Linewidth',linewi  ); hold on;
elseif resp_locked
    plot((ts_bef+1)*ones(1,100), linspace(-0.1,40,100), '--k','Linewidth',linewi ); hold on;
end
box off

ms_size = 1;
for wi = 1 : number_of_wi-1
    window_interval = (adaptor_on_ts+ 500+ time_step*(wi-1)):1:  (adaptor_on_ts+ 500+ time_step*(wi-1)+window_len_sliding);
    window_interval_center = (window_interval(1)+ window_interval(end))/2;
    if p_dv(wi) < p_thresh
        lF1 = plot(window_interval_center, t_dv(wi), 'o', 'MarkerSize', ms_size, 'MarkerFaceColor', color_red, 'MarkerEdgeColor', color_red); hold on;
    else
        lF1_nonsig = plot(window_interval_center, t_dv(wi), 'o', 'MarkerSize', ms_size, 'MarkerFaceColor', color_red_shade, 'MarkerEdgeColor', color_red_shade); hold on;
    end
    
    if p_dv_cond3(wi) < p_thresh
        lF2 = plot(window_interval_center, t_dv_cond3(wi), 'o', 'MarkerSize', ms_size, 'MarkerFaceColor', color_black, 'MarkerEdgeColor', color_black); hold on;
    else
        lF2_nonsig=plot(window_interval_center, t_dv_cond3(wi), 'o', 'MarkerSize', ms_size, 'MarkerFaceColor', color_black_shade, 'MarkerEdgeColor', color_black_shade); hold on;
    end
    
    if p_dv_cond4(wi) < p_thresh
        lF3 = plot(window_interval_center, t_dv_cond4(wi), 'o', 'MarkerSize', ms_size, 'MarkerFaceColor', color_purple, 'MarkerEdgeColor', color_purple); hold on;
    else
        lF3_nonsig =  plot(window_interval_center, t_dv_cond4(wi), 'o', 'MarkerSize', ms_size, 'MarkerFaceColor', color_purple_shade, 'MarkerEdgeColor', color_purple_shade); hold on;
    end
    
end
xlim([0 7001]) 
ylim([-10 10])
set(gca, 'ytick', [-10:5:10])
set(gca, 'tickdir', 'out')
set(gca, 'FontSize', fontsz)
set(gca, 'ticklength',[tlen1 tlen2])
ylabel('T-statistic', 'FontName', 'Helvetica', 'FontSize', fontsz)
xlabel('time (ms)',  'FontName', 'Helvetica', 'FontSize', fontsz)
ll2 = legend([ lF1 lF2 lF3], {'|d|', 'Adapt-See', 'Adapt-Believe'})
legend boxoff
set(ll2, 'Location', 'NorthWest')

