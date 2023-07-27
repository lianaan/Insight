clear all; close all;

colorz = [41 138 8; 152 191 100; 0 0 255; 137 195 255 ]/255;
colorz_shade_mid = (colorz + 4*repmat([228 228 228], 4, 1)/255)/5;
colorz_shade = (colorz + 3*repmat([228 228 228], 4, 1)/255)/4;

marginsa = [0.08 0.03 0.06 0.1]; %MARGINS = [LEFT RIGHT BOTTOM TOP]
guttera = [0.03 0.06];


load('params_all_models6_E2_Nsubj_22.mat')
load('alldata_E2.mat')
load('params_psych_curves_m1_m2_all.mat')
load('dv_val_mat.mat')
load('stims_set_load.mat')
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
stim_on_ts = 1035+3036; % RT is counted from stimulus on

adaptor_on_ts = 1035;

Nsubj = length(alldata);

window_len = 500;
window_start_post = 2000;

for sbji = 1: Nsubj
    
    data = alldata(sbji,:);
    
    if sbji<11
        load(['pupil_all_tr_cond_subj_', num2str(sbji),'.mat']);
    elseif sbji>=11
        load(['pupil_all_tr_cond_subj_', num2str(sbji+1),'.mat']);
    end
    
    pupil_all_tr_condE = pupil_all_tr_cond;
    
    valaa = [squeeze(pupil_all_tr_cond(1,:,1:end)); squeeze(pupil_all_tr_cond(2,:,1:end));...
        squeeze(pupil_all_tr_cond(3,:,1:end)); squeeze(pupil_all_tr_cond(4,:,1:end))];
    nanstd_valaa = nanstd(valaa);
    nanstd_valaa(nanstd_valaa<0.000001) = NaN;
    
    maxx_all = max(max([squeeze(pupil_all_tr_condE(1,:,:)); squeeze(pupil_all_tr_condE(2,:,:));...
        squeeze(pupil_all_tr_condE(3,:,:)); squeeze(pupil_all_tr_condE(4,:,:))]));
    
    pupil_all_tr_cond_div_by_max = pupil_all_tr_condE./maxx_all;
    
    for ci = 1:Ncond
        for ti = 1:Ntrials
            
            pupil_all_tr_condE(ci,ti,1:1:end) = squeeze(pupil_all_tr_cond_div_by_max(ci,ti,:));
            
            time_post_stim_1 = stim_on_ts+window_start_post;
            time_post_stim_1_all(sbji,ci,ti) = time_post_stim_1;
            if  time_post_stim_1+window_len-window_start_post<12000
                pupil_post_stim(sbji,ci,ti) = nanmedian(squeeze(pupil_all_tr_condE(ci,ti,time_post_stim_1:1:(time_post_stim_1+window_len))));
            else
                pupil_post_stim(sbji,ci,ti) = NaN;
            end
            
            time_pre_stim_1 = adaptor_on_ts+window_start_post;
            time_pre_stim_1_all(sbji,ci,ti) = time_pre_stim_1;
            if  time_pre_stim_1+window_len-window_start_post<12000
                pupil_pre_stim(sbji,ci,ti) = nanmedian(squeeze(pupil_all_tr_condE(ci,ti,time_pre_stim_1:1:(time_pre_stim_1+window_len))));
            else
                pupil_pre_stim(sbji,ci,ti) = NaN;
            end
            
            time_post_resp_1 =stim_on_ts+floor(1000*data(1,ci).resp_times(ti))+500;
            time_post_resp_1_all(sbji,ci,ti) = time_post_resp_1;
            if  time_post_resp_1+500<12000
                pupil_post_resp(sbji,ci,ti) = nanmedian(squeeze(pupil_all_tr_condE(ci,ti,time_post_resp_1:1:(time_post_resp_1+500))));
            else
                pupil_post_resp(sbji,ci,ti) = NaN;
            end
            
            time_pre_resp_1 =stim_on_ts+floor(1000*data(1,ci).resp_times(ti))-3500;
            time_pre_resp_1_all(sbji,ci,ti) = time_pre_resp_1;
            if  time_pre_resp_1 > adaptor_on_ts+1500 & time_pre_resp_1+window_len < stim_on_ts
                pupil_pre_resp(sbji,ci,ti) = nanmedian(squeeze(pupil_all_tr_condE(ci,ti,time_pre_resp_1:1:(time_pre_resp_1+window_len))));
            else
                pupil_pre_resp(sbji,ci,ti) = NaN;
            end
            
        end
        
        rt_all_tr(sbji,ci,:) = data(1,ci).resp_times;
        resp_all_tr(sbji,ci,:) = data(1,ci).resp;
        cf_all_tr(sbji,ci,:) = data(1,ci).conf;
    end
    pupil_all_tr_condE_all_sbj(sbji, 1:Ncond, 1:Ntrials, 1:12000)= pupil_all_tr_condE;
    
    paramz =  squeeze(params_psych_m2_all(sbji,:,:));
    
    for ci = 1:4
        
        bin_dist = binz(3)-binz(2);
        binz_rel = [paramz(ci,1) - 5.5* bin_dist+ [0:1:11]*bin_dist];
        binzz_saved_rel(sbji,ci,1:12) = binz_rel;
        for j = 1:(nbinz)
            indi = find(data(1,ci).stims>binz(j) & data(1,ci).stims<=binz(j+1) );
            if length(indi)>0
                pupil_all_post_stim(sbji,ci,j) = nanmedian(pupil_post_stim(sbji,ci,indi));
                pupil_all_pre_stim(sbji,ci,j) = nanmedian(pupil_pre_stim(sbji,ci,indi));
                
                pupil_all_post_resp(sbji,ci,j) = nanmedian(pupil_post_resp(sbji,ci,indi));
                pupil_all_pre_resp(sbji,ci,j) = nanmedian(pupil_pre_resp(sbji,ci,indi));
            else
                pupil_all_post_stim(sbji,ci,j) = nan;
                pupil_all_pre_stim(sbji,ci,j) = nan;
                pupil_all_post_resp(sbji,ci,j) = nan;
                pupil_all_pre_resp(sbji,ci,j) = nan;
            end
            rt_all(sbji,ci,j) = nanmedian(data(1,ci).resp_times(indi));
            
            indi_rel = find(data(1,ci).stims>binz_rel(j) & data(1,ci).stims<=binz_rel(j+1) );
            if length(indi_rel)>0
                pupil_all_rel(sbji,ci,j) = nanmedian(pupil_post_stim(sbji,ci,indi_rel));
            else
                pupil_all_rel(sbji,ci,j) = nan;
            end
        end
        
        
    end
    
end


%%

ts_length0 = 7200;
ts_length = 7200;
ts_bef = 4500;
ts_aft = 2000; ts_sell = ts_bef+ts_aft+1;

sigma_means = mean(params_psych_m2_all(:,:,2));
sigma_meansE = mean(reshape(params_psych_m2_all(:,:,2), Nsubj*Ncond,1));


for si = 1:Nsubj
    
    data = alldata(si,:);
    
    paramz =  squeeze(params_psych_m2_all(si,:,:));
    
    for ci = 1:Ncond
        
        indi_high_conf = data(1,ci).conf == 1;
        indi_low_conf = data(1,ci).conf == 0;
        
        number_of_low_conf_trials(si,ci) = sum(indi_low_conf);
        
        
        clear va_h; clear va_l;
        clear va_all; clear va_all_r_locked;
        va_h = squeeze(pupil_all_tr_condE_all_sbj(si,ci,indi_high_conf,1:ts_length));
        va_l = squeeze(pupil_all_tr_condE_all_sbj(si,ci,indi_low_conf,1:ts_length));
        
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
        
        if number_of_low_conf_trials(si,ci)>10 & number_of_low_conf_trials(si,ci)<111
            avg_pupil_trace_low_conf(si,ci,1:ts_length) = nanmedian(va_l);%nanmean(va_l);
            avg_pupil_trace_high_conf(si,ci,1:ts_length) = nanmedian(va_h);%nanmean(va_h);
            
        end
        avg_pupil_trace_all(si,ci,1:ts_length) = nanmedian(va_all);
        avg_pupil_trace_all_r_locked(si,ci,1:(ts_bef+ts_aft+1)) = nanmedian(va_all_r_locked);
        
        
    end
end

%% FIGURE 6 FINAL

colorz = [41 138 8; 152 191 100; 0 0 255; 137 195 255 ]/255;
colorz_shade_mid = (colorz + 3*repmat([254 254 254], 4, 1)/255)/4;
colorz_shade = (colorz + 8.5*repmat([254 254 254], 4, 1)/255)/10;

colorz_pink = [232 144 156]./255;
colorz_pink_shade  = (colorz_pink+ 3*[1 1 1])./4; %(colorz_pink+ 2*[1 1 1])./3;

colorz_yellow = [196 145 2]./255;
colorz_yellow_shade = (colorz_yellow+ 3*[1 1 1])./4; %(colorz_yellow+ 2*[1 1 1])./3;

window_sel = 500;

tlen1 = 0.014;
tlen2 = 0.014;
face_alpha_val = 0.4;
linewi = 1.1;
fontsz = 12;

time_sel = 1:1: ts_length;


h = figure(6)
set(gcf, 'Position', [100 100 800 650])
guttera = [0.08 0.09];
marginsa = [0.0800    0.1000    0.1400    0.1000];

tight_subplot(4,7, 1,[1 2 3 4], guttera, marginsa)
Npts_fill = 100;
x_fill1 = linspace(stim_on_ts+window_start_post, stim_on_ts+window_start_post+window_len, Npts_fill);
fill([x_fill1,x_fill1(end:-1:1)],  [0*ones(1,Npts_fill) fliplr(1*ones(1,Npts_fill))], colorz_pink_shade, 'EdgeColor', 'None'); hold on;


for ci = 1:4
    
    clear va_all;
    va_all = squeeze(avg_pupil_trace_all(:,ci,1:ts_length));
    
    fill([time_sel, time_sel(end:-1:1)],  [nanmean(va_all) - nanstd(va_all)./sqrt(size(va_all,1)), ...
        fliplr(nanmean(va_all) + nanstd(va_all)./sqrt(size(va_all,1))) ], colorz_shade(ci,:), 'EdgeColor', 'None'); hold on;
    plot(1:1:ts_length, nanmean(va_all,1),'Color', colorz(ci,:), 'Linewidth',2); hold on;
    set(gca, 'tickdir', 'out')
    xlim([0 ts_length])
    box off
    
end
plot(1:1:ts_end, zeros(1,ts_end), '-', 'Color', 'k'); hold on;
plot((1000)*ones(1,100), linspace(-200,1200,100), '--','Color', [0.5 0.5 0.5],'Linewidth',linewi ); hold on;
plot((1000+3000)*ones(1,100), linspace(-200,1200,100), '--k','Linewidth',linewi ); hold on;
plot((1000+3000+1000)*ones(1,100), linspace(-200,1200,100), '--k','Linewidth',linewi ); hold on;
box off
set(gca, 'tickdir', 'out')
set(gca, 'FontSize', fontsz)
set(gca, 'ticklength',[tlen1 tlen2])
xlim([0 7001])
ylim([-0.05 0.6])
ylabel('pupil area (a.u.)', 'FontName', 'Helvetica', 'FontSize', fontsz)
xlabel('time (ms)',  'FontName', 'Helvetica', 'FontSize', fontsz)

ts_len_r = 1:1:(ts_bef+ts_aft+1);

%%

tlen1 = 0.024;
tlen2 = 0.024;
tight_subplot(4,7, 1,[ 6 7], guttera, marginsa)
x_fill1 = linspace((ts_bef+500+1), (ts_bef+500+1+window_sel), Npts_fill);
fill([x_fill1,x_fill1(end:-1:1)],  [0*ones(1,Npts_fill) fliplr(1*ones(1,Npts_fill))], colorz_yellow_shade, 'EdgeColor', 'None'); hold on;
for ci = 1:4
    
    clear va_all_r;
    va_all_r = squeeze(avg_pupil_trace_all_r_locked(:,ci,1:ts_sell));
    
    fill([ts_len_r, ts_len_r(end:-1:1)],  [nanmean(va_all_r) - nanstd(va_all_r)./sqrt(size(va_all_r,1)), ...
        fliplr(nanmean(va_all_r) + nanstd(va_all_r)./sqrt(size(va_all_r,1))) ], colorz_shade(ci,:), 'EdgeColor', 'None'); hold on;
    
    plot(ts_len_r, nanmean(va_all_r,1),'Color', colorz(ci,:), 'Linewidth',2); hold on;
    set(gca, 'tickdir', 'out')
    xlim([3000 max(ts_len_r)-500])
    box off
end
plot(1:1:ts_sell, zeros(1,ts_sell), '-', 'Color', 'k'); hold on;
plot((ts_bef+1)*ones(1,100), linspace(-0.1,1,100), '--k','Linewidth',linewi ); hold on;
box off
set(gca, 'tickdir', 'out')
set(gca, 'FontSize', fontsz)
set(gca, 'ticklength',[tlen1 tlen2])
ylim([-0.05 0.6])
ylabel('pupil area (a.u.)', 'FontName', 'Helvetica', 'FontSize', fontsz)
xlabel('time (ms)',  'FontName', 'Helvetica', 'FontSize', fontsz)



%%

figure(6)
set(gcf, 'Position', [100 100 800 650])
time_sel_ci = (1035+0):1: ts_length;
time_sel_ci_1 = 1:1: length(time_sel_ci);
guttera_s = [0.03 0.02];
marginsa_s = marginsa;


eb_w = 0.4;%errorbar width
eb_t = 0.5; %errorbar line thickness
lw = 1.3; % line width
%msz_all = [5.5 5.5 5.5 5.5];
msz_all = [3.5 3.5 3.5 3.5];

y_max_bars = 0.38;

guttera2 = [0.02 0.07];
marginsa2 = marginsa;
grey_interval = (adaptor_on_ts+ window_start_post):1:(adaptor_on_ts+window_start_post+ window_len);

pink_interval = (stim_on_ts+ window_start_post):1:(stim_on_ts+ window_start_post+ window_len);
indi_sel = 1:1:22;


pupil_all_post_stim_SAVED = pupil_all_post_stim;
grand_average = nanmean(pupil_all_post_stim_SAVED(:));
for si = 1: Nsubj
    pupil_all_sj = squeeze(pupil_all_post_stim_SAVED(si,:,:));
    subj_average(si) = nanmean(pupil_all_sj(:));
    for ci = 1:Ncond
        
        pupil_all_post_stim(si,ci,:) = bsxfun(@minus, squeeze(pupil_all_post_stim_SAVED(si,ci,:))', subj_average(si)) + grand_average;
        
    end
end
%}
tight_subplot(4,9, 2,[4 5], [guttera_s(1) guttera2(2)], marginsa2)
x_fill1 = linspace(-0.4, 0.4, Npts_fill);
fill([x_fill1,x_fill1(end:-1:1)],  [-0.1*ones(1,Npts_fill) fliplr(0.5*ones(1,Npts_fill))], colorz_pink_shade, 'EdgeColor', 'None'); hold on;

for  ci = 1:Ncond
    
    h(ci)  = plot(binz_posE, nanmean(squeeze(pupil_all_post_stim(indi_sel,ci,:)),1),'o-','Color',colorz(ci,:), 'MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:), 'Linewidth',linewi ); hold on;
    errorbar(binz_posE, nanmean(squeeze(pupil_all_post_stim(indi_sel,ci,:)),1), nanstd(squeeze(pupil_all_post_stim(indi_sel,ci,:)),1)/sqrt(length(indi_sel)), 'Color',  colorz(ci,:),'Linestyle', 'none','Linewidth',1.5,'CapSize',0); hold on;
    
end
box off

xlim_min = -0.3;
xlim_max = 0.3;
ylim_min = 0;
ylim_max = 1.01;
xlim([xlim_min xlim_max])
plot(zeros(1,nbinz), linspace(0,1,nbinz), '--k', 'Linewidth', linewi); hold on;
box off
set(gca, 'tickdir', 'out')
set(gca, 'FontSize', fontsz)
set(gca, 'ticklength',[tlen1 tlen2])
set(gca, 'ytick', [0.0:0.1:0.5])
set(gca, 'xtick', [-0.3:0.1:0.3])
ylabel('pupil area (a.u.)', 'FontName', 'Helvetica', 'FontSize', fontsz)
xlabel('test stimulus speed clockwise (a.u.)',  'FontName', 'Helvetica', 'FontSize', fontsz)
ylim([0.00 0.25])
%ylim([-0.05 0.35])



% response-locked
figure(1)

eb_w = 0.4;%errorbar width
eb_t = 0.5; %errorbar line thickness
lw = 1; % line width
guttera2 = [0.02 0.07];
marginsa2 = marginsa;
yell_interval = (ts_bef+500+1):1: (ts_bef+500+1+500);


indi_sel = 1:1:22;
tight_subplot(4,9, 2,[8 9], [guttera_s(1) guttera2(2)], marginsa2)
x_fill1 = linspace(-0.4, 0.4, Npts_fill);
fill([x_fill1,x_fill1(end:-1:1)],  [0*ones(1,Npts_fill) fliplr(0.4*ones(1,Npts_fill))], colorz_yellow_shade, 'EdgeColor', 'None'); hold on;

pupil_all_post_resp_SAVED = pupil_all_post_resp;
grand_averageE = nanmean(pupil_all_post_resp_SAVED(:));
for si = 1: Nsubj
    pupil_all_sjE = squeeze(pupil_all_post_resp_SAVED(si,:,:));
    subj_averageE(si) = nanmean(pupil_all_sjE(:));
    for ci = 1:Ncond
        % adjustment of error bars as in Cousineau 2005
        
        %pupil_all_post_stim(si,ci,:) = bsxfun(@minus, squeeze(pupil_all_post_stim_SAVED(si,ci,:))', nanmean(squeeze(pupil_all_post_stim(si,ci,:)))) + nanmean(squeeze(pupil_all_post_stim(:,ci,:)));
        
        pupil_all_post_resp(si,ci,:) = bsxfun(@minus, squeeze(pupil_all_post_resp_SAVED(si,ci,:))', subj_averageE(si)) + grand_averageE;
       
    end
end

for  ci = 1:Ncond
    
    h(ci)  = plot(binz_posE, nanmean(squeeze(pupil_all_post_resp(indi_sel,ci,:)),1),'o-','Color',colorz(ci,:), 'MarkerSize', msz_all(ci),'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:), 'Linewidth',linewi ); hold on;
    errorbar(binz_posE, nanmean(squeeze(pupil_all_post_resp(indi_sel,ci,:)),1), nanstd(squeeze(pupil_all_post_resp(indi_sel,ci,:)),1)/sqrt(length(indi_sel)), 'Color',  colorz(ci,:),'Linestyle', 'none','Linewidth',1.5,'CapSize',0); hold on;
    
end
box off
xlim_min = -0.3;
xlim_max = 0.3;
ylim_min = 0;
ylim_max = 1.01;
xlim([xlim_min xlim_max])
plot(zeros(1,nbinz), linspace(0,1,nbinz), '--k', 'Linewidth', linewi); hold on;
box off
set(gca, 'tickdir', 'out')
set(gca, 'FontSize', fontsz)
set(gca, 'ticklength',[tlen1 tlen2])
ylabel('pupil area (a.u.)', 'FontName', 'Helvetica', 'FontSize', fontsz)
xlabel('test stimulus speed clockwise (a.u.)',  'FontName', 'Helvetica', 'FontSize', fontsz)
set(gca, 'ytick', [0:0.1:0.4])
set(gca, 'xtick', [-0.3:0.1:0.3])
ylim([0 0.25])


%%
psname = 'pupil_fig_new_error_barsE.pdf'
%print_pdf(psname)

%% two way repeated measures anova
clear vals; clear group; clear y;
vals = [];

for ji = 1:nbinz
    %vals = [vals; squeeze(pupil_all_post_stim_SAVED(:,3,ji)); squeeze(pupil_all_post_stim_SAVED(:,4,ji))];
    vals = [vals; squeeze(pupil_all_post_resp_SAVED(:,3,ji)); squeeze(pupil_all_post_resp_SAVED(:,4,ji))];% Adapt-See vs Adapt-Believe
end

Subjects = repmat([1:Nsubj],1,2*nbinz);%[[1:Nsubj] [1:Nsubj]];
%evidence, then Adapt-See vs Adapt-Believe

group = [[ones(1,(2*Nsubj)) 2*ones(1,(2*Nsubj)) 3*ones(1,(2*Nsubj)) 4*ones(1,(2*Nsubj)) 5* ones(1,(2*Nsubj))...
    6*ones(1,(2*Nsubj)) 7*ones(1,(2*Nsubj)) 8*ones(1,(2*Nsubj)) 9*ones(1,(2*Nsubj)) 10* ones(1,(2*Nsubj)) 11* ones(1,(2*Nsubj)) ]' ...
    [repmat([ones(1,Nsubj) 2*ones(1,Nsubj)],1,nbinz)]'];

y = vals;

g1 = group(:,1);
g2 = group(:,2);

[p,table,stats] = anovan(y,{g1,g2, Subjects}, 'random',3,'varnames', {'Evidence', 'Adapt_See_Vs_Adapt_Bel', 'Subject'});
[p_full,table_full,stats_full] = anovan(y,{g1,g2, Subjects}, 'random',3,'varnames', {'Evidence',  'Adapt_See_Vs_Adapt_Bel', 'Subject'}, 'model', 'full');

[p,table,stats] = anovan(y,{g1,g2}, 'random',2,'varnames', {'Evidence', 'Adapt_See_Vs_Adapt_Bel'}, 'model', 'full');
%%

varNames = cell(2*nbinz,1);
for i = 1 : 2*nbinz
    v = strcat('V',num2str(i));
    varNames{i,1} = v;
end

c1 = cell(1,11);
c1{1} = 'A'; c1{2} = 'B'; c1{3} = 'C'; c1{4} = 'D';  c1{5} = 'E'; c1{6} = 'F'; c1{7} = 'G'; c1{8} = 'H';
c1{9} = 'I'; c1{10} = 'J'; c1{11} = 'K';

for g1i = 1:nbinz
    g1_ra([2*g1i-1],1) = c1{g1i};
    g1_ra([2*g1i],1) = c1{g1i};
end
c2 = cell(1,2);
c2{1} = 'S'; c2{2} = 'B';
for g2i = 1:nbinz
    g2_ra([2*g2i-1],1)  = c2{1};
    g2_ra([2*g2i],1)  = c2{2};
end
y_ranova = [];
for vi = 1:22
    y_ranova = [y_ranova vals([22*(vi-1)+1:22*(vi-1)+22])]
end

pupil_area = array2table(y_ranova, 'VariableNames',varNames);
factorNames = {'Evidence', 'Adapt_See_Vs_Adapt_Bel'};
within = array2table([g1_ra g2_ra], 'VariableNames', factorNames);



rm = fitrm(pupil_area,'V1-V22~1','WithinDesign',within);
[ranovatblb] = ranova(rm, 'WithinModel','Evidence*Adapt_See_Vs_Adapt_Bel');

%Mrm1 = multcompare(rm,'Evidence','By','Adapt_See_Vs_Adapt_Bel','ComparisonType','dunn-sidak');


withins2 = within;
withins2.Evidence = categorical(withins2.Evidence);
withins2.Adapt_See_Vs_Adapt_Bel = categorical(withins2.Adapt_See_Vs_Adapt_Bel);


psname = 'avg_traces_across_all_version_sigma_min_set_0035_diff_sigma_max_set_012_EE.pdf'
%print_pdf(psname)

%}
%}

%%

load('GLMEs_sliding_window_stim_locked_.mat')

number_of_wi = length( lme_abs_stim_all);
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

%%

figure(6)
set(gcf, 'Position', [100 100 800 650])

window_len_sliding = 200;
time_step = 20;

stim_locked = 1;
resp_locked = 0;

tlen1 = 0.014;
tlen2 = 0.014;

color_red = [230 0 38]./255;
color_red_shade = (color_red + 7*repmat([254 254 254], 1, 1)/255)/8;
color_black = colorz(3,:);
color_black_shade = (color_black + 2*repmat([254 254 254], 1, 1)/255)/3;
color_purple = colorz(4,:);
color_purple_shade = (color_purple + 2*repmat([254 254 254], 1, 1)/255)/3;
p_thresh = 0.01;

tight_subplot(4,7, 3,[1 2 3 4], guttera, marginsa)
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


tight_subplot(4,7, 4,[1 2 3 4], guttera, marginsa)
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


%%

load('GLMEs_sliding_window_resp_locked_.mat')

number_of_wi = length( lme_abs_stim_all);
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

window_len_sliding = 200;
time_step = 20;


%%
stim_locked = 0;
resp_locked = 1;

color_red = [230 0 38]./255;
color_red_shade = (color_red + 7*repmat([254 254 254], 1, 1)/255)/8;
color_black = colorz(3,:);
color_black_shade = (color_black + 2*repmat([254 254 254], 1, 1)/255)/3;
color_purple = colorz(4,:);
color_purple_shade = (color_purple + 2*repmat([254 254 254], 1, 1)/255)/3;


tlen1 = 0.024;
tlen2 = 0.024;

tight_subplot(4,7, 3,[6 7], guttera, marginsa)
x_fill1 = linspace((ts_bef+500+1),  (ts_bef+500+1+window_sel), Npts_fill);
fill([x_fill1,x_fill1(end:-1:1)],  [-11*ones(1,Npts_fill) fliplr(42*ones(1,Npts_fill))], colorz_yellow_shade, 'EdgeColor', 'None'); hold on;
plot(1:1:ts_end, zeros(1,ts_end), '-', 'Color', 'k'); hold on;
if stim_locked
    plot((1000)*ones(1,100), linspace(-200,1200,100), '--','Color', [0.5 0.5 0.5],'Linewidth',2 ); hold on;
    plot((1000+3000)*ones(1,100), linspace(-200,1200,100), '--k','Linewidth',2 ); hold on;
    plot((1000+3000+1000)*ones(1,100), linspace(-200,1200,100), '--k','Linewidth',linewi ); hold on;
elseif resp_locked
    plot((ts_bef+1)*ones(1,100), linspace(-10,40,100), '--k','Linewidth',linewi ); hold on;
end
box off

ms_size = 1.5;
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
        lF2_nonsig=plot(window_interval_center,t_stim_cond3(wi), 'o', 'MarkerSize', ms_size, 'MarkerFaceColor', color_black_shade, 'MarkerEdgeColor', color_black_shade); hold on;
    end
    
    if p_stim_cond4(wi) < p_thresh
        lF3 = plot(window_interval_center, t_stim_cond4(wi), 'o', 'MarkerSize', ms_size, 'MarkerFaceColor', color_purple, 'MarkerEdgeColor', color_purple); hold on;
    else
        lF3_nonsig =  plot(window_interval_center, t_stim_cond4(wi), 'o', 'MarkerSize', ms_size, 'MarkerFaceColor', color_purple_shade, 'MarkerEdgeColor', color_purple_shade); hold on;
    end
    
end
xlim([3000 max(ts_len_r)-500])
set(gca, 'xtick', [3000:1000:6000])
ylim([-10 10])
set(gca, 'ytick', [-10:5:10])
set(gca, 'tickdir', 'out')
set(gca, 'FontSize', fontsz)
set(gca, 'ticklength',[tlen1 tlen2])
ylabel('T-statistic', 'FontName', 'Helvetica', 'FontSize', fontsz)


tight_subplot(4,7, 4,[6 7], guttera, marginsa)
x_fill1 = linspace((ts_bef+500+1),  (ts_bef+500+1+window_sel), Npts_fill);
fill([x_fill1,x_fill1(end:-1:1)],  [-11*ones(1,Npts_fill) fliplr(42*ones(1,Npts_fill))], colorz_yellow_shade, 'EdgeColor', 'None'); hold on;
plot(1:1:ts_end, zeros(1,ts_end), '-', 'Color', 'k'); hold on;
if stim_locked
    plot((1000+3000)*ones(1,100), linspace(-200,1200,100), '--k','Linewidth',2 ); hold on;
    plot((1000+3000+1000)*ones(1,100), linspace(-200,1200,100), '--k','Linewidth',2 ); hold on;
elseif resp_locked
    
    plot((ts_bef+1)*ones(1,100), linspace(-10,40,100), '--k','Linewidth',linewi ); hold on;
    
end
box off

ms_size = 1.5;
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
        lF2_nonsig=plot(window_interval_center,t_dv_cond3(wi), 'o', 'MarkerSize', ms_size, 'MarkerFaceColor', color_black_shade, 'MarkerEdgeColor', color_black_shade); hold on;
    end
    
    if p_dv_cond4(wi) < p_thresh
        lF3 = plot(window_interval_center, t_dv_cond4(wi), 'o', 'MarkerSize', ms_size, 'MarkerFaceColor', color_purple, 'MarkerEdgeColor', color_purple); hold on;
    else
        lF3_nonsig =  plot(window_interval_center, t_dv_cond4(wi), 'o', 'MarkerSize', ms_size, 'MarkerFaceColor', color_purple_shade, 'MarkerEdgeColor', color_purple_shade); hold on;
    end
    
end
ylabel('T-statistic', 'FontName', 'Helvetica', 'FontSize', fontsz)
xlabel('time (ms)',  'FontName', 'Helvetica', 'FontSize', fontsz)
xlim([3000 max(ts_len_r)-500])
set(gca, 'xtick', [3000:1000:6000])
ylim([-10 10])
set(gca, 'ytick', [-10:5:10])
set(gca, 'tickdir', 'out')
set(gca, 'FontSize', fontsz)
set(gca, 'ticklength', [tlen1 tlen2])

%%
psname = 'Fig6_pupil_ALL.pdf';


%% Pupil models
%fitlme
% pupil_post_stim ---> Nsubj*4*121
%pupil_valz = pupil_post_stim;
pupil_valz = pupil_post_resp;


pupils_for_lme = [];
cond_for_lme = [];
stim_strength_for_lme = [];
rt_for_lme = [];
subject_for_lme = [];
dec_var_for_lme = [];
conf_for_lme = [];
%load('dv_val_mat.mat')
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

tbl_abs_rt = table(rt_for_lme_ranking, pupils_for_lme_ranking,cond_for_lme, subject_for_lme,'VariableNames',{'RT','PupilArea','Condition', 'Subject'});
tbl_abs_rt.Condition = nominal(tbl_abs_rt.Condition);
tbl_abs_rt.Subject = nominal(tbl_abs_rt.Subject);

tbl_abs_rt_dv = table(rt_for_lme_ranking, dec_var_for_lme_ranking, pupils_for_lme_ranking,cond_for_lme, subject_for_lme,'VariableNames',{'RT', 'DecVar','PupilArea','Condition', 'Subject'});
tbl_abs_rt_dv.Condition = nominal(tbl_abs_rt_dv.Condition);
tbl_abs_rt_dv.Subject = nominal(tbl_abs_rt_dv.Subject);


lme_abs_stim = fitlme(tbl_abs_stim,'PupilArea ~ 1 + StimStrength+Condition + (1+StimStrength+Condition|Subject)');
lme_abs_dv = fitlme(tbl_abs_dv,'PupilArea ~ 1 + DecVar+Condition + (1+DecVar+Condition|Subject)');
lme_abs_rt = fitlme(tbl_abs_rt,'PupilArea ~ 1 + RT+Condition + (1+RT+Condition|Subject)');
lme_abs_rt_dv = fitlme(tbl_abs_rt_dv,'PupilArea ~ 1 + RT + DecVar + Condition + (1+RT+DecVar+Condition|Subject)');


%% only Adapt conditions

pupil_valz = pupil_post_resp;


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
    
    subject_for_lme = [subject_for_lme; i*ones(1, Ntrials*2)'];
    pupils_for_lme = [pupils_for_lme; squeeze(pupil_valz(i,3,:)); squeeze(pupil_valz(i,4,:))];
    cond_for_lme = [cond_for_lme; [3*ones(1,121) 4*ones(1,121)]'];
    stim_strength_for_lme = [ stim_strength_for_lme; [repmat(abs(stims_set_load),1,2)]'];
    
    rt_for_lme = [rt_for_lme; squeeze(rt_all_tr(i,3,:)); squeeze(rt_all_tr(i,4,:))];
    dec_var_for_lme = [dec_var_for_lme; squeeze(abs(dv_val(i,3,:))); squeeze(abs(dv_val(i,4,:)))];
    conf_for_lme = [conf_for_lme; squeeze(cf_all_tr(i,3,:)); squeeze(cf_all_tr(i,4,:))];
    
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

tbl_abs_rt = table(rt_for_lme_ranking, pupils_for_lme_ranking,cond_for_lme, subject_for_lme,'VariableNames',{'RT','PupilArea','Condition', 'Subject'});
tbl_abs_rt.Condition = nominal(tbl_abs_rt.Condition);
tbl_abs_rt.Subject = nominal(tbl_abs_rt.Subject);

tbl_abs_rt_dv = table(rt_for_lme_ranking, dec_var_for_lme_ranking, pupils_for_lme_ranking,cond_for_lme, subject_for_lme,'VariableNames',{'RT', 'DecVar','PupilArea','Condition', 'Subject'});
tbl_abs_rt_dv.Condition = nominal(tbl_abs_rt_dv.Condition);
tbl_abs_rt_dv.Subject = nominal(tbl_abs_rt_dv.Subject);



lme_abs_stim = fitlme(tbl_abs_stim,'PupilArea ~ 1 + StimStrength+Condition + (1+StimStrength+Condition|Subject)');
lme_abs_dv = fitlme(tbl_abs_dv,'PupilArea ~ 1 + DecVar+Condition + (1+DecVar+Condition|Subject)');
lme_abs_rt = fitlme(tbl_abs_rt,'PupilArea ~ 1 + RT+Condition + (1+RT+Condition|Subject)');
lme_abs_rt_dv = fitlme(tbl_abs_rt_dv,'PupilArea ~ 1 + RT + DecVar + Condition + (1+RT+DecVar+Condition|Subject)');


