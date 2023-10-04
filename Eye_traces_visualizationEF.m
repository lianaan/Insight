clear all; close all; clc

sbjind = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23]; 

preprocess_flag = 0; % as in Kret and Skaj-Shie, 2019; mixed with Filipowicz 2020

cond_orders = [[3 4 2 1];
    [ 3 4 2 1] ; 
    [ 3 4 2 1];
    [3 4 2 1]; 
    [3 4 1 2];
    [3 4 1 2]; 
    [3 4 2 1]; 
    [3 4 1 2];
    [3 4 2 1]; 
    [3 4 2 1];
    [3 4 1 2]; 
    [3 4 1 2]; 
    [3 4 1 2];  
    [3 4 2 1]  
    [3 4 1 2] 
    [3 4 2 1] 
    [3 4 2 1] 
    [3 4 2 1] 
    [3 4 1 2] 
    [3 4 1 2] 
    [3 4 2 1]  
    [3 4 2 1] 
    [3 4 2 1]]; 
load('stims_set_load.mat');
x_task_period_all_tr_all_blocks_all_subj = NaN(22,8, 61, 4, 9500);
y_task_period_all_tr_all_blocks_all_subj = NaN(22,8, 61, 4, 9500);


for si = 1:23
    si
    if si == 11
        continue
    end
    load(['alldata_eye_S',num2str(si), '.mat'])
    
    if si<11
        sii = si;
    elseif si >11
        sii = si-1;
    end
    resp_times = [alldata.resp_times];  % 484 trials
    conf_times = [alldata.conf_times];
    
    
    cii = cond_orders(si,:);
    block_cond = [cii(1) cii(1) cii(2) cii(2) cii(3) cii(3) cii(4) cii(4)];
    
    
    ms_rate_all = NaN(8,61, 960); % 8 blocks, 61 trials most (60 or 61) and 960 units of time as of now --- 960* 5 = 4800---consider extending later
    xpos_all_blocks = NaN(8, 450000);
    ypos_all_blocks = NaN(8, 450000);
    x_task_period_all_tr_all_blocks = NaN(8, 61, 4, 9500);
    y_task_period_all_tr_all_blocks = NaN(8, 61, 4, 9500);
    for block_no = 1:8
        
        screen_resolution = [1920 1080];
        screen_distance = 57; 
        screen_width = 59;
        screen_angle = 2*(180/pi)*(atan((screen_width/2) / screen_distance)) ; % total visual angle of screen
        screen_ppd = screen_resolution(1) / screen_angle;  % pixels per degree
        screen_fixposxy = screen_resolution .* [.5 .5];
        
        
        
        curr_dir = pwd;
        dirname = [curr_dir, '/', num2str(sbjind(si))]
        filez = dir([dirname]);
        flst = {filez.name};
        cd(dirname)
        all_data = dlmread([num2str(sbjind(si)),'Block', num2str(block_no), '_Z_A_B_C_D.txt']);
        pupil_dil = all_data(:,4);
        
        
        strtofind = ['block_no_',num2str(block_no),'.mat'];
        ix = regexp(flst,strtofind);
        ix =~ cellfun('isempty',ix);
        files = flst(ix);
        load(files{1})
        
        close all;
        
        figure
        set(gcf, 'Position', [100 100 500 400])
        marginsa = [0.08 0.08 0.06 0.1]; %MARGINS = [LEFT RIGHT BOTTOM TOP]
        guttera = [0.05 0.12];
        colorzz{1}=[210 105 30]/255;
        colorzz{2}=[107 142 35]/255;
        %1: StartTrial
        %2: Fix ON
        %3: Adaptor ON
        %4: Stim ON
        %5: stim off
        %6: Trial end
        ntrials_max = 500;
        times_pi = NaN(ntrials_max,6);
        clear task_per_times; clear task_per_times_cum;
        
        for pi = 2: 6 % task period index
            file_read = dlmread([num2str(sbjind(si)), 'Block',num2str(block_no), '_Z_A', num2str(pi),'TT.txt']);
            loc_in_times_pi = file_read(:,3); 
            times_pi(loc_in_times_pi, pi) = file_read(:,1);
        end
        
        times_pi(times_pi == 0)= NaN;
        
        task_per_times_raw = times_pi(:,3:end)- times_pi(:,2:end-1)+1;
        compl_trials = find(~isnan(task_per_times_raw(:,end)));
        task_per_times = task_per_times_raw(compl_trials,:);
        
        max_task_per_times = max(task_per_times); % Nan-pad up until this number
        min_task_per_times = min(task_per_times);
        
        task_per_times_cum(:,1) =  task_per_times(:,1);
        max_task_per_times_cum(1) = max_task_per_times(1);
        min_task_per_times_cum(1) = min_task_per_times(1);
        for i = 2: size(task_per_times,2)
            task_per_times_cum(:,i) = sum(task_per_times(:,1:i)')';
            max_task_per_times_cum(i) = sum(max_task_per_times(1:i));
            min_task_per_times_cum(i) = sum(min_task_per_times(1:i));
        end
  
        n_trialz = max(find(~isnan(times_pi(:,2))));
        n_trialz_end = find(~isnan(times_pi(:,6)));
        
        ts_lengths = [max_task_per_times 400]; % or min?
        
        xpos_all_tr = NaN(length(n_trialz_end),sum(max_task_per_times)); % min is 7202, max is 7207
        ypos_all_tr = NaN(length(n_trialz_end),sum(max_task_per_times));
        pupil_all_tr = NaN(length(n_trialz_end),sum(max_task_per_times));
        
        xpos_all_tr_cont = [];
        ypos_all_tr_cont = [];
        ts_cont_end = sum(min_task_per_times);%4800;
        x_task_period_all_tr = NaN(61, 4, max(max_task_per_times));
        y_task_period_all_tr = NaN(61, 4, max(max_task_per_times));
        
        for tti = 1:length(n_trialz_end)
            ti = n_trialz_end(tti);
            
            periods = find(~isnan(times_pi(ti,:))); 
            if length(periods) == 5
                length_periods = 4;
            else
                length_periods = length(periods);
            end
            
            for ppi = 1: length_periods
                min_min_b = NaN;
                min_min_e = NaN;
                [ts_begin_i,ts_begin_j] = find(all_data(:,1)== times_pi(ti,periods(ppi)));
                
                if isempty(ts_begin_i)
                    min_min_b = min(min(abs((all_data(:,1)- times_pi(ti,periods(ppi))))));
                    [ts_begin_i,ts_begin_j] = find(all_data(:,1)== times_pi(ti,periods(ppi))-min_min_b);
                    if isempty(ts_begin_i)
                        [ts_begin_i,ts_begin_j] = find(all_data(:,1)== times_pi(ti,periods(ppi))+min_min_b);
                    end
                end
                if ~isnan(times_pi(ti,periods(ppi)+1))
                    [ts_end_i,ts_end_j] = find(all_data(:,1)== times_pi(ti,periods(ppi)+1)-1);
                    
                    if isempty(ts_end_i)
                        min_min_e = min(min(abs((all_data(:,1)- times_pi(ti,periods(ppi)+1)))));
                        [ts_end_i,ts_end_j] = find(all_data(:,1)== times_pi(ti,periods(ppi)+1)-min_min_e);
                        if isempty(ts_end_i)
                            [ts_end_i,ts_end_j] = find(all_data(:,1)== times_pi(ti,periods(ppi)+1)+min_min_e);
                        end
                    end
                else
                    [ts_end_i,ts_end_j] = find(all_data == times_pi(ti,periods(ppi)+1)+ts_lengths(periods(ppi)-1));
                    if isempty(ts_end_i)
                        min_min_e = min(min(abs((all_data(:,1)- times_pi(ti,periods(ppi)+1)+ts_lengths(periods(ppi)-1)))));
                        [ts_end_i,ts_end_j] = find(all_data(:,1)== times_pi(ti,periods(ppi)+1)+ts_lengths(periods(ppi)-1)-min_min_e);
                        
                    end
                end
                
                
                length_of_eye_ts = ts_end_i- ts_begin_i +1;
                
                xpos = (all_data(ts_begin_i: ts_end_i, 2) - screen_fixposxy(1)); %/screen_ppd
                ypos = (all_data(ts_begin_i: ts_end_i, 3) - screen_fixposxy(2));
                pupil = pupil_dil(ts_begin_i: ts_end_i);
                
                xx_per_tr = [];
                if length_of_eye_ts >1000
                    xx_per_tr(1:length_of_eye_ts,1) = xpos;
                    xx_per_tr(1:length_of_eye_ts,2) = ypos;
                    
                    v_per_tr = vecvel(xx_per_tr,1000,2);
                    v_per_tr_mean(ppi,1:2) = nanmean(v_per_tr(1:end-1,:));
                    v_per_tr_sem(ppi,1:2) =  nanstd(v_per_tr(1:end-1,:))/sqrt(length(v_per_tr(1:end-1,:))-1);
                    v_per_tr_mean_across_tr(tti,ppi) = sqrt(v_per_tr_mean(ppi,1).^2+ v_per_tr_mean(ppi,2).^2);
                else
                    v_per_tr_mean_across_tr(tti,ppi) = NaN;
                end
                
                if ppi == 1
                    prevv = 1;
                else
                    prevv = max_task_per_times_cum(ppi-1);
                end
                
                
                xpos_all_tr(tti,prevv:prevv+ length(xpos)-1) = xpos';
                ypos_all_tr(tti,prevv:prevv+length(ypos)-1) = ypos; 
                pupil_all_tr(tti,prevv:prevv+length(pupil)-1) = pupil;
               
                x_task_period_all_tr(tti,ppi,1:length(xpos) ) = xpos./screen_ppd; 
                y_task_period_all_tr(tti,ppi,1:length(ypos) ) = ypos./screen_ppd;
                
            end
            if tti == 1
                xpos_all_tr_cont = xpos_all_tr(tti,1:ts_cont_end);
                ypos_all_tr_cont = ypos_all_tr(tti,1:ts_cont_end);
                pupil_all_tr_cont = pupil_all_tr(tti,1:ts_cont_end);
            elseif tti > 1
                xpos_all_tr_cont = [xpos_all_tr_cont xpos_all_tr(tti,1:ts_cont_end)-(-xpos_all_tr(tti-1,ts_cont_end)+xpos_all_tr(tti,1))+0.00001]; % changed on May 10 2023
                ypos_all_tr_cont = [ypos_all_tr_cont ypos_all_tr(tti,1:ts_cont_end)-(-ypos_all_tr(tti-1,ts_cont_end)+ypos_all_tr(tti,1))+0.00001];
                pupil_all_tr_cont = [pupil_all_tr_cont pupil_all_tr(tti,1:ts_cont_end)];
            end
            % concatenate trials   
        end
        
        % pre-process and remove baseline based on fixation
        %calculate this threshold based on the entire pupil time series
        %---edges that glue trials wont affect the values substantially
        % remove artifacts with high speed as in as in Kret and Skaj-Shie, 2019
        
        dilation_speed = max(abs(pupil_all_tr_cont(2:end-1)- pupil_all_tr_cont(1:end-2)), ...
            abs(pupil_all_tr_cont(3:end)- pupil_all_tr_cont(2:end-1)));
        MAD = nanmedian(abs(dilation_speed-nanmedian(dilation_speed)));
        if MAD <0.00001 % do not allow for MAD to be 0
            MAD = 1;
        end
        n_multl = 16;
        %Kret and Skaj-Shie, 2019 --they say case by case basis to decide
        threshold = nanmedian(dilation_speed)+ n_multl* MAD;
        
        task_lengthss = quantile(task_per_times_cum,0.9); %ignore really long trials with weird pupil artifacts
        task_length = task_lengthss(end); % will be a bit over 7000 ms
        
        for tti =  1:length(n_trialz_end)
            pupil_sel = pupil_all_tr(tti,1:end);
            pupil_dil = pupil_sel;
            dilation_speed_tr = max(abs(pupil_dil(2:end-1)- pupil_dil(1:end-2)), ...
                abs(pupil_dil(3:end)- pupil_dil(2:end-1)));
            prop_removed(tti) = sum(dilation_speed_tr>threshold)/ length(pupil_dil);
            pupil_dil(dilation_speed_tr>threshold) = NaN;
            %first identify gaps 
            to_interp = double(isnan(pupil_dil));
            t_pN = []; %vector of transitions from pupil to NaN
            t_Np = []; %vector of transitions from NaN to pupil
            t_pN = find([to_interp 0]==1 & [0 to_interp]==0);
            t_Np = find([to_interp 0]==0 & [0 to_interp]==1);
            lengths_of_NaN_stretches = t_Np - t_pN;
            
            if preprocess_flag == 1 & sum(lengths_of_NaN_stretches)<1000
                pupil_dil2 = pupil_dil;
                if t_pN(end) <  7100
                    end_of_interp = length(t_pN);
                elseif t_pN(end) >7100
                    end_of_interp = length(t_pN)-1;
                end
                
                for it = 1: end_of_interp
                    if t_pN(it)>20
                        pupil_dil2(t_pN(it):t_Np(it)-1) = interp1(t_pN(it)-20:1:t_pN(it)-2,pupil_dil(t_pN(it)-20:1:t_pN(it)-2),t_pN(it):1:t_Np(it)-1,'linear', 'extrap');
                    else
                        pupil_dil2(t_pN(it):t_Np(it)-1) = interp1(t_Np(it):1:t_Np(it)+20,pupil_dil(t_Np(it):1:t_Np(it)+20),t_pN(it):1:t_Np(it)-1,'linear', 'extrap');
                    end
                end
                
                pupil_dil4 = NaN(1,task_per_times_cum(tti,4));
                if t_Np(end) >7100
                    pupil_to_filt = pupil_dil2(1:t_pN(end)-1);
                else
                    pupil_to_filt = pupil_dil2;
                end
                
                % Urai 2017 Butterworth filter of order 2, bandpass between 0.01
                % and 10 Hz - 2 steps, one low pass and one high pass
                fs =   1000; % frequency of signal, 1000 Hz
                fcutlow = 0.01;   %low cut frequency in Hz
                fcuthigh = 10;   %high cut frequency in Hz
                % implementation of the bandpass filter, as a high, and
                % then a low
                
                % filter the pupil timecourse twice
                % first, get rid of slow drift
                [b,a]      = butter(2, fcutlow / fs, 'high');
                pupil_dil3 = filtfilt(b,a, pupil_to_filt); % filter with zero lag
                
                
                % also get rid of fast instrument noise
                [b,a] = butter(2, fcuthigh / fs, 'low');
                pupil_filtered = filtfilt(b,a, pupil_dil3);
                pupil_dil4(1:length(pupil_filtered)) = pupil_filtered;
                % did not downsample 
                
            else
                pupil_dil4 = pupil_dil;
            end
            % pick the last part of the fixation period to substract
            avgg = nanmean(pupil_dil4(task_per_times(tti,1)-400: task_per_times(tti,1)));
            pupilE = pupil_dil4 - avgg;
            pupil_all_trE(tti,1:length(pupil_dil4)) = pupilE;
        end
        %%
        %save(['pupil_all_tr_mean_',num2str(sbjind(si)),'_block_no_',num2str(block_no), '.mat'], 'pupil_all_trE','pupil_all_tr','task_length','max_task_per_times','min_task_per_times','-mat')
        %save(['vel_all_',num2str(sbjind(si)),'_block_no_',num2str(block_no),'.mat'], 'v_per_tr_mean_across_tr', '-mat')

        %% concatenated over all trials
        xpos_all_tr_cont = xpos_all_tr_cont/screen_ppd; % to convert to degreees
        ypos_all_tr_cont = ypos_all_tr_cont/screen_ppd;
        
        xpos_all_tr_contE = xpos_all_tr_cont; % to convert to degreees
        ypos_all_tr_contE = ypos_all_tr_cont;
        
        if ismember(si, [3 4 5 6 7]) 
            %interpolate where we have NaN
            to_interp_x = double(isnan(xpos_all_tr_cont));
            t_pN = []; %vector of transitions from pupil to NaN
            t_Np = []; %vector of transitions from NaN to pupil
            t_pN = find([to_interp_x 0]==1 & [0 to_interp_x]==0);
            t_Np = find([to_interp_x 0]==0 & [0 to_interp_x]==1);
            lengths_of_NaN_stretches_x = t_Np - t_pN;
            for txi = 1: length(t_pN)
                xpos_all_tr_contE(t_pN(txi):t_Np(txi)-1) =  xpos_all_tr_cont(t_pN(txi)-1);  %interp1( [t_pN(txi):t_Np(txi)-1], [xpos_all_tr_cont(t_pN(txi)-1): xpos_all_tr_cont(t_Np(txi))], t_pN(txi):t_Np(txi)-1 );
                ypos_all_tr_contE(t_pN(txi):t_Np(txi)-1) =  ypos_all_tr_cont(t_pN(txi)-1); %mean([ ypos_all_tr_cont(t_pN(txi)-1) ypos_all_tr_cont(t_Np(txi)) ]);
            end
            
        end
        % if participants 3,4,5,6,7 - where we did not use the heuristic filter extra, do Kalman filter
        filter2_size=7; %size of the 2nd Kalman filter in BMD reduced + threshold
        if ismember(si, [3 4 5 6 7])
            
            sigz = 0.01; % approx. based on BMD paper, Mihali et al, 2017
            sigx = 0.025; % aprox. based on BMD paper
            
            xka = [xpos_all_tr_contE' ypos_all_tr_contE'];
            
            
            T=size(xka,1);
            xf=zeros(T,2);
            xs=zeros(T,2);
            
            %Kalman filter 1
            Ez =sigz^2*eye(2);
            Ex =sigx^2*eye(2);
            P=(-Ez+sqrt(Ez.^2+4*Ez*Ex))/2;
            K=(P+Ez)*inv(P+Ez+Ex);
            
            xf(1,:)=xka(1,:);
            for t = 2:length(xka)
                xf(t,:)=(xf(t-1,:)'+K*(xka(t,:)-xf(t-1,:))')';
            end
            
            %backward part, Kalman smoother 1
            K=P*inv(P+Ez);
            xs(T,:)=xf(T,:);
            for t = (T-1):-1:1
                xs(t,:)=(xf(t,:)'+K*(xs(t+1,:)-xf(t,:))')';
            end
            
            %Kalman filter 2
            Ez=eye(2);
            Ex=filter2_size*eye(2);
            
            xf=zeros(T,2);
            P=(-Ez+sqrt(Ez.^2+4*Ez*Ex))/2;
            K=(P+Ez)*inv(P+Ez+Ex);
            
            xf(1,:)=xs(1,:);
            for t = 2:length(xka)
                xf(t,:)=xf(t-1,:)+(xs(t,:)-xf(t-1,:))*K;
            end
            %backward, Kalman smoother 2
            xs=zeros(T,2);
            Ks=P*inv(P+Ez);
            
            xs(T,:)=xf(T,:);
            for t = (T-1):-1:1
                xs(t,:)=xf(t,:)+(xs(t+1,:)-xf(t,:))*Ks;
            end
    
        end
        %%

        time_per_file = 20*ts_cont_end;
        count_files = floor( length(xpos_all_tr_cont)/time_per_file)+1;
        % maybe divide into some number of trials at a time
        
        cd ..
        root_dir = pwd;
        addpath(root_dir)
        xx_all = [];
 
        if ismember(si, [3 4 5 6 7])
            xx_all(:,1) = xs(:,1);
            xx_all(:,2) = xs(:,2);
        else
            xx_all(:,1) = xpos_all_tr_cont;
            xx_all(:,2) = ypos_all_tr_cont;
        end
        xpos_all_blocks(block_no,1: length(xpos_all_tr_cont)) = xpos_all_tr_cont;
        ypos_all_blocks(block_no,1: length(xpos_all_tr_cont)) = ypos_all_tr_cont;
        v = vecvel(xx_all,1000,2);
        clear sac;

        [sac, radius] = microsacc(xx_all,v,6,6); % changed to a min dur of 6 seconds for a microsaccade
        %[sac, radius] = microsacc(xx_all,v,6,12); % also tried a min dur of 12 seconds
        sac_all_blocks_all_subj{si, block_no} = sac;
   
  
        insta_vel = sqrt((xx_all(2:end,1)-xx_all(1:end-1,1)).^2+(xx_all(2:end,2)-xx_all(1:end-1,2)).^2)/0.001; % to convert ms to sec and get deg/sec
        % to compute the distribution of drift velocities, exclude the
        % periods the eye is in microsaccade state
        
        total_sac_times = [];
        for sis = 1: size(sac,1)
            total_sac_times = [total_sac_times sac(sis,1): sac(sis,2)];
            total_sac_times_excl = [total_sac_times sac(sis,1)-26: sac(sis,2)+26]; %cf Murakami et al 2005
        end
        
        insta_vel_exclude_microsacc =  insta_vel(setdiff(1:1:length(insta_vel),total_sac_times_excl));
        insta_vel_exclude_microsacc_exclude_artifacts = insta_vel_exclude_microsacc;
   
        insta_vel_exclude_microsacc_exclude_artifacts (insta_vel_exclude_microsacc>quantile(insta_vel_exclude_microsacc,0.95)) = NaN; % changed as of july 11
        insta_vel_microsacc =  insta_vel(total_sac_times);
        nanmean(insta_vel_exclude_microsacc_exclude_artifacts)
        nanstd(insta_vel_exclude_microsacc_exclude_artifacts)
        nanmedian(insta_vel_microsacc);%  19
        nanmean(insta_vel_microsacc);%  31.7251 
        
        sac_all_blocks_all_subj{sii, block_no} = sac;
        if size(sac,1) >0
            sac_size = sqrt(sac(:,6).^2+ sac(:,7).^2);
        else
            sac_size = NaN;
        end
        sac_small = sac(sac_size<1, :);
        sac_large = sac(sac_size>=1, :);
        sac_small_all_blocks_all_subj{sii, block_no} = sac_small;
        sac_large_all_blocks_all_subj{sii, block_no} = sac_large;
  
        sac_rate_all_blocks_all_subj(sii, block_no)= size(sac,1)/sum(~isnan(xpos_all_tr_cont)/1000);
        sac_small_rate_all_blocks_all_subj(sii, block_no) = size(sac_small,1)/sum(~isnan(xpos_all_tr_cont)/1000);
        sac_large_rate_all_blocks_all_subj(sii, block_no) = size(sac_large,1)/sum(~isnan(xpos_all_tr_cont)/1000);
        
     
        drift_vel_mean(sii,block_no) = nanmean(insta_vel_exclude_microsacc_exclude_artifacts);
        drift_vel_std(sii,block_no) = nanstd(insta_vel_exclude_microsacc_exclude_artifacts);
        microsacc_vel_mean(sii,block_no) = nanmean(insta_vel_microsacc);
        microsacc_vel_std(sii,block_no) = nanstd(insta_vel_microsacc);
        
        drift_vel_all(sii,block_no, 1: length(insta_vel_exclude_microsacc_exclude_artifacts)) = insta_vel_exclude_microsacc_exclude_artifacts;
        
        %% time course analysis
        trial_lengths = repmat(4800, [size(n_trialz_end,1) 1]);
        Ntr = size(n_trialz_end,1);
        ms_rate = nan(Ntr, ceil(max(trial_lengths)/5));
        t_passed = 0;
        for tr_ind = 1: Ntr 
            for tii = 1: floor(trial_lengths(tr_ind)/5)
                
                t_start = t_passed + 5*(tii-1) +1;
                t_end = t_passed+ 5*(tii-1)+51;
                
                if size(sac,1) >0
                    ms_start = find(sac(:,1)>t_start);
                    ms_end = find(sac(:,2)<t_end);
                    
                    ms_rate(tr_ind,tii) = length(intersect(ms_start, ms_end))/50*1000;
                else
                    ms_start = NaN;
                    ms_end = NaN;
                    ms_rate(tr_ind,tii) = NaN;
                end
                
            end
            t_passed = t_passed + trial_lengths(tr_ind);

        end
        ms_rate_all(block_no,1:size(ms_rate,1), 1:size(ms_rate,2)) = ms_rate;
        
        x_task_period_all_tr_all_blocks(block_no, 1:size(x_task_period_all_tr,1), 1:size(x_task_period_all_tr,2),...
            1:size(x_task_period_all_tr,3) ) = x_task_period_all_tr;
        y_task_period_all_tr_all_blocks(block_no, 1:size(y_task_period_all_tr,1), 1:size(y_task_period_all_tr,2),...
            1:size(y_task_period_all_tr,3) ) = y_task_period_all_tr;
        
        ms_rate_avg (block_no,1:960) = nanmean(ms_rate,1);
     
        ms_rate_avg_fix (block_no) = mean(ms_rate_avg(1:1000/5)); 
        ms_rate_avg_adaptor(block_no)  = mean(ms_rate_avg(1000/5+1:1000/5+3000/5)); 
        ms_rate_avg_stim(block_no)  = mean(ms_rate_avg(4000/5+1:4000/5+800/5)); 
 
        
    end
    %%
    x_task_period_all_tr_all_blocks_all_subj(sii, 1:8,1:61,1:4,1:9500) = x_task_period_all_tr_all_blocks;
    y_task_period_all_tr_all_blocks_all_subj(sii, 1:8,1:61,1:4,1:9500) = y_task_period_all_tr_all_blocks;
    
   
    
    %% organize by condition
    
    for ci = 1: 4
        indiz = find(block_cond == ci);
        
        x_per_cond(ci,1:900000) = [xpos_all_blocks(indiz(1),1:end) xpos_all_blocks(indiz(2),1:end)];
        y_per_cond(ci,1:900000) = [ypos_all_blocks(indiz(1),1:end) ypos_all_blocks(indiz(2),1:end)];
        ms_rate_per_cond(ci,1:122,1:960)  = [squeeze(ms_rate_all(indiz(1),:,:)); squeeze(ms_rate_all(indiz(2),:,:))];
        x_task_period_all_tr_all_cond(ci,1:122,1:4,1:9500) = [squeeze(x_task_period_all_tr_all_blocks(indiz(1),:,:,:)); squeeze(x_task_period_all_tr_all_blocks(indiz(2),:,:,:))];
        y_task_period_all_tr_all_cond(ci,1:122,1:4,1:9500) = [squeeze(y_task_period_all_tr_all_blocks(indiz(1),:,:,:)); squeeze(y_task_period_all_tr_all_blocks(indiz(2),:,:,:))];
        
        xx_all(1:900000,1) = x_per_cond(ci,:);
        xx_all(1:900000,2) = y_per_cond(ci,:);
        
        
        x_std_cond(ci) = nanstd(squeeze(xx_all(:,1)));
        y_std_cond(ci) = nanstd(squeeze(xx_all(:,2)));
        
        drift_vel_mean_cond(sii,ci) = (drift_vel_mean(sii,indiz(1))+  drift_vel_mean(sii,indiz(2)))/2;
        drift_vel_std_cond(sii,ci) = (drift_vel_std(sii,indiz(1))+  drift_vel_std(sii,indiz(2)))/2;
    end
    
    
    x_per_cond_ALL(sii,1:4,1:900000) = x_per_cond;
    y_per_cond_ALL(sii,1:4,1:900000) = y_per_cond;
    x_std_cond_ALL(sii,1:4) = x_std_cond;
    y_std_cond_ALL(sii,1:4) = y_std_cond;
    for ci = 1: 4
        x_task_period_all_tr_all_cond_all_subj(sii, ci,1:122,1:4,1:9500) =  x_task_period_all_tr_all_cond(ci,1:122,1:4,1:9500);
        y_task_period_all_tr_all_cond_all_subj(sii, ci,1:122,1:4,1:9500) =  y_task_period_all_tr_all_cond(ci,1:122,1:4,1:9500);
    end
   
   
    
    %% microsaccade rate analysis
   
    for ci = 1: 4
        
        ms_rate_fixE = squeeze(ms_rate_per_cond(ci,:,1:200));
        ms_rate_fix_cond(ci) = nanmedian(ms_rate_fixE(:));
        ms_rate_adaptorE = squeeze(ms_rate_per_cond(ci,:,201:800));
        ms_rate_adaptor_cond(ci) = nanmedian(ms_rate_adaptorE(:));
        ms_rate_stimE = squeeze(ms_rate_per_cond(ci,:,801:960));
        ms_rate_stim_cond(ci) = nanmedian(ms_rate_stimE(:));
    end
    
    
    ms_rate_fix_cond_subj_ALL(sii,1:4) = ms_rate_fix_cond;
    ms_rate_adaptor_cond_subj_ALL(sii,1:4) = ms_rate_adaptor_cond;
    ms_rate_stim_cond_subj_ALL(sii,1:4) = ms_rate_stim_cond;
    
    ms_rate_per_cond_ALL(sii,1:4,1:122,1:960) =  ms_rate_per_cond;
    %% plot per participant
    
    
end % end of participant loop

%% figures ms rate across all participants
figure
set(gcf, 'Position',[100 100 700 200])
colorz = [41 138 8; 152 191 100; 0 0 255; 137 195 255 ]/255;

marginsa2 = [0.09 0.08 0.08 0.1];
linewi = 1.1;
fontsz = 12;
set(gca, 'FontSize', fontsz)
for ci = 1: 4
    tight_subplot(1,3,1,1, guttera, marginsa2)
    bar(ci,nanmean(ms_rate_fix_cond_subj_ALL(:,ci)), 'FaceColor', colorz(ci,:), 'EdgeColor', colorz(ci,:)); hold on;
    errorbar(ci,nanmean(ms_rate_fix_cond_subj_ALL(:,ci)), nanstd(ms_rate_fix_cond_subj_ALL(:,ci))/sqrt(22),'Color', 'k','Capsize',0, 'Linewidth', linewi); hold on;
    if ci == 4
        title('Fixation')
    end
    box off
    set(gca, 'tickdir', 'out')
    set(gca, 'xtick', 1:1:4)
    set(gca, 'xticklabels', {})
    ylabel('ms rate (/sec)')
    set(gca, 'FontSize', fontsz)
    
    tight_subplot(1,3,1,2,guttera, marginsa2)
    bar(ci,nanmean(ms_rate_adaptor_cond_subj_ALL(:,ci)), 'FaceColor', colorz(ci,:), 'EdgeColor', colorz(ci,:)); hold on;
    errorbar(ci,nanmean(ms_rate_adaptor_cond_subj_ALL(:,ci)), nanstd(ms_rate_adaptor_cond_subj_ALL(:,ci))/sqrt(22),'Color', 'k','Capsize',0, 'Linewidth', linewi); hold on;
    if ci == 4
        title('Adaptor')
    end
    box off
    set(gca, 'tickdir', 'out')
    set(gca, 'xtick', 1:1:4)
    set(gca, 'xticklabels', {})
    set(gca, 'FontSize', fontsz)
    
    tight_subplot(1,3,1,3,guttera, marginsa2)
    bar(ci,nanmean(ms_rate_stim_cond_subj_ALL(:,ci)), 'FaceColor', colorz(ci,:), 'EdgeColor', colorz(ci,:)); hold on;
    errorbar(ci,nanmean(ms_rate_stim_cond_subj_ALL(:,ci)), nanstd(ms_rate_stim_cond_subj_ALL(:,ci))/sqrt(22),'Color', 'k','Capsize',0, 'Linewidth', linewi); hold on;
    if ci == 4
        title('Stimulus')
    end
    box off
    set(gca, 'tickdir', 'out')
    set(gca, 'xtick', 1:1:4)
    set(gca, 'xticklabels', {})
    set(gca, 'FontSize', fontsz)
end


%%
psname = 'microsaccade_rate_per_task_period_all_participants.pdf'
%print_pdf(psname)


%%
%x_per_cond_ALL(si,ci, 1:900000)
for si = 1: 22
    for ci = 1: 4
        x_per_sbj_per_cond_std(si,ci) = nanstd(squeeze(x_per_cond_ALL(si,ci,:)));
        y_per_sbj_per_cond_std(si,ci) = nanstd(squeeze(y_per_cond_ALL(si,ci,:)));
        
        x_per_sbj_per_cond_mean(si,ci) = nanmean(squeeze(x_per_cond_ALL(si,ci,:)));
        y_per_sbj_per_cond_mean(si,ci) = nanmean(squeeze(y_per_cond_ALL(si,ci,:)));
    end
end

%% planned post-hoc two-sided t-tests

[h,p,ci,stats]=ttest(nanmean(squeeze(x_per_cond_ALL(:,3,:)),2),nanmean(squeeze(x_per_cond_ALL(:,4,:)),2) )
% p = 0.216
[h,p,ci,stats]=ttest(nanmean(squeeze(y_per_cond_ALL(:,3,:)),2),nanmean(squeeze(y_per_cond_ALL(:,4,:)),2) )
% p = 0.790

[h,p,ci,stats]=ttest(nanmean(squeeze(x_per_cond_ALL(:,1,:)),2),nanmean(squeeze(x_per_cond_ALL(:,2,:)),2) )
% p = 0.982
[h,p,ci,stats]=ttest(nanmean(squeeze(y_per_cond_ALL(:,1,:)),2),nanmean(squeeze(y_per_cond_ALL(:,2,:)),2) )
% p = 0.332

[h,p,ci,stats]=ttest(nanmean(squeeze(x_per_cond_ALL(:,1,:)),2),nanmean(squeeze(x_per_cond_ALL(:,3,:)),2) )
% p = 0.661
[h,p,ci,stats]=ttest(nanmean(squeeze(y_per_cond_ALL(:,1,:)),2),nanmean(squeeze(y_per_cond_ALL(:,3,:)),2) )
% p = 0.299

[h,p,ci,stats]=ttest(nanmean(squeeze(x_per_cond_ALL(:,2,:)),2),nanmean(squeeze(x_per_cond_ALL(:,4,:)),2) )
% p = 0.484
[h,p,ci,stats]=ttest(nanmean(squeeze(y_per_cond_ALL(:,2,:)),2),nanmean(squeeze(y_per_cond_ALL(:,4,:)),2) )
% p = 0.779

%%
%close all;
figure(1)
set(gcf, 'Position', [100 100 700 700])

marginsa = [0.09 0.08 0.06 0.08]; % left right bottom top
guttera= [0.02 0.07];
Cond_list = { 'No-Adapt-See', 'No-Adapt-Believe','Adapt-See', 'Adapt-Believe'};
linewi_set = [3 2.5 2 1.5];

pts =linspace(-3, 3, 101);
for ci = 1:4
    figure(1)
    tight_subplot(4,4,1,ci, guttera, marginsa)
    set(gca, 'FontSize', fontsz)
    x_per_cond_all_sbj_ALL = squeeze(x_per_cond_ALL(:,ci,:));
    y_per_cond_all_sbj_ALL = squeeze(y_per_cond_ALL(:,ci,:));
    
    N = histcounts2(x_per_cond_all_sbj_ALL(:),y_per_cond_all_sbj_ALL(:), pts, pts);
    climz = [0 40000];
    imagesc(pts, pts,  N', climz); hold on;

    axis equal
    set(gca, 'xlim', pts([1 end]), 'ylim',pts([1 end]), 'YDir', 'normal')
    xlim([-3 3]); ylim([-3 3]); 
    plot(linspace(pts(1), pts(end),10), zeros(1,10),'--r'); hold on;
    plot( zeros(1,10), linspace(pts(1), pts(end),10),'--r'); hold on;
    color_circle = [1 0 0];
    h = circle(0,0,3, color_circle); hold on;
    title(Cond_list{ci})
    xlabel('x position (deg)', 'FontName', 'Helvetica', 'FontSize', fontsz)
    ylabel('y position (deg)', 'FontName', 'Helvetica', 'FontSize', fontsz)
    set(gca, 'tickdir', 'out')
    
    
    figure(1)
    tight_subplot(4,4,2,2,guttera, marginsa)
    set(gca, 'FontSize', fontsz)
    plot(nanmean(x_per_cond_all_sbj_ALL(:)), nanmean(y_per_cond_all_sbj_ALL(:)), 'ok', 'MarkerFaceColor', 'k'); hold on;
    x_std_cond(ci) = nanstd(x_per_cond_all_sbj_ALL(:));
    y_std_cond(ci) = nanstd(y_per_cond_all_sbj_ALL(:));
    x_neg = nanstd(x_per_cond_all_sbj_ALL(:));
    x_pos = nanstd(x_per_cond_all_sbj_ALL(:));
    y_neg = nanstd(y_per_cond_all_sbj_ALL(:));
    y_pos = nanstd(y_per_cond_all_sbj_ALL(:));
    plot(linspace(pts(1), pts(end),10), zeros(1,10),'--r'); hold on;
    plot( zeros(1,10), linspace(pts(1), pts(end),10),'--r'); hold on;
    color_circle = [1 0 0];
    h = circle(0,0,3, color_circle); hold on;
    errorbar(nanmean(x_per_cond_all_sbj_ALL(:)), nanmean(y_per_cond_all_sbj_ALL(:)), y_neg, y_pos, x_neg, x_pos, 'Color',colorz(ci,:), 'Linewidth',linewi_set(ci), 'Capsize',0); hold on;
    axis equal
    xlim([-3 3]); ylim([-3 3]);
    box off
    xlabel('x position (deg)', 'FontName', 'Helvetica', 'FontSize', fontsz)
    ylabel('y position (deg)', 'FontName', 'Helvetica', 'FontSize', fontsz)
    set(gca, 'tickdir', 'out')
    
end

cond_orders_sel = cond_orders([1:10 12:23], :);
sac_all_blocks_all_subj_22 = sac_all_blocks_all_subj([1:10 12:23],:);
sac_small_all_blocks_all_subj_22 = sac_small_all_blocks_all_subj; 
sac_large_all_blocks_all_subj_22 = sac_large_all_blocks_all_subj; 
for si = 1: 22
    cii = cond_orders_sel(si,:);
    block_cond = [cii(1) cii(1) cii(2) cii(2) cii(3) cii(3) cii(4) cii(4)];
    
    for ci = 1:4
        indiz = find(block_cond == ci);
        sac_all_cond_all_sbj{si,ci} = [cell2mat(sac_all_blocks_all_subj_22(si, indiz(1))); cell2mat(sac_all_blocks_all_subj_22(si, indiz(2)))];
        sac_rate_all_cond_all_subj(si,ci) = 1/2*(sac_rate_all_blocks_all_subj(si,indiz(1))+sac_rate_all_blocks_all_subj(si,indiz(2)));
        
        sac_small_all_cond_all_sbj{si,ci} = [cell2mat(sac_small_all_blocks_all_subj_22(si, indiz(1))); cell2mat(sac_small_all_blocks_all_subj_22(si, indiz(2)))];
        sac_small_rate_all_cond_all_subj(si,ci) = 1/2*(sac_small_rate_all_blocks_all_subj(si,indiz(1))+sac_small_rate_all_blocks_all_subj(si,indiz(2)));
        
        sac_large_all_cond_all_sbj{si,ci} = [cell2mat(sac_large_all_blocks_all_subj_22(si, indiz(1))); cell2mat(sac_large_all_blocks_all_subj_22(si, indiz(2)))];
        sac_large_rate_all_cond_all_subj(si,ci) = 1/2*(sac_large_rate_all_blocks_all_subj(si,indiz(1))+sac_large_rate_all_blocks_all_subj(si,indiz(2)));
        
        
        sac_rate_all_cond_all_subj(si,ci) = 1/2*(sac_rate_all_blocks_all_subj(si,indiz(1))+sac_rate_all_blocks_all_subj(si,indiz(2)));
        sac_small_all_cond_all_sbj{si,ci} = [cell2mat(sac_small_all_blocks_all_subj_22(si, indiz(1))); cell2mat(sac_small_all_blocks_all_subj_22(si, indiz(2)))];
        sac_large_all_cond_all_sbj{si,ci} = [cell2mat(sac_large_all_blocks_all_subj_22(si, indiz(1))); cell2mat(sac_large_all_blocks_all_subj_22(si, indiz(2)))];
    end
    
end
%%

figure(1)
for ci = 1: 4
    tight_subplot(4,4,3,ci,guttera, marginsa)
    set(gca, 'FontSize', fontsz)
    for si = 1:22
        sacE = sac_small_all_cond_all_sbj{si,ci};
        if size(sac_small_all_cond_all_sbj{si,ci},1)>0
            plot(sqrt(sacE(:,6).^2+ sacE(:,7).^2),sacE(:,3), 'o','MarkerSize',0.5, 'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:) ); hold on; % 'EdgeColor', colorz(ci,:)
        end
    end
    box off
    set(gca, 'tickdir', 'out')
    xlabel('amplitude (deg)', 'FontName', 'Helvetica', 'FontSize', fontsz)

    if ci == 1
        ylabel('peak velocity (deg/sec)', 'FontName', 'Helvetica', 'FontSize', fontsz)
    else
        set(gca, 'yticklabels', [])
    end

    xlim([0 1.1]); ylim([0 110]);
 
end
%%
figure(1)
pars_scatter = 0.15;
mszi = 3;

tight_subplot(4,4,4,2,guttera, marginsa)
set(gca, 'FontSize', fontsz)
for ci = 1: 4
   
    if ci <=2
        plot((ci)*ones(1, 22)- pars_scatter+ 2*pars_scatter*rand(1,22),  sac_small_rate_all_cond_all_subj(:,ci), 'o','MarkerSize',mszi, 'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:)); hold on;
        bar(ci,mean(sac_small_rate_all_cond_all_subj(:,ci)), 'FaceColor','none', 'EdgeColor', colorz(ci,:), 'LineWidth',2); hold on;
        ebi(ci) = errorbar(ci,mean(sac_small_rate_all_cond_all_subj(:,ci)),std(sac_small_rate_all_cond_all_subj(:,ci))/sqrt(22), 'Color', 'k', 'Linewidth',1.5, 'Capsize', 0); hold on;
    elseif ci >=3
        plot((ci+1)*ones(1, 22)- pars_scatter+ 2*pars_scatter*rand(1,22),  sac_small_rate_all_cond_all_subj(:,ci), 'o','MarkerSize',mszi, 'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:)); hold on;
        bar(ci+1,mean(sac_small_rate_all_cond_all_subj(:,ci)), 'FaceColor','none', 'EdgeColor', colorz(ci,:), 'LineWidth',2); hold on;
        ebi(ci) = errorbar(ci+1,mean(sac_small_rate_all_cond_all_subj(:,ci)),std(sac_small_rate_all_cond_all_subj(:,ci))/sqrt(22), 'Color', 'k', 'Linewidth',1.5, 'Capsize', 0); hold on;
    end
end
set(gca, 'tickdir', 'out')
ylabel('microsaccade rate (/s)', 'FontName', 'Helvetica', 'FontSize', fontsz)
set(gca, 'xticklabels', [])
box off
xlim([0.4 5.6])

%%
psname = 'Eye_plot_with_everything.pdf'
%print_pdf(psname)
%% microsaccades planned post-hoc two-sided t-tests
[h,p,ci,stats]=ttest(sac_small_rate_all_cond_all_subj(:,3),sac_small_rate_all_cond_all_subj(:,4) )
% p = 0.355
[h,p,ci,stats]=ttest(sac_small_rate_all_cond_all_subj(:,1),sac_small_rate_all_cond_all_subj(:,2) )
% p = 0.648
[h,p,ci,stats]=ttest(sac_small_rate_all_cond_all_subj(:,1),sac_small_rate_all_cond_all_subj(:,3) )
% p = 0.545
[h,p,ci,stats]=ttest(sac_small_rate_all_cond_all_subj(:,2),sac_small_rate_all_cond_all_subj(:,4) )
% p = 0.441
%% bar plot for several measures, as long as they have transformed into 22 by 4 in the correct order
vals_bar_plot = sac_small_rate_all_cond_all_subj;
%vals_bar_plot = drift_vel_mean_cond;
%vals_bar_plot = drift_vel_std_cond;
close all;
figure(101)
set(gcf, 'Position', [100 100 460 300])
pars_scatter = 0.15;
mszi = 5;
guttera_bar = [0.09 0.09];
marginsa_bar = [0.100    0.100    0.100    0.130]; %[LEFT RIGHT BOTTOM TOP]

tight_subplot(1,1,1,1,guttera_bar, marginsa_bar)
set(gca, 'FontSize', fontsz)
for ci = 1: 4
    
    if ci <=2
        b_leg(ci) =plot((ci)*ones(1, 22)- pars_scatter+ 2*pars_scatter*rand(1,22),  vals_bar_plot(:,ci), 'o','MarkerSize',mszi, 'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:)); hold on;
        bar(ci,mean(vals_bar_plot(:,ci)), 'FaceColor','none', 'EdgeColor', colorz(ci,:), 'LineWidth',2); hold on;
        ebi(ci) = errorbar(ci,mean(vals_bar_plot(:,ci)),std(vals_bar_plot(:,ci))/sqrt(22), 'Color', 'k', 'Linewidth',1.5, 'Capsize', 0); hold on;
    elseif ci >=3
        b_leg(ci) = plot((ci+1)*ones(1, 22)- pars_scatter+ 2*pars_scatter*rand(1,22),  vals_bar_plot(:,ci), 'o','MarkerSize',mszi, 'MarkerFaceColor', colorz(ci,:), 'MarkerEdgeColor', colorz(ci,:)); hold on;
        bar(ci+1,mean(vals_bar_plot(:,ci)), 'FaceColor','none', 'EdgeColor', colorz(ci,:), 'LineWidth',2); hold on;
        ebi(ci) = errorbar(ci+1,mean(vals_bar_plot(:,ci)),std(vals_bar_plot(:,ci))/sqrt(22), 'Color', 'k', 'Linewidth',1.5, 'Capsize', 0); hold on;
    end
end
set(gca, 'tickdir', 'out')
ylabel('microsaccade rate (/s)', 'FontName', 'Helvetica', 'FontSize', fontsz)
%ylabel('drift vel  ', 'FontName', 'Helvetica', 'FontSize', fontsz)
%ylabel('drift vel std ', 'FontName', 'Helvetica', 'FontSize', fontsz)
box off
lege = legend([b_leg(1),b_leg(2), b_leg(3), b_leg(4)], {'No-Adapt-See','No-Adapt-Believe','Adapt-See', 'Adapt-Believe' }, 'FontSize',12)
legend boxoff
set(lege, 'Position', [0.7013    0.7950    0.2457    0.1650])
set(gca, 'xtick', [])
set(gca, 'xticklabels', [])
xlim([0.4 5.6])
%%
psname = 'sacc_small_rate_robustness.pdf'
%psname = 'drift_vel_mean_robustness.pdf'
%psname = 'drift_vel_std_robustness.pdf'
%psname = 'microsacc_vel_robustness.pdf'
%print_pdf(psname)
%%

load('psych_curves_fitting_m2_201_E2.mat')
MAE_strength =  -(mu_est_all(indi_sel,3) - mu_est_all(indi_sel,1));
%%
[r,p] = corr(sac_small_rate_all_cond_all_subj(:,3), MAE_strength, 'type', 'Spearman')
[r,p] = corr(drift_vel_std_cond(:,3), MAE_strength, 'type', 'Spearman')
%%
psname = 'fixation_stability_across_all_task_periods_all_participantsEFF.pdf'
%print_pdf(psname)




%%

%x_task_period_all_tr_all_cond_all_subj(sii, ci,1:122,1:4,1:9500)
for si = 1: 22
    for ci = 1: 4
        for pi = 1:4
            
            xsel = squeeze(x_task_period_all_tr_all_cond_all_subj(si, ci,1:122,pi,1:9500));
            ysel = squeeze(y_task_period_all_tr_all_cond_all_subj(si, ci,1:122,pi,1:9500));
            
            x_per_sbj_per_cond_per_period_std(si,ci,pi) = nanstd(xsel(:));
            y_per_sbj_per_cond_per_period_std(si,ci,pi) = nanstd(ysel(:));
            
            x_per_sbj_per_cond_per_period_mean(si,ci,pi) = nanmean(xsel(:));
            y_per_sbj_per_cond_per_period_mean(si,ci,pi) = nanmean(ysel(:));
        end
    end
end
%%
for pi = 1:4
    
    [h,p_per_period_x_mean(pi),t, stats_mean_x(pi)]=ttest(x_per_sbj_per_cond_per_period_mean(:,1,pi), x_per_sbj_per_cond_per_period_mean(:,2,pi))
    [h,p_per_period_x_std(pi),t, stats_std_x(pi )]=ttest(x_per_sbj_per_cond_per_period_std(:,3,pi), x_per_sbj_per_cond_per_period_std(:,4,pi))
    
    [h,p_per_period_y_mean(pi),t, stats_mean_y(pi)]=ttest(y_per_sbj_per_cond_per_period_mean(:,2,pi), y_per_sbj_per_cond_per_period_mean(:,4,pi))
    [h,p_per_period_y_std(pi),t, stats_std_y(pi)]=ttest(y_per_sbj_per_cond_per_period_std(:,3,pi), y_per_sbj_per_cond_per_period_std(:,4,pi))
   
end


%%

%% two-way repeated-measures anova
%22 by 4 variables of interest
%x_per_sbj_per_cond_mean; y_per_sbj_per_cond_mean
%sac_small_rate_all_cond_all_subj
%x_per_sbj_per_cond_std, y_per_sbj_per_cond_std
%drift_vel_std_cond
clear vals; clear group; clear y;
Nsubj = 22;
vals = [];
vals0 = x_per_sbj_per_cond_mean;
vals = [vals0(:,1)' vals0(:,2)' vals0(:,3)' vals0(:,4)']';


Subjects = repmat([1:Nsubj],1,4)';
group = [[ones(1,2*Nsubj) 2*ones(1,2*Nsubj)]'...
    [repmat([ones(1,Nsubj) 2*ones(1,Nsubj)],1,2)]'];

y = vals;

g1 = group(:,1);
g2 = group(:,2);
%%

[p1,table1,stats1] = anovan(y,{g1,g2, Subjects}, 'random',3,'varnames', {'No_Adapt_Vs_Adapt', 'See_Vs_Bel', 'Subject'});
[p2_full,table2_full,stats2_full] = anovan(y,{g1,g2, Subjects}, 'random',3,'varnames', {'No_Adapt_Vs_Adapt',  'See_Vs_Bel', 'Subject'}, 'model', 'full');
[p3,table3,stats3] = anovan(y,{g1,g2}, 'random',2,'varnames', {'No_Adapt_Vs_Adapt', 'See_Vs_Bel'}, 'model', 'full');
%%
%indi_sel = 1:1:22;
indi_sel = [1:1:3 5:1:22];  % exclude P4 bc very low microsacc rate, close to 0 and unusually high drift vel
MAE_strength = -(mu_est_all(indi_sel,3) - mu_est_all(indi_sel,1));
%MAE_strength = -(mu_est_all(indi_sel,3));



figure;
set(gcf,'Position',[100 100 600 300])
guttera3 = [0.09 0.08];
marginsa3 = [0.120    0.100    0.140    0.1000]; %[LEFT RIGHT BOTTOM TOP]
tight_subplot(1,2,1,1,guttera3, marginsa3)
plot(sac_small_rate_all_cond_all_subj(indi_sel,3),MAE_strength, 'o', 'MarkerFaceColor','k','MarkerEdgeColor', 'k'); hold on;
box off
brob = robustfit(sac_small_rate_all_cond_all_subj(indi_sel,3),MAE_strength);
plot(sac_small_rate_all_cond_all_subj(indi_sel,3), brob(1)+ brob(2)*sac_small_rate_all_cond_all_subj(indi_sel,3), '-k', 'Linewidth',1); hold on;
set(gca, 'FontSize',12)
set(gca, 'tickdir', 'out')
xlabel('microsaccade rate (/s)')
ylabel('MAE strength')
[rho1,p1] = corr(sac_small_rate_all_cond_all_subj(indi_sel,3),MAE_strength,'type', 'Spearman');
title(['\rho = ', num2str(rho1, '%.2f'), ', p = ', num2str(p1, '%.2f')] )



tight_subplot(1,2,1,2,guttera3, marginsa3 )
plot(drift_vel_std_cond(indi_sel,3),MAE_strength, 'o', 'MarkerFaceColor','k','MarkerEdgeColor', 'k'); hold on;
box off
brob = robustfit(drift_vel_std_cond(indi_sel,3),MAE_strength);
plot(drift_vel_std_cond(indi_sel,3), brob(1)+ brob(2)*drift_vel_std_cond(indi_sel,3), '-k', 'Linewidth',1); hold on;
set(gca, 'FontSize',12)
set(gca, 'tickdir', 'out')
%xlim([0 20])
xlabel('standard deviation drift velocity (deg/s)')
ylabel('MAE strength')
[rho2,p2] = corr(drift_vel_std_cond(indi_sel,3),MAE_strength,'type', 'Spearman');
title(['\rho = ', num2str(rho2, '%.2f'), ', p = ', num2str(p2, '%.2f')] )

%%
psname = 'corrs_robustness_exclude_P4.pdf'; 
%print_pdf(psname)



