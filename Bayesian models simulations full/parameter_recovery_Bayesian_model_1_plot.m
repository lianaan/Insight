clear all; close all;

mu_A = -0.055;

for indi_1 = 1:6
    
    for indi_2 = 1:6
        
        for indi_3 = 1:6
            
            for repi = 1:3
                
                load(['sim_data_M1_',num2str(indi_1),'_', num2str(indi_2),'_', num2str(indi_3), '_', num2str(repi) '.mat']);
                params_in(indi_1, indi_2, indi_3, repi,1:4) = params_input;
                nll_in(indi_1, indi_2, indi_3, repi) = -LL_sum;
                
                
                load(['FITS_sim_data_M1_',num2str(indi_1),'_', num2str(indi_2),'_', num2str(indi_3), '_', num2str(repi) '.mat']);
                params_out_N(indi_1, indi_2, indi_3, repi,1:4) = [mu_A params_fit];
                nll_out_N(indi_1, indi_2, indi_3, repi) =  nll;
                params_out(indi_1, indi_2, indi_3, repi,1:4) = [mu_A params_fit];
                nll_out(indi_1, indi_2, indi_3, repi) =  nll;
                
            end
            
        end
        
    end
    
end



%%

indi_diff = 1:1:216;

params2_set = log([ 0.01 0.02 0.04 0.005 0.1 0.12]);
params3_set = [0.85 0.9 0.95 0.8 0.7 0.55];
params4_set = [-0.15 -0.1 -0.05 0 0.05 -0.2];

figure
marginsa = [0.11 0.03 0.09 0.07]; %MARGINS = [LEFT RIGHT BOTTOM TOP]
guttera = [0.09 0.09];
set(gcf, 'Position', [100 100 600 240])

mszs = 6*[ 4.5 4 3.5   6  3 2 ];
colorz = [252 78 42; 253 141 60; 254 178 76;  177 0 38; 254 217 118; 255 237 160]/255;
transpa = 0.5;

fontsz = 12;

tight_subplot(1,3,1,1, guttera, marginsa)
for indi_1 = 1:6
    
    muL_in = reshape(params_in(indi_1, :, :,:,4),1,108);
    muL_out = reshape(params_out(indi_1, :, :,:,4),1,108);
    scatter(muL_in, muL_out,mszs(indi_1), 'o',   'MarkerFaceColor', colorz(indi_1,:),'MarkerEdgeColor',colorz(indi_1,:),'MarkerFaceAlpha', transpa, 'MarkerEdgeAlpha', transpa); hold on;
    
    if indi_1 ==6
        axis equal
        xlim([-0.22 0.07])
        ylim([-0.22 0.07])
        plot(sort(params4_set),sort(params4_set),'--k'); hold on;
        box off
        set(gca, 'tickdir', 'out')
        xlabel('\mu_{likelihood} ground truth', 'FontName', 'Helvetica', 'FontSize', fontsz)
        ylabel('\mu_{likelihood} recovered', 'FontName', 'Helvetica', 'FontSize', fontsz)
    end
end

tight_subplot(1,3,1,2, guttera, marginsa)
for indi_1 = 1:6
    
    noise_in = reshape(exp(params_in(indi_1, :, :,:,2)),1,108);
    noise_out = reshape(exp(params_out(indi_1, :, :,:,2)),1,108);
    scatter(noise_in,noise_out,mszs(indi_1), 'o',  'MarkerFaceColor', colorz(indi_1,:),'MarkerEdgeColor',colorz(indi_1,:),'MarkerFaceAlpha', transpa, 'MarkerEdgeAlpha', transpa); hold on;
    if indi_1 ==6
        axis equal
        xlim([0.001 0.15])
        ylim([0.001 0.15])
        plot(exp(sort(params2_set)),exp(sort(params2_set)), '--k'); hold on;
        box off
        set(gca, 'tickdir', 'out')
        xlabel(' \sigma_{encoding} ground truth', 'FontName', 'Helvetica', 'FontSize', fontsz)
        ylabel(' \sigma_{encoding} recovered', 'FontName', 'Helvetica', 'FontSize', fontsz)
    end
end



tight_subplot(1,3,1,3, guttera, marginsa)
for indi_1 = 1:6
    
    kconf_in = reshape(params_in(indi_1, :, :,:,3),1,108);
    kconf_out = reshape(params_out(indi_1, :, :,:,3),1,108);
    scatter(kconf_in, kconf_out,mszs(indi_1), 'o', 'MarkerFaceColor', colorz(indi_1,:),'MarkerEdgeColor',colorz(indi_1,:),'MarkerFaceAlpha', transpa, 'MarkerEdgeAlpha', transpa); hold on;
    if indi_1 ==6
        axis equal
        xlim([0.500 0.9999])
        ylim([0.500 0.9999])
        box off
        plot(sort(params3_set),sort(params3_set), '--k'); hold on;
        set(gca, 'tickdir', 'out')
        xlabel('k_{confidence} ground truth', 'FontName', 'Helvetica', 'FontSize', fontsz)
        ylabel('k_{confidence} recovered', 'FontName', 'Helvetica', 'FontSize', fontsz)
    end
end

%%
psname = 'Parameter_recovery_M1_enhanced_fin_3.pdf'
print_pdf(psname)

%%

% corr for noise

params_in_noise = [];
params_out_noise = [];

for indi_1 = 1:6
    
    params_in_noise = [params_in_noise repmat(exp(params2_set(indi_1)),1,6*6*3)];
    
    for indi_2 = 1:6
        
        for indi_3 = 1:6
            
            for repi = 1:3
                
                params_out_noise = [params_out_noise exp(params_out(indi_1, indi_2, indi_3, repi,2))];
                
            end
            
        end
        
    end
    
end
%%
[r_noise,p_noise]= corr(params_in_noise', params_out_noise', 'type', 'Spearman')
%r_noise = 0.97

%%
params_in_muL = [];
params_out_muL = [];

for indi_3 =1:6
    
    params_in_muL = [params_in_muL repmat((params4_set(indi_3)),1,6*6*3)];
    
    for indi_1 = 1:6
        
        for indi_2 = 1:6
            
            for repi = 1:3
                
                params_out_muL = [params_out_muL (params_out(indi_1, indi_2, indi_3, repi,4))];
                
            end
            
        end
        
    end
    
end
%%
[r_muL,p_muL]= corr(params_in_muL', params_out_muL', 'type', 'Spearman')
%r_muL = 0.98
%%
params_in_kconf = [];
params_out_kconf = [];

for indi_2 = 1:6
    params_in_kconf = [params_in_kconf repmat((params3_set(indi_2)),1,6*6*3)];
    
    for indi_1 = 1:6
        
        for indi_3 = 1:6
            
            for repi = 1:3
                
                params_out_kconf = [params_out_kconf (params_out(indi_1, indi_2, indi_3, repi,3))];
                
            end
            
        end
        
    end
    
end
%%
[r_kconf,p_kconf]= corr(params_in_kconf', params_out_kconf', 'type', 'Spearman')
%r_kconf = 0.85
