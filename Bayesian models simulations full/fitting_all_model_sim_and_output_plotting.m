clear all; close all;
%%
%{
modi = 1; % fit Bayesian model 1
for i1 = 1:4
    for i2 = 1:3
        for i3 = 1:3
            for i4 = 1:5
                mpr_indices = [i1 i2 i3 i4];
                Fitting_modell_sim(mpr_indices,modi);
            end
        end
    end
end
%}

%%
curr_dir = pwd;
cd ([pwd, '/fits Afull/'])

mu_A = -0.055;
mu_like_set = [0 -0.05 -0.1 -0.15];
sigma_enc_set = [0.01 0.05 0.1];
sigma_A_set = [0.01 0.03 0.05];

Ntrials = 121;%1000;
prob_right = 0.5;
k_confidence = 0.9;
conf_in = k_confidence;


for i1 = 1:4
    for i2 = 1:3
        for i3 = 1:3
            for i4 = 1:5
                
                
                mu_like_in(i1,i2,i3,i4) = mu_like_set(i1);
                sigma_enc_in(i1,i2,i3,i4) = sigma_enc_set(i2);
                sigma_A_in(i1,i2,i3,i4) = sigma_A_set(i3);
                runi_in(i1,i2,i3,i4) = i4;
                
                load(['Afits_sim_full_model_1', '_mu_like_', num2str(mu_like_set(i1)), '_sigma_enc_', num2str(sigma_enc_set(i2)), '_sigma_A_', num2str(sigma_A_set(i3)), '_rep_no_', num2str(i4) '.mat'])
                mu_A_fit(i1,i2,i3,i4) = mu_A;
                sigma_all_fit(i1,i2,i3,i4) = params_fit_best(2);
                conf_fit(i1,i2,i3,i4) = params_fit_best(3);
                mu_like_fit(i1,i2,i3,i4) = params_fit_best(4);
            end
        end
    end
end

%%

marginsa = [0.08 0.1 0.05 0.05]; %MARGINS = [LEFT RIGHT BOTTOM TOP]
guttera = [0.1 0.07];
linewi = 1.1;

close all;
figure
set(gcf, 'Position', [100 100 760 380])
set(gca, 'FontSize',14)


colorz = [41 138 8; 152 191 100; 0 0 255; 137 195 255 ]/255;
colorz_sigma_enc = [255 35 0; 242 120 30;  254 196 79; 255 227 149]/255;
colorz_sigma_A = [44 127 184; 127 205 187; 237 248 177]/255;
msz_all = [6.3 4.4 3.7];
tight_subplot(1,3,1, 1, guttera, marginsa)
for i1 = 1: 4 % for the mu_likelihood, modi
    for i3 = 1:3 % for the sigma A
        pl1(i3) = plot(reshape(squeeze(mu_like_in(i1,:,i3,:)),15,1), reshape(squeeze(mu_like_fit(i1,:,i3,:)),15,1), 'o', 'MarkerSize', msz_all(i3), 'MarkerFaceColor', colorz_sigma_A(i3,:),'MarkerEdgeColor', colorz_sigma_A(i3,:),'Linewidth', linewi); hold on; %'o-', 'MarkerFaceColor', colorz(kk,:),'MarkerEdgeColor', colorz(kk,:),
    end
end
plot(mu_like_in_lin, mu_like_in_lin, '--k'); hold on;
box off
axis equal
set(gca, 'tickdir', 'out')
la1 = legend([pl1(1) pl1(2) pl1(3)], {' \sigma_{A} = 0.01', ' \sigma_{A} = 0.03', ' \sigma_{A} = 0.05'}); hold on;
legend boxoff
set(la1, 'Position', [0.1387    0.7559    0.0917    0.1658])
xlim([-0.16 0.02]); ylim([-0.16  0.02])
xlabel('\mu_{likelihood} in')
ylabel('\mu_{likelihood} out')
mulk = get(gca,'XTickLabel');
set(gca,'XTickLabel',mulk,'FontName','Helvetica','fontsize',12)

tight_subplot(1,3,1, 2, guttera, marginsa)
set(gca, 'FontSize',14)
for i2 = 1:3 % for the sigma encoding
    for i3 = 1:3 % for the sigma A
        loglog(reshape(squeeze((sigma_enc_in(:,i2,i3,:))),20,1), reshape(squeeze(exp(sigma_all_fit(:,i2,i3,:))),20,1), 'o', 'MarkerSize', msz_all(i3), 'MarkerFaceColor', colorz_sigma_A(i3,:),'MarkerEdgeColor', colorz_sigma_A(i3,:),'Linewidth', linewi); hold on; %'o-', 'MarkerFaceColor', colorz(kk,:),'MarkerEdgeColor', colorz(kk,:),
    end
end
loglog((sigma_enc_in_lin),(sigma_enc_in_lin) , '--k'); hold on;
plot((sigma_enc_in_lin),(sigma_enc_in_lin) , '--k'); hold on;
p1=plot(0.125, exp(-3.97), '^', 'Color', colorz(1,:)); hold on;
p2=plot(0.125, exp(-4.03), '^', 'Color', colorz(2,:)); hold on;
p3=plot(0.125, exp(-3.13), '^', 'Color', colorz(3,:)); hold on;
p4=plot(0.125, exp(-2.5530), '^', 'Color', colorz(4,:)); hold on;
box off
axis equal
set(gca, 'tickdir', 'out')
xlabel(' \sigma_{encoding} in')
ylabel(' \sigma_{encoding} out')
xlim([exp(-5) 0.12]); ylim([exp(-5) 0.12])
e = get(gca,'XTickLabel');
set(gca,'XTickLabel',e,'FontName','Helvetica','fontsize',12)

tight_subplot(1,3,1, 3, guttera, marginsa)
set(gca, 'FontSize',18)
for i2 = 1:3 % for the sigma encoding 
    for i3 = 1:3 % for the sigma A
        pl3(i2)= loglog(reshape(squeeze((sigma_A_in(:,i2,i3,:))),20,1), reshape(squeeze(exp(sigma_all_fit(:,i2,i3,:))),20,1), 'o', 'MarkerSize', msz_all(i2), 'MarkerFaceColor', colorz_sigma_enc(i2,:),'MarkerEdgeColor', colorz_sigma_enc(i2,:),'Linewidth', linewi); hold on; %'o-', 'MarkerFaceColor', colorz(kk,:),'MarkerEdgeColor', colorz(kk,:),
    end
end
loglog((sigma_A_in_lin),(sigma_A_in_lin) , '--k'); hold on;
plot(log(sigma_enc_in_lin), -3.97 *ones(1, length(sigma_enc_in_lin)), '--', 'Color', colorz(1,:)); hold on;
plot(log(sigma_enc_in_lin), -4.03 *ones(1, length(sigma_enc_in_lin)), '--', 'Color', colorz(2,:)); hold on;
plot(log(sigma_enc_in_lin), -3.13 *ones(1, length(sigma_enc_in_lin)), '--', 'Color', colorz(3,:)); hold on;
plot(log(sigma_enc_in_lin), -2.5530 *ones(1, length(sigma_enc_in_lin)), '--', 'Color', colorz(4,:)); hold on;
box off
axis equal
set(gca, 'tickdir', 'out')
la3 = legend([pl3(1) pl3(2) pl3(3)], {' \sigma_{encoding} = 0.01', ' \sigma_{encoding} = 0.05', ' \sigma_{encoding} = 0.1'}); hold on;
legend boxoff
set(la3, 'Position', [0.7487    0.7559    0.0917    0.1658])
xlabel(' \sigma_{A}  in')
ylabel(' \sigma_{encoding} out')
xlim([exp(-5) 0.11]); ylim([exp(-5) 0.11])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Helvetica','fontsize',12)
%%
psname = 'param_rec_full_model.pdf'
print_pdf(psname)
%%
[r,p] = corr(mu_like_in_lin', mu_like_fit_lin', 'type', 'Spearman')
% r = 0.88, p = 3.4 * 10^(-59)


%%
figure

colorbar
color1 = [208 150 150];
color2 = [243 199 167];
color3 = [252 235 195];

[a_enc,b_enc,c_enc]=unique(sigma_enc_in_lin);
[a_A,b_A,c_A]=unique(sigma_A_in_lin);
for ii1 = 1:3
    indi_enc = c_enc == ii1;
    for jj1 = 1:3
        indi_A =  c_A == jj1;
        sigma_fit_2d(ii1,jj1) = mean(sigma_fit_lin(indi_enc & indi_A));
    end
end
imagesc(log(a_enc),log(a_A),sigma_fit_2d'); colorbar
title('log sigma fit')
set(gca, 'YDir', 'Normal')
xlabel('log sigma enc')
ylabel('log sigma A')
psname = 'sigma_fit_with_sigma_enc_and_sigma_A.pdf'
%print_pdf(psname)
