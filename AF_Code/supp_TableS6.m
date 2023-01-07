
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script will produce Table S6 of the Supplementary Information
%
% It will also produce Monte Carlo approximations of the distribution of
% the various estimators under the null (no trend) and compare with the
% test statistic obtained from the original data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% (c) Mikkel Bennedsen (2023)
%
% This code can be used, distributed, and changed freely. Please cite Bennedsen,
% Hillebrand, and Koopman (2022): "Is there evidence of a trend in the CO2 airborne fraction?".
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NB1: The data are taken from van Marle et al. (2022): 
%      "New land-use-change emissions indicate a declining CO2 airborne fraction", Nature 603, 450–454 (2022)
%
% NB2: Sen's estimator is implemented using the "ktaub" function 
%      written by Jeff Burkey (2013). This code can be freely distributed, but 
%      please see the license for this code in "Functions/ktaub/license.txt.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('Data');
addpath(genpath('Functions'));
%% Init
filenam = 'Data/Marle_et_al_Nature_AirborneFraction_Datasheet.xlsx';

tit_str = {'GCP-raw','GCP-filter','H&N-raw','H&N-filter','New-raw','New-filter'};

alpha = 0.05;
wantplot = 0;

B = 1e5; % number of bootstrap replications
%% Load data
dat = xlsread(filenam,6);

%% Make data
t = dat(:,1);
AF = dat(:,2);
AF_corr = dat(:,4);

AF2 = dat(:,2+4); % HN
AF2_corr = dat(:,4+4);

AF3 = dat(:,2+8); % GCP
AF3_corr = dat(:,4+8);

%% Do analysis
sen_b = nan(B,6);
bhat_b = nan(B,6);
b2hat = nan(6,1);
senhat = nan(6,1);
sen_p = nan(6,1);
bhat_p = nan(6,1);
for i = 1:6
    disp(i/6);
    if i == 1 % Data: GCP (raw)
        y = AF3;
    elseif i == 2 % Data: GCP (filter)
        y = AF3_corr;
    elseif i == 3 % Data: H&N (raw)
        y = AF2;
    elseif i == 4 % Data: H&N (filter)
        y = AF2_corr;
    elseif i == 5 % Data: New (raw)
        y = AF;
    elseif i == 6 % Data: New (filter)
        y = AF_corr;
    end 
    datain = [t,y];
    
    for b = 1:B
        %% Permute original data to obtain data without trend
        indx = randperm(length(y)); 
        y_b = y(indx);

        %% Sen's estimator
        datain_b = [t,y_b];
    
        [taub, tau, h, sig, Z, S, sigma, sen, n, senplot, CIlower, CIupper, D, Dall, C3, nsigma] = ktaub(datain_b, alpha, wantplot);
    
        sen_b(b,i) = sen;

        %% OLS estimator of slope b in y = a + b*t
        X = [ones(length(y_b),1),t-t(1)];
        bhat = (X'*X)\X'*y_b;
        
        bhat_b(b,i) = bhat(2);

    end
    
    %% Calculate estimates on original data
    [taub, tau, h, sig, Z, S, sigma, sen, n, senplot, CIlower, CIupper, D, Dall, C3, nsigma] = ktaub(datain, alpha, wantplot);
    
    X = [ones(length(y),1),t-t(1)];
    bhat = (X'*X)\X'*y;

    senhat(i) = sen;
    b2hat(i) = bhat(2);
    
    %% Test H0: b=0
    sen_p(i,1)  = mean(abs(sen_b(:,i))>abs(senhat(i))); % Two-sided: H1: b!=0
    bhat_p(i,1) = mean(abs(bhat_b(:,i))>abs(b2hat(i)));

    sen_p(i,2)  = mean((sen_b(:,i))>(senhat(i))); % one-sided: H1: b>0
    bhat_p(i,2) = mean((bhat_b(:,i))>(b2hat(i)));

    sen_p(i,3)  = mean((sen_b(:,i))<(senhat(i))); % one-sided: H1: b<0
    bhat_p(i,3) = mean((bhat_b(:,i))<(b2hat(i)));
    

    fig1 = figure(1);
    subplot(3,2,i);
    histogram(sen_b(:,i), 'Normalization', 'pdf'), hold on
    ax = gca;

    plot(senhat(i)*ones(100,1),linspace(0,ax.YLim(2),100),'r-','LineWidth',2), hold on

    
    text(ax.XLim(1)*0.95,ax.YLim(2)*0.95,['p-val (two-sided): ',num2str(sen_p(i,1))],'FontSize',8);
    if senhat(i)<0
        text(ax.XLim(1)*0.95,ax.YLim(2)*0.85,['p-val (one-sided): ',num2str(sen_p(i,3))],'FontSize',8);
    else
        text(ax.XLim(1)*0.95,ax.YLim(2)*0.85,['p-val (one-sided): ',num2str(sen_p(i,2))],'FontSize',8);
    end
    title(['Sens estimator: ',tit_str{i}],'FontSize',8);

    fig2 = figure(2);
    subplot(3,2,i);
    histogram(bhat_b(:,i), 'Normalization', 'pdf'), hold on
    ax = gca;

    plot(b2hat(i)*ones(100,1),linspace(0,ax.YLim(2),100),'r-','LineWidth',2), hold on

    

    text(ax.XLim(1)*0.95,ax.YLim(2)*0.95,['p-val (two-sided): ',num2str(bhat_p(i,1))],'FontSize',8);
    if b2hat(i)<0
        text(ax.XLim(1)*0.95,ax.YLim(2)*0.85,['p-val (one-sided): ',num2str(bhat_p(i,3))],'FontSize',8);
    else
        text(ax.XLim(1)*0.95,ax.YLim(2)*0.85,['p-val (one-sided): ',num2str(bhat_p(i,2))],'FontSize',8);
    end
    title(['OLS estimator: ',tit_str{i}],'FontSize',8);
end


%% Output table
disp(' ');
disp(' Bootstrapped p-values: Table S6 in Supplementary Material')
disp('        Mann-Kendall tests    Regression test')
disp('      MK        MK+       MK-     Slope1')
disp([sen_p,bhat_p(:,1)])


