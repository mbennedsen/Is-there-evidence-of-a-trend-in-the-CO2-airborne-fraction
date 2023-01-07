%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script will run the estimation behind, and plot, Figure 1 of the main paper.
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc; clear; close all;
addpath('Data');
addpath(genpath('Functions'));
%% Init
filenam = 'Data/Marle_et_al_Nature_AirborneFraction_Datasheet.xlsx';

tit_str = {'GCP-raw','GCP-filter','H&N-raw','H&N-filter','New-raw','New-filter'};

alpha = 0.05;
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
pVals = nan(6,7);
for i = 1:6
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
        
    
    %% OLS (no breaks): y = a + b*t
    X = [ones(length(y),1),t-t(1)];
    bhat = (X'*X)\X'*y;
    
    fig2 = figure(2);
    subplot(3,2,i);
    plot(t,y,'k-','LineWidth',1.5), hold on
    plot(t,X*bhat,'b-.','LineWidth',1.5), hold on
    

    %% OLS  w. break in intercept: y = a1 + a2*I(t) + b*t 
    if mod(i,2) == 0
        br_date = 1990;
    else
        br_date = 1988;
    end
    X = [ones(length(y),1),cumsum(t==br_date),t-t(1)];
    bhat = (X'*X)\X'*y;
    
    plot(t,X*bhat,'r-','LineWidth',1.5), hold on
    %xlabel('Year','FontSize',6);
    ylabel('Airborne fraction (unitless)','FontSize',5);
    if i == 1
        lgd = legend('Data','Linear trend','Linear trend+break','Location','NorthEast');
        lgd.FontSize = 4;
        legend('boxoff');
    end
    title(['Data: ',tit_str{i}],'FontSize',6);
    axis([t(1)-1,t(end)+1,0.2,0.8]);
    set(gca,'FontSize',5);
end

