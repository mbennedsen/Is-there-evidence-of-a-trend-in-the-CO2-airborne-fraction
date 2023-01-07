%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script will detect break points in the data using (1) the method
% proposed in the paper, where the sum of squared errors (SSE) are minimized,
% and (2) a profile likelihood approach, suggested by an anonymous referee, see
% footnote 1 in the Supplementary Information.
%
% Model 1: Break in intercept (level).
% Model 2: Break in trend (slope).
% Model 3: Break in both intercept and trend.
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

tau_start = 11; % Take out tau_start-1 time periods in the beginning and in the end of the sample to avoid end-point degeneracies.

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
SSE = nan(length(t),6,3);
logL = nan(length(t),6,3);
BIC = nan(length(t),6,3);
opt_break_point_likelihood = nan(6,3);
opt_break_point_sse = nan(6,3);
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

    T = length(t);
    
    for j = tau_start:(T-tau_start+1)
        br_date = t(j);

    
        %% Model 1: %% OLS  w. break in intercept: y = a1 + a2*I(t) + b*t 
        X = [ones(length(y),1),cumsum(t==br_date),t-t(1)];
        bhat = (X'*X)\X'*y;
        e_hat = y-X*bhat;
    
        SSE(j,i,1)  = sum(e_hat.^2);
        logL(j,i,1) = -0.5*length(t)*log(SSE(j,i,1)/length(t));
        BIC(j,i,1) = -2*logL(j,i,1) + length(bhat)*log(length(t));
           
        %% Model 2: OLS w. break in trend (not intercept): y = a + b1*t + b2*I(t)*(t-tau)
        X = [ones(length(y),1),t-t(1),cumsum(t==br_date).*(t-br_date+1)];
        bhat = (X'*X)\X'*y;
        e_hat = y-X*bhat;
        
        SSE(j,i,2)  = sum(e_hat.^2);
        logL(j,i,2) = -0.5*length(t)*log(SSE(j,i,2)/length(t));
        BIC(j,i,2) = -2*logL(j,i,2) + length(bhat)*log(length(t));
           
        %% Model 3: OLS w. break in intercept+trend: y = a1 + a2*I(t) + b1*t + b2*I(t)*(t-tau)
        X = [ones(length(y),1),cumsum(t==br_date),t-t(1),cumsum(t==br_date).*(t-br_date+1)];
        bhat = (X'*X)\X'*y;
        e_hat = y-X*bhat;
        
        SSE(j,i,3) = sum(e_hat.^2);
        logL(j,i,3) = -0.5*length(t)*log(SSE(j,i,3)/length(t));
        BIC(j,i,3) = -2*logL(j,i,3) + length(bhat)*log(length(t));
    
    end

    fig1 = figure(1);
    subplot(3,2,i);
    
    %% Find breakpoint
    for j = 1:3
        tmp = SSE(:,i,j);
        [val,indx] = min(tmp);
        opt_break_point_sse(i,j) = t(indx);

        tmp = logL(:,i,j);
        [val,indx] = max(tmp);
        opt_break_point_likelihood(i,j) = t(indx);


        plot(t,tmp,'LineWidth',2), hold on
    end
    title(['Profile likelihood: ',tit_str{i}],'FontSize',8);
    xlabel('Breakpoint (\tau)','FontSize',8);
    ylabel('Likelihood','FontSize',8);
    if i == 1
        legend({'Break in intercept','Break in trend','Break in trend+intercept'},'Fontsize',6,'Location','NorthWest');
        legend('boxoff');
    end


end

%% Print out results
disp(' ')
disp(' --- Optimal breakpoint based on SSE ---')
disp('       Model1      Model2      Model3')
disp(opt_break_point_sse);

disp(' ')
disp(' --- Optimal breakpoint based on profile likelihood --- ')
disp('       Model1      Model2      Model3')
disp(opt_break_point_likelihood);
