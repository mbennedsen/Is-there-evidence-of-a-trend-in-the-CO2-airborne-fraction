%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script will produce Table S5 of the Supplementary Information
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

SSE   = nan(6,3); 
logL  = nan(6,3); 
BIC   = nan(6,3);
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

    %% Model 1: OLS  w. break in intercept: y = a1 + a2*I(t) + b*t 
    if mod(i,2) == 0
        br_date = 1990;
    else
        br_date = 1988;
    end
    X = [ones(length(y),1),cumsum(t==br_date),t-t(1)];
    bhat = (X'*X)\X'*y;
    e_hat = y-X*bhat;

    SSE(i,1)  = sum(e_hat.^2);
    logL(i,1) = -0.5*length(t)*log(SSE(i,1)/length(t));
    BIC(i,1) = -2*logL(i,1) + length(bhat)*log(length(t));
     
    %% Model 2: OLS w. break in trend (not intercept): y = a + b1*t + b2*I(t)*(t-tau+1)
    % Test H0: b2=0
    if mod(i,2) == 0
        br_date = 1990;
    else
        br_date = 1988;
    end
    X = [ones(length(y),1),t-t(1),cumsum(t==br_date).*(t-br_date+1)];
    bhat = (X'*X)\X'*y;
    e_hat = y-X*bhat;

    SSE(i,2)  = sum(e_hat.^2);
    logL(i,2) = -0.5*length(t)*log(SSE(i,2)/length(t));
    BIC(i,2) = -2*logL(i,2) + length(bhat)*log(length(t));
    
    
    
    %% Model 3: OLS w. break in intercept+trend: y = a1 + a2*I(t) + b1*t + b2*I(t)*(t-tau+1)
    % Test H0: b2 = 0
    if mod(i,2) == 0
        br_date = 1990;
    else
        br_date = 1988;
    end
    X = [ones(length(y),1),cumsum(t==br_date),t-t(1),cumsum(t==br_date).*(t-br_date+1)];
    bhat = (X'*X)\X'*y;
    e_hat = y-X*bhat;


    SSE(i,3)  = sum(e_hat.^2);
    logL(i,3) = -0.5*length(t)*log(SSE(i,3)/length(t));
    BIC(i,3) = -2*logL(i,3) + length(bhat)*log(length(t));


end

%% Output Table:
disp(' ');
disp('Model comparison: Table S5 in Supplementary Information');

disp(' ');
disp(' log-likelihoods:')
disp(logL')

disp(' BIC:')
disp(BIC')

disp(' Likelihood ratio test statistics:')
disp(-2*(logL(:,1:2) - logL(:,3))')


