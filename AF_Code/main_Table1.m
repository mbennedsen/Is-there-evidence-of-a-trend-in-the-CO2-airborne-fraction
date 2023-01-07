%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script will run the hypothesis testing behind, and output the p-values 
% of, Table 1 of the main paper.
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
% NB2: The Mann-Kendall test is implemented using the "Mann_Kendall" function 
%      written by Simone Fatichi (2009). This code can be freely distributed, but 
%      please see the license for this code in "Functions/Mann_Kendall/license.txt.
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
        
    [h2,sig2]=Mann_Kendall(y,alpha);    

    pVals(i,1) = sig2;
    if i == 5 || i == 6 % For the "New" data, the slope is negative, so adjust the p-values to get the correct for MK- and MK+
        pVals(i,2) = 1-sig2/2;
        pVals(i,3) = sig2/2;
    else
        pVals(i,2) = sig2/2;
        pVals(i,3) = 1-sig2/2;
    end
    %% OLS (no breaks): y = a + b*t
    % Test: H0: b=0.
    X = [ones(length(y),1),t-t(1)];
    bhat = (X'*X)\X'*y;
    e_hat = y-X*bhat;
    s2 = e_hat'*e_hat/(length(e_hat) - length(bhat));
    SIG = s2*inv(X'*X);
    
    EstCov = hac(X,y,'display','off','intercept',false);
    se_HAC = sqrt(EstCov(2,2));
    
    pVals(i,4) = 2*normcdf(-abs(bhat(2))/se_HAC);
    

    %% OLS  w. break in intercept: y = a1 + a2*I(t) + b*t 
    % Test: H0: b=0.
    if mod(i,2) == 0
        br_date = 1990;
    else
        br_date = 1988;
    end
    X = [ones(length(y),1),cumsum(t==br_date),t-t(1)];
    bhat = (X'*X)\X'*y;
    e_hat = y-X*bhat;

    EstCov = hac(X,y,'display','off','intercept',false);
    se_HAC = sqrt(EstCov(3,3));

    pVals(i,5) = 2*normcdf(-abs(bhat(3))/se_HAC);

    
    %% OLS w. break in intercept+trend: y = a1 + a2*I(t) + b1*t + b2*I(t)*(t-tau+1)
    % Test H0: b2 = 0
    if mod(i,2) == 0
        br_date = 1990;
    else
        br_date = 1988;
    end
    X = [ones(length(y),1),cumsum(t==br_date),t-t(1),cumsum(t==br_date).*(t-br_date+1)];
    bhat = (X'*X)\X'*y;
    e_hat = y-X*bhat;
    
    EstCov = hac(X,y,'display','off','intercept',false);
    se_HAC = sqrt(EstCov(4,4));
    
    pVals(i,6) = 2*normcdf(-abs(bhat(4))/se_HAC);
    
    %% OLS w. break in intercept+trend: y = a1 + a2*I(t) + b1*t + b2*I(t)*(t-tau+1)
    % Test H0: a2 = 0
    if mod(i,2) == 0
        br_date = 1990;
    else
        br_date = 1988;
    end
    X = [ones(length(y),1),cumsum(t==br_date),t-t(1),cumsum(t==br_date).*(t-br_date+1)];
    bhat = (X'*X)\X'*y;
    e_hat = y-X*bhat;
    
    EstCov = hac(X,y,'display','off','intercept',false);
    se_HAC = sqrt(EstCov(2,2));
    
    pVals(i,7) = 2*normcdf(-abs(bhat(2))/se_HAC);
    

end

%% Print output to screen
disp(' ');
disp('Table 1 from the main paper (p-values):')
disp(pVals);
