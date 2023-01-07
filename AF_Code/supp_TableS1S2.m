%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script will produce Tables S1 and S2 of the Supplementary Information
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
Table_S1 = nan(6,6);
Table_S2  = nan(3,6);
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
    
    
    %% OLS for simple model: y = a + b*t
    % Test H0: b = 0
    if mod(i,2) == 0
        br_date = 1990;
    else
        br_date = 1988;
    end
    X = [ones(length(y),1),t-t(1)];
    bhat = (X'*X)\X'*y;
    e_hat = y-X*bhat;
    
    EstCov = hac(X,y,'display','off','intercept',false);
    
    % a
    Table_S1(1,i) = bhat(1);
    Table_S1(2,i) = sqrt(EstCov(1,1));
    Table_S1(3,i) = bhat(1)/sqrt(EstCov(1,1));

    % b
    Table_S1(4,i) = bhat(2);
    Table_S1(5,i) = sqrt(EstCov(2,2));
    Table_S1(6,i) = bhat(2)/sqrt(EstCov(2,2));
    z = Table_S1(6,i);

    Table_S2(1,i) = 2*normcdf(-abs(z)); % Two-sided
    Table_S2(2,i) = 1-normcdf(z); % One-sided, H1: b>0
    Table_S2(3,i) = normcdf(z); % One-sided, H1: b<0

end

%%
disp(' ');
disp(' Simple model y = a + b*t: Regression results (Table S1 in Supplementary Information)')
disp(Table_S1);
disp(' ');

disp(' Simple model y = a + b*t: p-values (Table S2 in Supplementary Information)')
disp(Table_S2');
