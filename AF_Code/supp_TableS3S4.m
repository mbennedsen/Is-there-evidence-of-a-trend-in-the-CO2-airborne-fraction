%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script will produce Tables S3 and S4 of the Supplementary Information
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
TableS3     = nan(12,6);
TableS4    = nan(9,6);
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
    
    %% OLS w. break in intercept+trend: y = a1 + a2*I(t) + b1*t + b2*I(t)*(t-tau+1)
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
    
    % a1
    TableS3(1,i) = bhat(1);
    TableS3(2,i) = sqrt(EstCov(1,1));
    TableS3(3,i) = bhat(1)/sqrt(EstCov(1,1));

    % b1
    TableS3(4,i) = bhat(3);
    TableS3(5,i) = sqrt(EstCov(3,3));
    TableS3(6,i) = bhat(3)/sqrt(EstCov(3,3));

    % a2
    TableS3(7,i) = bhat(2);
    TableS3(8,i) = sqrt(EstCov(2,2));
    TableS3(9,i) = bhat(2)/sqrt(EstCov(2,2));

    % b2
    TableS3(10,i) = bhat(4);
    TableS3(11,i) = sqrt(EstCov(4,4));
    TableS3(12,i) = bhat(4)/sqrt(EstCov(4,4));

    
    %% OLS w. break in intercept (not trend): y = a1 + a2*I(t) + b*t
    if mod(i,2) == 0
        br_date = 1990;
    else
        br_date = 1988;
    end
    X = [ones(length(y),1),cumsum(t==br_date),t-t(1)];
    bhat = (X'*X)\X'*y;
    e_hat = y-X*bhat;
    
    EstCov = hac(X,y,'display','off','intercept',false);
    se_HAC = sqrt(EstCov(2,2));
    
    % a1
    TableS4(1,i) = bhat(1);
    TableS4(2,i) = sqrt(EstCov(1,1));
    TableS4(3,i) = bhat(1)/sqrt(EstCov(1,1));

    % b1
    TableS4(4,i) = bhat(3);
    TableS4(5,i) = sqrt(EstCov(3,3));
    TableS4(6,i) = bhat(3)/sqrt(EstCov(3,3));

    % a2
    TableS4(7,i) = bhat(2);
    TableS4(8,i) = sqrt(EstCov(2,2));
    TableS4(9,i) = bhat(2)/sqrt(EstCov(2,2));


end

disp(' ');
disp(' Model with break in both intercept and slope: Table S3 in Supplementary Information')
disp(TableS3);
disp(' ');


disp(' Model with break in intercept only: Table S4 in Supplementary Information')
disp(TableS4);



