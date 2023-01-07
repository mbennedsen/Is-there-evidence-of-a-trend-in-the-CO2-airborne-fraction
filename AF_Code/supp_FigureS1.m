%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script will produce Figure S1 of the Supplementary Information
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
% NB2: This script uses the built-in state space model (ssm) functions from
%      MATLAB.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('Data');
addpath(genpath('Functions'));
%% Init
filenam = 'Data/Marle_et_al_Nature_AirborneFraction_Datasheet.xlsx';

tit_str = {'GCP-raw','GCP-filter','H&N-raw','H&N-filter','New-raw','New-filter'};

options = optimset('TolX',1e-15,'TolFun',1e-15);
options = optimset(options,'MaxFunEvals',10000,'Display','off');

p0 = [0;0];      % Starting values for parameters in ML estimation.
yes_diffuse = 1; % =1: Diffuse initial of Kalman filter.
yes_refine = 1;  % =1 Use MATLAB's built-in function to find good starting values for parameters.
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
    
    %% Build state space model:
    if yes_diffuse == 1
        Mdl = dssm(@(params)helpfct_SS_v01(params,yes_diffuse));
    else
        Mdl = ssm(@(params)helpfct_SS_v01(params,yes_diffuse));
    end

    %% Estimate SS model
    if yes_refine == 1
        %% Find good initial values with refine function
        refine_output = refine(Mdl,y,p0);
        logL = cell2mat({refine_output.LogLikelihood})';
        [~,maxLogLIndx] = max(logL);
        refinedParams0 = refine_output(maxLogLIndx).Parameters;
        
        [EstMdl,estParams0,~,logL,output_est] = estimate(Mdl,y,refinedParams0,'Options',options); 
    else
        %[EstMdl,estParams0,EstParamCov,logL,output_est] = estimate(Mdl,y,params0,'Options',Opts);
        [EstMdl,estParams0,~,logL,output_est] = estimate(Mdl,y,p0,'Options',options);
    end

    %% Get estimate of trend using all data y1, y2, ..., yT.
    [mu_est,logL,Output_sm]   = smooth(EstMdl,y);

    %% Construct standard errors on mu_est
    stVal = 2; % (Lose one observation because of diffuse initialization)
    mu_est_cov = cell2mat({Output_sm.SmoothedStatesCov})';
    mu_est_cov = [nan(stVal-1,1);mu_est_cov];
    

    std_mu_est = nan(size(mu_est));
    %cov_smooth= nan(size(mu_est,2),size(mu_est,2),size(mu_est,1));
    for j = 1:size(mu_est,1)
        %cov_smooth(:,:,j) = mu_est_cov( ((j-1)*size(mu_est,2)+1):j*size(mu_est,2),: );
        std_mu_est(j,:) =  sqrt(mu_est_cov(j))' ;
    end

    figure(1);
    subplot(3,2,i)
    plot(t,y,'k-','LineWidth',2), hold on
    plot(t(stVal:end),mu_est(stVal:end),'r-','LineWidth',2), hold on
    plot(t(stVal:end),mu_est(stVal:end) + 1.96*std_mu_est(stVal:end),'r--','LineWidth',1), hold on
    plot(t(stVal:end),mu_est(stVal:end) - 1.96*std_mu_est(stVal:end),'r--','LineWidth',1), hold on
    
    ylabel('Airborne fraction (unitless)','FontSize',5);
    if i == 1
        lgd = legend('Data','Estimated trend','95p CB','Location','NorthWest');
        lgd.FontSize = 4;
        legend('boxoff');
    end
    title(['Data: ',tit_str{i}],'FontSize',6);
    %axis([t(1)-1,t(end)+1,0.2,0.8]);
    set(gca,'FontSize',5);
end

