%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script will produce Figure S2 of the Supplementary Information
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

tau_start = 11;

p0 = [0;0];      % Starting values for parameters in ML estimation.
p00 = [0];
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

    sig = exp(estParams0(1:2));
    
    T = length(y);

    %% Estimate dummy coefficient: One for each period
    dummy_t = nan(T,1);
    for j = 1:T 
        %% Set up state space system
        [A,B,C,D,Mean0,Cov0,StateType] = helpfct_SS_wDummy_v01(sig,j,T,yes_diffuse);
        Mdl = ssm(A,B,C,D,'Mean0',Mean0,'Cov0',Cov0,'StateType',StateType);  

        %% GLS estimate of coefficient in front of dummy variable
        [mu_est2,logL,output] = filter(Mdl,y);
        dummy_est = cell2mat({output.ForecastedStates})';
        dummy_est = dummy_est(end,2);
        dummy_std = cell2mat({output.ForecastedStatesCov})';
        dummy_std = sqrt(dummy_std(end,2));

        dummy_t(j)   = dummy_est/dummy_std;
    end

    %% plot t-stats
    figure(1);
    subplot(3,2,i)
    stem(t,dummy_t), hold on
    plot(t,1.96*ones(length(t),1),'k-'), hold on
    plot(t,-1.96*ones(length(t),1),'k-'), hold on
    
    ylabel('t-stat','FontSize',5);
    title(['Data: ',tit_str{i}],'FontSize',6);
    %axis([t(1)-1,t(end)+1,0.2,0.8]);
    set(gca,'FontSize',5);
end

