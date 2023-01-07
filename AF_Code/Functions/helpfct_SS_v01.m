function [A,B,C,D,Mean0,Cov0,StateType] = helpfct_SS_v01(params,yesDiffuse)

%% Load in parameters
sigx1  = exp(params(1));
sigy1  = exp(params(2));

%% Set mean and covariance for initial states. 
StateType = [2]; % 0: Stationary; 1: unity; 2: non-stationary.

if yesDiffuse == 1
    Mean0 = 0; % Mean of initial states
    Cov0  = blkdiag(Inf); % Cov. matrix of initial states
else
    Mean0 = 0; % Mean of initial states
    Cov0  = blkdiag(1e5); % Cov. matrix of initial states
end
%% Construct transition matrices
A = 1;

%% Matrix C in obs. eq.
C = 1;

%% Construct error structure in state eq.
B =  sigx1;

%% Construct error structure in obs eq.
D = sigy1;


