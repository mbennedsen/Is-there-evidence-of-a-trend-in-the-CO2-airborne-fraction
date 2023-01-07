function [A,B,C,D,Mean0,Cov0,StateType] = helpfct_SS_wDummy_v01(sig,indx,T,yesDiffuse)

%% Load in parameters
sigx1     = sig(1);%exp(params(1));
sigy1     = sig(2);%exp(params(2));
%dummyval  = params(1);

%% Set mean and covariance for initial states. 
StateType = [2;2]; % 0: Stationary; 1: unity; 2: non-stationary.

if yesDiffuse == 1
    Mean0 = [0;0]; % Mean of initial states
    Cov0  = blkdiag(Inf,Inf); % Cov. matrix of initial states
else
    Mean0 = [0;0]; % Mean of initial states
    Cov0  = blkdiag(1e5,1e5); % Cov. matrix of initial states
end
%% Construct transition matrices
A = [];
for i = 1:T
    if i == indx
        A = [A;{[1,1;0,1]}];
    else
        A = [A;{[1,0;0,1]}];
    end
end
%% Matrix C in obs. eq.
C = [1,0];

%% Construct error structure in state eq.
B =  [sigx1,0;
       0,0];

%% Construct error structure in obs eq.
D = sigy1;

%% The matrices not defined to be time-varying must be "forced" to be time-varying (Matlab idiosyncracity) 
%A = repmat({A},T,1);
B = repmat({B},T,1);
C = repmat({C},T,1);
D = repmat({D},T,1);
