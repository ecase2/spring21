%%% macro semester 2, quarter 1 %%%
%%% problem set 1, question 3 %%%
%%% Emily Case Jan 30 %%%
%%% stuck :) 

% clean workspace
clc; clear

% set shooting method convergence variables
tol             = 1e-3;  % convergence tolerance
numcands        = 3e4;  % size of consumption grid increments
        %%% so danny uses consumption grid instead 
        %%% of capital 
traj_periods    = 500;   % number of periods to test for convergence to SS


%% Define parameters
sigma   = 1;
alpha   = 1/3;
beta    = 0.99^(1/12);
delta   = 0.01;
T       = 12;
D       = 1;

%% Calculate model's steady-state 

% capital, from Euler equation
k_ss = (((1/beta) - 1 + delta)/alpha)^(1/(alpha-1));

% consumption, from the resource constraint
c_ss = k_ss^alpha - delta*k_ss; % NOTE: assuming D=0 in the steady-state

%% shooting method %% 

%%% make a grid of consumption as "candidates" to try for the initial jump
%%% down from the steady state. jump in t=0, then from t=1 to T-1 follow
%%% trajectories from original stuff. then shock at time T, then back to
%%% trajectories for the original one. find the optimal c "jump" which
%%% brings it back to the steady state. 

%%% before making a grid, i want to try it with just one guess, and see if
%%% I can accurately code that before going through everything. 

c = 3.5; % temporary guess for initial jump
k=k_ss;
c_traj = ones(1,100);
k_traj = ones(1,100);
c_traj(1) = c;
k_traj(1) = k_ss;

% before the shock, after the initial jump:
for t=1:T-1
    cnext = c*beta^(-1/sigma)*(1-delta+alpha*k^(1-alpha))^(-1/sigma);
    knext = k^alpha -c +(1-delta)*k; % no D yet
    
    c_traj(t+1) = cnext;
    k_traj(t+1) = knext;
    
    % update:
    c = cnext;
    k=knext;
end

% for when the shock happens in T=12: 
    %%% note that i have to use T+1 bc trajectory started at 1 and not at
    %%% 0. 
c_traj(T+1) = c_traj(T)*beta^(1/sigma)*(1-delta+alpha*k_traj(T)^(1-alpha))^(-1/sigma);
k_traj(T+1) = k_traj(T)^alpha -c_traj(T) +(1-delta)*k_traj(T) -D;

% for after the shock:
c = c_traj(T+1);
k = k_traj(T+1);
for i=T+1:length(c_traj)
    
    cnext = c*beta^(1/sigma)*(1-delta+alpha*k^(1-alpha))^(-1/sigma);
    knext = k^alpha -c +(1-delta)*k; % no D 
    
    c_traj(i+1) = cnext;
    k_traj(i+1) = knext;
    
    % update:
    c = cnext;
    k=knext;
end

