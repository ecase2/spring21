%% macro semester 2, quarter 1 %%
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

k = 170;

cmin = 3;
cmax = 4;
cguess = (cmin+cmax)/2;

cnext = cguess*beta^(1/sigma)*(1-delta+alpha*k^(1-alpha))^(-1/sigma)
knext = k^alpha -cguess + (1-delta)*k

% trying with very small grids:

% initialize grid 
kgrid = 160:1:200;
cgrid = 1:.1:5;

for i=1:length(kgrid)
    % consumption
    cgrid (1) = cguess
    cnext = cgrid(i)*beta^(1/sigma)*(1-delta+alpha*kgrid(i)^(1-alpha))^(-1/sigma)
    cgrid(i+1) = cnext
    
    kgrid(i+1) = kgrid(i)^alpha -cgrid(i) + (1-delta)*kgrid(i)
    
end
%%% this doesn't work and i do not know why