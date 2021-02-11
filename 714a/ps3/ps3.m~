%{
Macro quarter 3 problem set 3
Emily Case
%}

%% preparing data %% 
[x,xt] = xlsread('macropset3q3.xls','matlab');
[I,it] = xlsread('macropset3q3.xls','inv');

% add name vectors later % 

x = log(x); 
I = log(I);

xhp = hpfilter(x);
Ihp = hpfilter(I);

x= x-xhp;
I= I-Ihp;

x = x(:,2:end);% get rid of date column
I = I(:,2);

% initialize some stuff
delta = 0.025;
k = 0*I;

% assume the first capital is 0, which is already input
for t = 1:length(k)-1 % input k values 
    k(t+1) = (1-delta)*k(t) + delta * I(t); % linearized LOM
end


%% question 4 %% 

% more parameters % 
beta = 0.99;
alpha = 1/3;
sigma = 1;
phi = 1;
Abar = 1;
taubarL = 0;
taubarI = 0;
Gbar = 1/3;

% we can figure out a_t from the log lin 
% split up into individual vectors:
y = x(:,1);
c = x(:,2);
l = x(:,3);
T = length(y);
I = I(end-T+1:end);
k = k(end-T+1:end);

% solve for a_t
a = y - alpha*k -(1-alpha)*l;

% we need some bars for the next one which we do with solver 

% get steady state values:
[Ybar,Cbar,Kbar,Lbar] = calc_ss(alpha, sigma, phi, Gbar, Abar, taubarL, taubarI, beta, delta);
Ibar = delta*Kbar;

% solve for g_t
g = (y - (Cbar/Ybar)*c - (Ibar/Ybar)*I )/Gbar;

% solve for tau hat_L,t
tauhatL = alpha*k-alpha*l -phi*l -sigma*c;

% find persistence parameters %
rhoa = a(1:end-1)\a(2:end);
rhog = g(1:end-1)\g(2:end);
rhoL = tauhatL(1:end-1)\tauhatL(2:end);
rhoI0 = 0; % guess 

%% question 5 %% 
    % blanchard kahn method % 
    
[A,B] = LOM(rhoI0,rhoa,rhog,rhoL,alpha,delta,sigma,phi,Gbar,Abar,taubarI,taubarL,beta,Ybar,Kbar,Cbar,Lbar)
    




