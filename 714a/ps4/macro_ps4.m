%{ 
Macro quarter 3 
problem set 4 
Emily Case
%}

clear, clc;
cd '/Users/emilycase/Desktop/spring21/714a/ps4'

%% setting things up %% 

% parameters % 
rho = 1+ 1e-9; % dividing by 0 later will be an issue.
theta = 5;
Nk = 20;
K = 100000;
A = exp(normrnd(0,1,K,Nk));
W = 1; % initialize wage, we solve for this later

tol = 1e-6;

%% problem 4 %% 

% initial guess for sik, which was (Pik/Pk)^(1-theta), aka a ratio of
% prices
s0 = (1/Nk) + 0*A; 

maxiter = 1e10; 
iter = 1;
s = s0;
diff = 999;

while (diff>tol) && (iter<maxiter) 
    % solve pik from our recursive formula for producers' prices
    pik = (W./A).*(1 - 1./((1-theta) + s.*(theta-rho) ) );
    
    % solve pk from combining FOC with Ck definition
    pk = ( sum(pik.^(1-theta),2)).^(1/(1-theta));
    
    % use these to solve for new s 
    s = (pik./pk).^(1-theta); 
    
    % rho doesn't converge/the code takes forever so we create a new guess
    % which is a weighted average.
    tune = 0.6;
    s0 = s; % reestablish this 
    s = s0*tune + s*(1-tune);
    
    % check for convergence between s and s0
    diff = sum(sum(abs(s-s0)));
    
    % update iteration
    iter = iter+1;
end 

%% problem 5 %% 

p = ((1/K)*sum(pk.^(1-rho))).^(1/(1-rho));
rw = W/p; % real wage
C = rw; 
disp(['C = ' num2str(C)])

piks = W./A; % s is for social planners problems, aka first best
pks = (sum(piks.^(1-theta),2)).^(1/(1-theta)); 
p = ((1/K)*sum(pks.^(1-rho))).^(1/(1-rho));
rws = W/p;
Cs = rws;
disp(['Cs = ' num2str(Cs)])
