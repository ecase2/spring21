%{
This file is for macro quarter 3 problem set 3, blanchard kahn method.
(for problem 5)
%}

% main function:
function [rhoI,tauhatI] = BKFP(rhoI0,rhoa,rhog,rhoL,a,g,tauhatL,c,k,alpha,delta,sigma,phi ...
        ,Gbar,Abar,taubarI,taubarL,beta,Ybar,Kbar,Cbar,Lbar);
% Solves for fixed point solution to htauT by continually applying
% Blanchard-Kahn.

% initialize
diff = 999;
maxiter = 1e6;
iter = 1;

while (iter<maxiter)&&(diff>1e-6) 
    
    % get values using BK_calcrho function (below)
    [rhoI,tauhatI] = BK_calcrho(rhoI0,rhoa,rhog,rhoL,a,g,tauhatL,c,k,alpha,delta,sigma,phi ...
        ,Gbar,Abar,taubarI,taubarL,beta,Ybar,Kbar,Cbar,Lbar);
    
    % update:
    diff = abs(rhoI0 - rhoI); % how different is rhoI from the guess we made?
    iter = iter + 1;
    rhoI0 = rhoI; 
end
end



function tauhatI = BK_tauhatI(rhoI0,rhoa,rhog,rhoL,a,g,tauhatL,c,k,alpha,...
    delta,sigma,phi ,Gbar,Abar,taubarI,taubarL,beta,Ybar,Kbar,Cbar,Lbar)
% local function to calculate tauhatI

% retrieve matrices A and B using LOM function and get Q and Lambda
[A,B] = LOM(rhoI0,rhoa,rhog,rhoL,alpha,delta,sigma,phi,Gbar,Abar,taubarI,taubarL,beta,Ybar,Kbar,Cbar,Lbar);
[Q, Lambda] = eig(A);
iQ = inv(Q);

% we want to retrieve the eigenvalue that is greater than 1, the
% "explosive" one
if abs(Lambda(1,1))>1;
    sel = 1; % then the first one is >1
else 
    sel = 2; % the second one must be >1
end 

% define matrix C as in typed up answers
C = iQ*B;
lm = diag(Lambda); % vector of e values instead of matrix
lam = lm(sel); % pulls out the e value that is >1

% create matrix with shocks on the diagonal
rho = eye(4) .* [ rhoa, rhog, rhoL, rhoI0];

% now use this e value in Theta
Theta = -lam^(-1) * C(sel,:) * inv(eye(4) - lam^(-1) * rho ); 

% iQ*X = Theta* Z, so we can now solve for tauhat:
tauhatI = ([k, c]*iQ(sel,:)' -([a, g, tauhatL]*Theta(1:3)'))/Theta(4);
end




function [rhoI,tauhatI] = BK_calcrho(rhoI0,rhoa,rhog,rhoL,a,g,tauhatL,c,k,alpha,delta,sigma,phi ...
        ,Gbar,Abar,taubarI,taubarL,beta,Ybar,Kbar,Cbar,Lbar)
% this mostly just makes the top global function work a little bit
% smoother.
    
    % use the previous function to calculate tauhatI for a specific guess
    tauhatI = BK_tauhatI(rhoI0,rhoa,rhog,rhoL,a,g,tauhatL,c,k,alpha,...
    delta,sigma,phi ,Gbar,Abar,taubarI,taubarL,beta,Ybar,Kbar,Cbar,Lbar);
    
    % OLS regression to get rhoI
    rhoI = tauhatI(1:end-1)\tauhatI(2:end); 
end

