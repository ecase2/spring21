function [gdpa, gdpg, gdptauL, gdptauI] = counterfact(rhoI,rhoa,rhog,rhoL,a,g,...
    tauhatL,tauhatI,c,k,alpha,delta,sigma,phi,Gbar,Abar,taubarI,...
    taubarL,beta,Ybar,Kbar,Cbar,Lbar)


% recall iQ*X = Theta* Z, so we can create vectors:
v = [k c]; 
z = [a g tauhatL tauhatI];

% retrieve matrices A and B using LOM function and get Q and Lambda
% note that this is the same process as in BKFP.m, however now we have found the
% ACTUAL rhoI (where previously we used a guess, rhoI0)
    [A,B] = LOM(rhoI,rhoa,rhog,rhoL,alpha,delta,sigma,phi,Gbar,Abar,...
        taubarI,taubarL,beta,Ybar,Kbar,Cbar,Lbar);
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
    rho = eye(4) .* [ rhoa, rhog, rhoL, rhoI];

    % now use this e value in Theta
    Theta = -lam^(-1) * C(sel,:) * inv(eye(4) - lam^(-1) * rho ); 

% now new stuff... we want to look at gdp switching on and off each
% wedge/shock, so iterate through the four:

vi = v; % duplicates it so we can replace values 
l = 0*k; % initialize labor (length of the other ones)

for i=1:4
    switch i
        case 1
            z = [a 0*a 0*a 0*a]; % have to do 0*a bc a is a vector 
        case 2
            z = [0*a g 0*a 0*a];
        case 3
            z = [0*a 0*a tauhatL 0*a];
        case 4
            z = [0*a 0*a 0*a tauhatI];
    end
    
    for j = 1:length(k)
        
        % we can solve for c (see write up). note that this z has different
        % things switched on and off.
        vi(j,2) = (1/iQ(sel,2))*( -1*iQ(sel,1)*vi(j,1) + Theta*z(j,:)');
        
        % solve for labor:
        l(j) = ((-sigma)/(phi + alpha)) * vi(j,2) + ...
            ((alpha)/(phi + alpha))* vi(j,1) + ...
            ((-1*taubarL)/((alpha + phi)*(1-taubarL)))*z(j,3);
        
        if j<length(vi(:,1))
            % tomorrow's capital using today's capital, today's
            % consumption, and today's labor.
            vi(j+1,1) = (1-delta)*vi(j,1) + delta*(Ybar/(delta*Kbar))*(z(j,1)...
                + alpha*vi(j,1) + (1-alpha)*l(j) - (Cbar/Ybar)*vi(j,2) - Gbar*z(j,2));
        end
    end 
    
    % from a,k,l calculate y
    y = z(:,1) + alpha*vi(:,1) +(1-alpha)*l;
    switch i
        case 1
            gdpa = y;
        case 2
            gdpg = y;
        case 3
            gdptauL = y;
        case 4
            gdptauI = y;
    end
end 

end

