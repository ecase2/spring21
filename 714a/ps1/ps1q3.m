%%% macro semester 2, quarter 1 %%%
%%% problem set 1, question 4 %%%
%%% Emily Case Jan 30 %%%
%%% (with lots of help and guidance from Danny Edgel)

%% TRYING AGAIN %%
% clean workspace
clc; clear

% set directory 
cd '/Users/emilycase/Desktop/spring21/714a'

% set shooting method convergence variables 
tol = 1e-3;  % convergence tolerance

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

%% Shooting method
% create a vector of candidates for the initial c jump
traj = 400; % trajectory length 
mindist = 1;

% BEFORE SHOCK
while mindist > 1e-2
    if mindist == 1
        min_c = 0; 
    else
        min_c = .01*c_ss + .99*min_c;
    end
    cposs = min_c:.0001:c_ss; % need to expand this later 
    cgrid = ones(length(cposs), traj);%%% each row is the trajectory for different 
                                      %%% c jump possibilities
    cgrid(:,1) = cposs; % enter in the first jumps

    % initialize the k trajectory similarly
    kgrid = ones(length(cposs), traj);
    kgrid(:,1) = k_ss; % because k won't change from initial jump

    for i = 1:(traj-1) %(T-1) % fills in next columns of grid from t=0 to t= T-1
        fprintf('\nLast dist = %d; filling in pre-shock, %d of %d...',mindist, i,traj-1)

        %%% note that 1 is already filled with the jump
        for j = 1:length(cposs) % goes through rows

            % capital before the shock:
            kgrid(j,i+1) = kgrid(j,i)^alpha - cgrid(j,i) + (1-delta)*kgrid(j,i);

            if i == T; kgrid(j,i+1) = kgrid(j,i+1) - D; end

            % consumption before the shock:
            %cgrid(j,i+1) = cgrid(j,i)*beta^(1/sigma)...
            %    *(alpha*kgrid(j,i+1)^(alpha-1) + 1-delta)^(1/sigma);
            cgrid(j,i+1) = (beta*(alpha*kgrid(j,i+1)^(alpha-1)+1-delta)) ...
                ^(1/sigma)*cgrid(j,i);
        end
    end
    % find distance between points and steady state:
    dist = ones(length(cposs),1);
    for i =1:length(cposs)

        dist(i) = norm(abs([kgrid(i, traj), cgrid(i, traj)] - [k_ss, c_ss]));
    end 

    [mindist, minind] = min(dist)
end

% adding in steady state:
k_path = [k_ss, kgrid(minind,:)];
c_path = [c_ss, cgrid(minind,:)];

%% PLOTS - obviously need to refine. how to do time on x axis?
figure(1)
    plot(0:traj,k_path); yline(k_ss,'--','Color','red')
    title('Capital Trajectory'); xlabel('t'); ylabel('K_t')
    saveas(gcf,'capital.png')
    
figure(2)
    plot(0:traj,c_path);  yline(c_ss,'--','Color','red')
    title('Consumption Trajectory'); xlabel('t'); ylabel('C_t')
    saveas(gcf,'consumption.png')

%% original attempt that didn't really work:
%{
% SHOCK 
for j = 1:length(cposs)
    
    % capital before at shock:
    kgrid(j,T+1) = kgrid(j,T)^alpha - cgrid(j,T) + (1-delta)*kgrid(j,T) - D;
    
    % consumption before at shock:
    %cgrid(j,T+1) = cgrid(j,T)*beta^(1/sigma)...
    %    *(1-delta + alpha*kgrid(j,T+1)^(alpha-1))^(1/sigma);
    
    cgrid(j,T+1) = (beta*(alpha*kgrid(j,T+1)^(alpha-1)+1-delta)) ...
        ^(1/sigma)*cgrid(j,T);
    
end

% AFTER SHOCK - THIS ONE TAKES FOREVER TO RUN?????
for i = (T+1):(traj-1) % fills in next columns of grid from t=13 to end
    %%% note that 1 is already filled with the jump
    fprintf('\nfilling in trajectory, %d of %d...',i,traj)
    
    for j = 1:length(cposs) % goes through rows
        
        knext = kgrid(j,i)^alpha - cgrid(j,i) + (1-delta)*kgrid(j,i);
        knext = kgrid(j,i)^alpha - cgrid(j,i) + (1-delta)*kgrid(j,i);
        %cnext = cgrid(j,i)*beta^(1/sigma)...
        %    *(1-delta + alpha*kgrid(j,i+1)^(alpha-1))^(1/sigma);
        
        cnext = (beta*(alpha*knext^(alpha-1)+1-delta))^(1/sigma)*cgrid(j,i);
        
        % check if the next period's values are real (code slows way down
        % if they're complex); if not, force convergence at last value
        if isreal(knext) && isreal(cnext)
            % capital after the shock:
            kgrid(j,i+1) = knext;

            % consumption after the shock:
            cgrid(j,i+1) = cnext;
        else
            kgrid(j,i+1) = kgrid(j,i);
            cgrid(j,i+1) = cgrid(j,i);
        end
        
        
    end
end
%}
