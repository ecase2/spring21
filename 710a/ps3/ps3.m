%{
ECON 710 (Quarter 3 metrics) 
Problem set 3 
Question 3 
Emily Case
with assistance from Danny Edgel
%}

clc, clear;



%% load in the data and clean it up 

data = readtable('AK91.csv');

% split up some variables into vectors:
Y = data.lwage;
X1 = data.educ;


% create dummy variables in X2:
x.sob = dummyvar(data.sob);
count1 = sum(x.sob);
x.sob = x.sob(:, count1>0);

x.yob = dummyvar(data.yob);
count2 = sum(x.yob);
x.yob = x.yob(:, count2>0);

% need to generate X_2 matrix of dummies (state birth and year birth). note
% that we have to omit the first dummies
X2 = [ones(size(X1)), x.yob(:,2:end), x.sob(:,2:end)];

X = [X1, X2];

% since qob is our instrument, we can also create our Z matrix:
x.qob = dummyvar(data.qob);
Z = [x.qob(:,2:4), X2];


%% calculate estiimator and its heteroskedasticity 

% estimate beta (from Hansen)
betahat = (X'*Z/(Z'*Z)*Z'*X)^(-1)*(X'*Z/(Z'*Z)*Z'*Y);

    
% estimate matrices in Hansen's book
n = length(X1);
Qzz = (1/n)*(Z'*Z);
Qxz = (1/n)*(X'*Z);
ehat = Y-X*betahat;
Ohat = 0*Qzz; % initialize 

for i=1:n % fill in the Ohats
    Ohat =  Ohat + (1/n)*Z(i,:)'*Z(i,:)*ehat(i)^2; 
end

% standard error
Vb = ((Qxz/Qzz*Qxz')\(Qxz/Qzz*Ohat/Qzz*Qxz')/(Qxz/Qzz*Qxz'))/n;
% have matlab show results
betahat(1) % should be 0.108
sqrt(Vb(1,1)) % should be 0.02



