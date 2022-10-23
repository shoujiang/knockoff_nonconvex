%% demo of cv of nonconvex linear regression
%
clear;
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);
n = 500;
p = 100;
X = randn(n,p);             % generate a random design matrix
X = bsxfun(@rdivide, X, sqrt(sum(X.^2,1))); % normalize predictors
X = [ones(size(X,1),1) X];  % add intercept
b = zeros(p+1,1);           % true signal:
b(2:6) = 3;                 % first 5 predictors are 3
b(7:11) = -3;               % next 5 predictors are -3
y = X*b+randn(n,1);         % response vector

%%
% Sparse regression at a fixed tuning parameter value
penalty = 'SCAD';           % set penalty function to lasso
penparam = 3.7;
% penidx = ...                % leave intercept unpenalized
%     [false; true(size(X,2)-1,1)];
% lambdastart = ...           % find the maximum tuning parameter to start
%     max(lsq_maxlambda(sum(X(:,penidx).^2),-y'*X(:,penidx),penalty,penparam));
% display(lambdastart);

% lambda = 0.9*lambdastart;   % tuning parameter value

% Setting lambda %
lambda_max =norm(X'*y,'inf'); 
lambda_min = lambda_max * 0.001;
nlambda=10;
k = 5;
% for i=1:nlambda
%     Lambda1(i)=lambda_max*(lambda_min/lambda_max)^(i/nlambda);
%     lambda=Lambda1(i);
%     fprintf('iteration times:%d\n',i);
% end
%  i=1:nlambda;
% Lambda = lambda_max*(lambda_min/lambda_max).^(i./nlambda)
Lambda = 1:10;
opt_lambda = nonconvex_knockoffs.CV_NonConvex_LR(X,y,Lambda,k,penalty)

betahat = ...               
    lsq_sparsereg(X,y,opt_lambda,'penalty',penalty);

figure;                     % plot penalized estimate
bar(0:length(betahat)-1,betahat);
xlabel('j');
ylabel('\beta_j');
xlim([-1,length(betahat)]);
title([penalty '(' num2str(penparam) '), \lambda=' num2str(Lambda,2)]);