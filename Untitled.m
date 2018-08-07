% Random seed for controlled experiment
for C = 1:1
    
    rng(C)
    % Generate X's
    n = 500; X = sort(10 * rand(n, 1));
    
    % True model (choose a kernel)
    covfunc = {@covSEiso};
%     covfunc = {@covMaterniso, 1};
%     covfunc = {@covMaterniso, 5};
    
    % Set parameters
    sigma = 1; noiseVariance = 0.001^2; signalVariance = 1;    
    hyp.cov = [log(sigma); log(signalVariance)];
    
    % Generate samples
    Ktmp = feval(covfunc{:}, hyp.cov, X);
    K = 0.5 * (Ktmp + Ktmp') + eye(n) * (noiseVariance);
    
    y = mvnrnd(zeros(n, 1), K); y = y(:);
    
    % Gaussian process regression
    [YPred, param] = GPWrapper(X, y(:), X);
    covfunc = param.covfunc;
    hyp = param.hyp;
    
    % Generate replicates from alternate model
    % covfunc = {@covSEiso};
    % hyp.cov = param.hyp.cov;
    % hyp.cov = [log(0.5); log(signalVariance)];
    
    % Generate replicate sample(s) from fitted model
    Ktmp = feval(covfunc{:}, hyp.cov, X);
    K = 0.5 * (Ktmp + Ktmp') + eye(n) * (exp(hyp.lik)^2);
    
    nRep = 1;
    y_rep = zeros(n, nRep);
    for rep = 1:nRep
        y_rep(:, rep) = mvnrnd(zeros(n, 1), K); y = y(:);
    end
    
    % PPC
    % countZeroCrossing = @(Y)(sum(abs(diff(bsxfun(@ge, Y, mean(Y))))));       
    % normFirstDer = @(X)(sum(diff(Y).^2)); 
    % normSecondDer = @(X)(sum(diff(diff(Y)).^2));
    % chiSquare = @(X)(diag(X' * inv(K) * X));
    % sctest = @(X)(mylbqtest(X, K));
        
    % find plug-in p-value analytically (for chi-square and LBQ)
    pvalList(C) = 1 - chi2cdf(y' * inv(K) * y, n);
%     [U, S] = eig(K); [~, pvalList(C)] = lbqtest(diag(diag(S).^-0.5) * U' * y);
    
    % find plug-in p-value using y_rep
    % testStatistic = chiSquare;
    % pvalList(C) = 2 * min(mean(testStatistic(y_rep) >  testStatistic(y)), mean(testStatistic(y_rep) <  testStatistic(y)));    
    fprintf('%d %0.6f\n', C, pvalList(C));
end
 
% TODO 
% Check implementation of two-sided p-value (OK http://onlinelibrary.wiley.com/doi/10.1002/bimj.4710320615/abstract)

function [YPred, hyp] = GPWrapper_MC(XTrain, YTrain, XTest, varargin)
% The function returns Bernoulli probabilities at testing points XTest
% given training samples XTrain (samples x features) and labels YTrain
%
% Author: Sohan Seth, sseth@inf.ed.ac.uk
 
YTrain = YTrain(:);
YTrain(YTrain == 0) = -1; % gp takes {+1, -1} as class labels
% normalize training and testing data
[XTrain, MU, SIGMA]  = zscore(XTrain);
XTest = bsxfun(@rdivide, bsxfun(@minus, XTest, MU), SIGMA);
 
meanfunc = @meanConst;
hyp.mean = mean(YTrain);
 
%covfunc = {@covMaterniso, 5};
%ell = median(pdist(XTrain)); % median intersample distance
%fprintf('using Matern kernel with median intersample distance\n')
%sf = 1; hyp.cov = log([ell; sf]);
 
covfunc = @covSEiso;
ell = median(pdist(XTrain)); % median intersample distance
fprintf('using Gaussian kernel with median intersample distance\n')
sf = 1; hyp.cov = log([ell; sf]);
 
%covfunc = @covSEard;
%ell = ones(size(XTrain, 2), 1);
%fprintf('using Gaussian kernel with identity covariance\n')
%sf = 1; hyp.cov = log([ell; sf]);
 
likfunc = @likErf;
 
hyp = minimize(hyp, @gp, -40, @infEP, meanfunc, covfunc, likfunc, XTrain, YTrain);
[~, ~, ~, ~, YPred] = gp(hyp, @infEP, meanfunc, covfunc, likfunc, XTrain, YTrain, XTest, ones(length(XTest), 1));
YPred = exp(YPred);
 
if nargout == 2
    hyp.covfunc = covfunc;
    hyp.meanfunc = meanfunc;
end
end