function [YPred, param] = GPWrapper(XTrain, YTrain, XTest, param)

% A wrapper for Gaussian Process regression
% Uses Gaussian Process Toolbox http://www.gaussianprocess.org/gpml/code/matlab/doc/
%
% Author: Sohan Seth, sseth@inf.ed.ac.uk
        
% normalize training and testing data
[XTrain, MU, SIGMA]  = zscore(XTrain);
XTest = bsxfun(@rdivide, bsxfun(@minus, XTest, MU), SIGMA);

meanfunc = @meanConst;
hyp.mean = mean(YTrain);

% meanfunc = {@meanSum, {@meanLinear, @meanConst}};
% hyp.mean = zeros(size(X, 2) + 1, 1);

covfunc = {@covSEiso};
ell = median(pdist(XTrain))^2; 
sf = 1; hyp.cov = log([ell; sf]);

% likfunc = @likGauss;
likfunc = @likLaplace; 
fprintf('likelihood %s\n', func2str(likfunc))
hyp.lik = log(0.01);

% inference = @infExact;
inference = @infVB;
fprintf('inference %s\n', func2str(inference))
%figure; imagesc(covMaterniso(3, hyp.cov, XTrain)), colorbar

rng default
hyp = minimize(hyp, @gp, -100, inference, meanfunc, covfunc, likfunc, XTrain, YTrain);
YPred = gp(hyp, inference, meanfunc, covfunc, likfunc, XTrain, YTrain, XTest);

if nargout == 2
    param.covfunc = covfunc;
    param.hyp = hyp;
end