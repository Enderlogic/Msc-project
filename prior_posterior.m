%% LengthScale
if prior.lengthScale.mu < min(sample_hyp.cf{1, 1}.lengthScale)
    pxl = 1.2 * prior.lengthScale.mu - 0.2 * max(sample_hyp.cf{1, 1}.lengthScale) : 1.2 * (max(sample_hyp.cf{1, 1}.lengthScale) - prior.lengthScale.mu) / 100 : max(sample_hyp.cf{1, 1}.lengthScale);
elseif prior.lengthScale.mu > max(sample_hyp.cf{1, 1}.lengthScale)
    pxl = min(sample_hyp.cf{1, 1}.lengthScale) : 1.2 * (prior.lengthScale.mu - min(sample_hyp.cf{1, 1}.lengthScale)) / 100 : 1.2 * prior.lengthScale.mu - 0.2 * min(sample_hyp.cf{1, 1}.lengthScale);
else
    pxl = min(sample_hyp.cf{1, 1}.lengthScale) : (max(sample_hyp.cf{1, 1}.lengthScale) - min(sample_hyp.cf{1, 1}.lengthScale)) / 100 : max(sample_hyp.cf{1, 1}.lengthScale);
end
pyl = lognpdf(pxl, prior.lengthScale.mu, prior.lengthScale.s2);
figure('rend','painters','pos',[10 10 500 400])
[xl, nl] = histnorm(sample_hyp.cf{1, 1}.lengthScale, 40, 'plot');
hold on
plot(pxl, max(xl) / max(pyl) * pyl, 'r')
vline(hyp_data.lengthScale, 'g', 'true value', 0.1)
% vline(prior.lengthScale.mu, '--r', 'mean of prior', 0.2)
vline(mean(sample_hyp.cf{1, 1}.lengthScale), '--', 'mean of posterior', 0.3)
legend('posterior', 'prior')
title(['Prior and posterior distribution of length scale in ' model_cov{1} ' kernel'])
ylabel('Probability')
%% Magnitude
if prior.magnSigma2.mu < min(sample_hyp.cf{1, 1}.magnSigma2)
    pxm = 1.2 * prior.magnSigma2.mu - 0.2 * max(sample_hyp.cf{1, 1}.magnSigma2) : 1.2 * (max(sample_hyp.cf{1, 1}.magnSigma2) - prior.magnSigma2.mu) / 400 : max(sample_hyp.cf{1, 1}.magnSigma2);
elseif prior.magnSigma2.mu > max(sample_hyp.cf{1, 1}.magnSigma2)
    pxm = min(sample_hyp.cf{1, 1}.magnSigma2) : 1.2 * (prior.magnSigma2.mu - min(sample_hyp.cf{1, 1}.magnSigma2)) / 400 : 1.2 * prior.magnSigma2.mu - 0.2 * min(sample_hyp.cf{1, 1}.magnSigma2);
else
    pxm = min(sample_hyp.cf{1, 1}.magnSigma2) : (max(sample_hyp.cf{1, 1}.magnSigma2) - min(sample_hyp.cf{1, 1}.magnSigma2)) / 400 : max(sample_hyp.cf{1, 1}.magnSigma2);
end
pym = lognpdf(pxm, log(prior.magnSigma2.mu), prior.magnSigma2.mu);
figure('rend','painters','pos',[10 10 500 400])
[xm, nm] = histnorm(sample_hyp.cf{1, 1}.magnSigma2, 40, 'plot');
hold on
plot(pxm, max(xm) / max(pym) * pym, 'r')
vline(hyp_data.magnitude, 'g', 'true value', 0.1)
% vline(prior.magnSigma2.mu, '--r', 'mean of prior', 0.2)
vline(mean(sample_hyp.cf{1, 1}.magnSigma2), '--', 'mean of posterior', 0.3)
legend('posterior', 'prior')
title(['Prior and posterior distribution of magnitude in ' model_cov{1} ' kernel'])
ylabel('Probability')

%% Sigma
if prior.sigma2.mu < min(sample_hyp.lik.sigma2)
    pxs = 1.2 * prior.sigma2.mu - 0.2 * max(sample_hyp.lik.sigma2) : 1.2 * (max(sample_hyp.lik.sigma2) - prior.sigma2.mu) / 100 : max(sample_hyp.lik.sigma2);
elseif prior.sigma2.mu > max(sample_hyp.lik.sigma2)
    pxs = min(sample_hyp.lik.sigma2) : 1.2 * (prior.sigma2.mu - min(sample_hyp.lik.sigma2)) / 100 : 1.2 * prior.sigma2.mu - 0.2 * min(sample_hyp.lik.sigma2);
else
    pxs = min(sample_hyp.lik.sigma2) : (max(sample_hyp.lik.sigma2) - min(sample_hyp.lik.sigma2)) / 100 : max(sample_hyp.lik.sigma2);
end
figure('rend','painters','pos',[10 10 500 400])
[xs, ns] = histnorm(sample_hyp.lik.sigma2, 40, 'plot');
pys = lognpdf(pxs, log(prior.sigma2.mu), prior.sigma2.mu);
hold on
plot(pxs, max(xs) / max(pys) * pys, 'r')
vline(hyp_data.sigma, 'g', 'true value', 0.1)
% vline(prior.sigma2.mu, '--r', 'mean of prior', 0.2)
vline(mean(sample_hyp.lik.sigma2), '--', 'mean of posterior', 0.3)
legend('posterior', 'prior')
title(['Prior and posterior distribution of ¦Ò when using ' model_cov{1} ' kernel'])
ylabel('Probability')

%% Period
if strcmp(model_cov{1}, 'Periodic')
    if prior.period.mu < min(sample_hyp.cf{1, 1}.period)
        pxp = 1.2 * prior.period.mu - 0.2 * max(sample_hyp.cf{1, 1}.period) : 1.2 * (max(sample_hyp.cf{1, 1}.period) - prior.period.mu) / 400 : max(sample_hyp.cf{1, 1}.period);
    elseif prior.period.mu > max(sample_hyp.cf{1, 1}.period)
        pxp = min(sample_hyp.cf{1, 1}.period) : 1.2 * (prior.period.mu - min(sample_hyp.cf{1, 1}.period)) / 400 : 1.2 * prior.period.mu - 0.2 * min(sample_hyp.cf{1, 1}.period);
    else
        pxp = min(sample_hyp.cf{1, 1}.period) : (max(sample_hyp.cf{1, 1}.period) - min(sample_hyp.cf{1, 1}.period)) / 400 : max(sample_hyp.cf{1, 1}.period);
    end
    pyp = lognpdf(pxp, log(prior.period.mu), prior.period.mu);
    figure('rend','painters','pos',[10 10 500 400])
    [xp, np] = histnorm(sample_hyp.cf{1, 1}.period, 200, 'plot');
    hold on
    plot(pxm, max(xp) / max(pyp) * pyp, 'r')
    vline(hyp_data.period, 'g', 'true value', 0.1)
%     vline(prior.period.mu, '--r', 'mean of prior', 0.2)
    vline(mean(sample_hyp.cf{1, 1}.period), '--', 'mean of posterior', 0.3)
    legend('posterior', 'prior')
    title(['Prior and posterior distribution of period in ' model_cov{1} ' kernel'])
    ylabel('Probability')
end