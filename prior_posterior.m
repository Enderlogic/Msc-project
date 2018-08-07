%% SE kernel
if strcmp(model_cov, 'SE')
    % length scale
    if prior.SE.lengthScale.sh < min(sample_hyp.cf{1, 1}.lengthScale)
        pxl = 1.2 * prior.SE.lengthScale.sh - 0.2 * max(sample_hyp.cf{1, 1}.lengthScale) : 1.2 * (max(sample_hyp.cf{1, 1}.lengthScale) - prior.SE.lengthScale.sh) / 100 : max(sample_hyp.cf{1, 1}.lengthScale);
    elseif prior.SE.lengthScale.sh > max(sample_hyp.cf{1, 1}.lengthScale)
        pxl = min(sample_hyp.cf{1, 1}.lengthScale) : 1.2 * (prior.SE.lengthScale.sh - min(sample_hyp.cf{1, 1}.lengthScale)) / 100 : 1.2 * prior.SE.lengthScale.sh - 0.2 * min(sample_hyp.cf{1, 1}.lengthScale);
    else
        pxl = min(sample_hyp.cf{1, 1}.lengthScale) : (max(sample_hyp.cf{1, 1}.lengthScale) - min(sample_hyp.cf{1, 1}.lengthScale)) / 100 : max(sample_hyp.cf{1, 1}.lengthScale);
    end
%     pyl = lognpdf(pxl, prior.SE.lengthScale.sh, prior.SE.lengthScale.s2);
    pyl = gampdf(pxl, prior.SE.lengthScale.sh, prior.SE.lengthScale.is);

    figure('rend','painters','pos',[10 10 500 400])
    [xl, ~] = histnorm(sample_hyp.cf{1, 1}.lengthScale, 40, 'plot');
    hold on
    plot(pxl, max(xl) / max(pyl) * pyl, 'r')
    vline(hyp_data.lengthScale, 'g', 'true value', 0.1)
    vline(prior.SE.lengthScale.sh, '--r', 'mean of prior', 0.2)
    vline(mean(sample_hyp.cf{1, 1}.lengthScale), '--', 'mean of posterior', 0.3)
    legend('posterior', 'prior')
    title(['Prior and posterior distribution of length scale in ' model_cov{1} ' kernel'])
    ylabel('Probability')
    
    % magnitude
    if prior.SE.magnSigma2.sh < min(sample_hyp.cf{1, 1}.magnSigma2)
        pxm = 1.2 * prior.SE.magnSigma2.sh - 0.2 * max(sample_hyp.cf{1, 1}.magnSigma2) : 1.2 * (max(sample_hyp.cf{1, 1}.magnSigma2) - prior.SE.magnSigma2.sh) / 400 : max(sample_hyp.cf{1, 1}.magnSigma2);
    elseif prior.SE.magnSigma2.sh > max(sample_hyp.cf{1, 1}.magnSigma2)
        pxm = min(sample_hyp.cf{1, 1}.magnSigma2) : 1.2 * (prior.SE.magnSigma2.sh - min(sample_hyp.cf{1, 1}.magnSigma2)) / 400 : 1.2 * prior.SE.magnSigma2.sh - 0.2 * min(sample_hyp.cf{1, 1}.magnSigma2);
    else
        pxm = min(sample_hyp.cf{1, 1}.magnSigma2) : (max(sample_hyp.cf{1, 1}.magnSigma2) - min(sample_hyp.cf{1, 1}.magnSigma2)) / 400 : max(sample_hyp.cf{1, 1}.magnSigma2);
    end
%     pym = lognpdf(pxm, prior.SE.magnSigma2.sh, prior.SE.magnSigma2.sh);
    pym = gampdf(pxm, prior.SE.magnSigma2.sh, prior.SE.magnSigma2.is);
    figure('rend','painters','pos',[10 10 500 400])
    [xm, ~] = histnorm(sample_hyp.cf{1, 1}.magnSigma2, 40, 'plot');
    hold on
    plot(pxm, max(xm) / max(pym) * pym, 'r')
    vline(hyp_data.magnSigma2, 'g', 'true value', 0.1)
    vline(prior.SE.magnSigma2.sh, '--r', 'mean of prior', 0.2)
    vline(mean(sample_hyp.cf{1, 1}.magnSigma2), '--', 'mean of posterior', 0.3)
    legend('posterior', 'prior')
    title(['Prior and posterior distribution of magnitude in ' strjoin(model_cov) ' kernel'])
    ylabel('Probability')
    
    % Sigma
    if prior.sigma2.sh < min(sample_hyp.lik.sigma2)
        pxs = 1.2 * prior.sigma2.sh - 0.2 * max(sample_hyp.lik.sigma2) : 1.2 * (max(sample_hyp.lik.sigma2) - prior.sigma2.sh) / 100 : max(sample_hyp.lik.sigma2);
    elseif prior.sigma2.sh > max(sample_hyp.lik.sigma2)
        pxs = min(sample_hyp.lik.sigma2) : 1.2 * (prior.sigma2.sh - min(sample_hyp.lik.sigma2)) / 100 : 1.2 * prior.sigma2.sh - 0.2 * min(sample_hyp.lik.sigma2);
    else
        pxs = min(sample_hyp.lik.sigma2) : (max(sample_hyp.lik.sigma2) - min(sample_hyp.lik.sigma2)) / 100 : max(sample_hyp.lik.sigma2);
    end
    figure('rend','painters','pos',[10 10 500 400])
    [xs, ~] = histnorm(sample_hyp.lik.sigma2, 40, 'plot');
%     pys = lognpdf(pxs, prior.sigma2.sh, prior.sigma2.sh);
    pys = gampdf(pxs, prior.sigma2.sh, prior.sigma2.is);
    hold on
    plot(pxs, max(xs) / max(pys) * pys, 'r')
    vline(hyp_data.sigma2, 'g', 'true value', 0.1)
    vline(prior.sigma2.sh, '--r', 'mean of prior', 0.2)
    vline(mean(sample_hyp.lik.sigma2), '--', 'mean of posterior', 0.3)
    legend('posterior', 'prior')
    title(['Prior and posterior distribution of 考 when using ' strjoin(model_cov) ' kernel'])
    ylabel('Probability')
elseif strcmp(model_cov, 'Periodic')
    % length scale
    if prior.Periodic.lengthScale.sh < min(sample_hyp.cf{1, 1}.lengthScale)
        pxl = 1.2 * prior.Periodic.lengthScale.sh - 0.2 * max(sample_hyp.cf{1, 1}.lengthScale) : 1.2 * (max(sample_hyp.cf{1, 1}.lengthScale) - prior.Periodic.lengthScale.sh) / 100 : max(sample_hyp.cf{1, 1}.lengthScale);
    elseif prior.Periodic.lengthScale.sh > max(sample_hyp.cf{1, 1}.lengthScale)
        pxl = min(sample_hyp.cf{1, 1}.lengthScale) : 1.2 * (prior.Periodic.lengthScale.sh - min(sample_hyp.cf{1, 1}.lengthScale)) / 100 : 1.2 * prior.Periodic.lengthScale.sh - 0.2 * min(sample_hyp.cf{1, 1}.lengthScale);
    else
        pxl = min(sample_hyp.cf{1, 1}.lengthScale) : (max(sample_hyp.cf{1, 1}.lengthScale) - min(sample_hyp.cf{1, 1}.lengthScale)) / 100 : max(sample_hyp.cf{1, 1}.lengthScale);
    end
%     pyl = lognpdf(pxl, prior.Periodic.lengthScale.sh, prior.Periodic.lengthScale.s2);
    pyl = gampdf(pxl, prior.Periodic.lengthScale.sh, prior.Periodic.lengthScale.is);

    figure('rend','painters','pos',[10 10 500 400])
    [xl, ~] = histnorm(sample_hyp.cf{1, 1}.lengthScale, 40, 'plot');
    hold on
    plot(pxl, max(xl) / max(pyl) * pyl, 'r')
    vline(hyp_data.lengthScale, 'g', 'true value', 0.2)
%     vline(prior.Periodic.lengthScale.sh, '--r', 'mean of prior', 0.2)
%     vline(mean(sample_hyp.cf{1, 1}.lengthScale), '--', 'mean of posterior', 0.3)
%     title(['Prior and posterior distribution of length scale in ' model_cov{1} ' kernel'])

    x_label = 'Length scale';
    y_label = 'Relative probability';
    xlabel(x_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
    ylabel(y_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
    legend({'posterior', 'prior'}, 'FontSize',15, 'Interpreter', 'latex', 'Location', 'best')
    legend boxoff
    box off
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    
    % magnitude
    if prior.Periodic.magnSigma2.sh < min(sample_hyp.cf{1, 1}.magnSigma2)
        pxm = 1.2 * prior.Periodic.magnSigma2.sh - 0.2 * max(sample_hyp.cf{1, 1}.magnSigma2) : 1.2 * (max(sample_hyp.cf{1, 1}.magnSigma2) - prior.Periodic.magnSigma2.sh) / 400 : max(sample_hyp.cf{1, 1}.magnSigma2);
    elseif prior.Periodic.magnSigma2.sh > max(sample_hyp.cf{1, 1}.magnSigma2)
        pxm = min(sample_hyp.cf{1, 1}.magnSigma2) : 1.2 * (prior.Periodic.magnSigma2.sh - min(sample_hyp.cf{1, 1}.magnSigma2)) / 400 : 1.2 * prior.Periodic.magnSigma2.sh - 0.2 * min(sample_hyp.cf{1, 1}.magnSigma2);
    else
        pxm = min(sample_hyp.cf{1, 1}.magnSigma2) : (max(sample_hyp.cf{1, 1}.magnSigma2) - min(sample_hyp.cf{1, 1}.magnSigma2)) / 400 : max(sample_hyp.cf{1, 1}.magnSigma2);
    end
%     pym = lognpdf(pxm, prior.Periodic.magnSigma2.sh, prior.Periodic.magnSigma2.sh);
    pym = gampdf(pxm, prior.Periodic.magnSigma2.sh, prior.Periodic.magnSigma2.is);
    figure('rend','painters','pos',[10 10 500 400])
    [xm, ~] = histnorm(sample_hyp.cf{1, 1}.magnSigma2, 40, 'plot');
    hold on
    plot(pxm, max(xm) / max(pym) * pym, 'r')
    vline(hyp_data.magnSigma2, 'g', 'true value', 0.1)
%     vline(prior.Periodic.magnSigma2.sh, '--r', 'mean of prior', 0.2)
%     vline(mean(sample_hyp.cf{1, 1}.magnSigma2), '--', 'mean of posterior', 0.3)
%     title(['Prior and posterior distribution of magnitude in ' strjoin(model_cov) ' kernel'])

    x_label = 'Magnitude';
    y_label = 'Relative probability';
    xlabel(x_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
    ylabel(y_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
    legend({'posterior', 'prior'}, 'FontSize',15, 'Interpreter', 'latex', 'Location', 'best')
    legend boxoff
    box off
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    
    % period
    if prior.Periodic.period.sh < min(sample_hyp.cf{1, 1}.period)
        pxp = 1.2 * prior.Periodic.period.sh - 0.2 * max(sample_hyp.cf{1, 1}.period) : 1.2 * (max(sample_hyp.cf{1, 1}.period) - prior.Periodic.period.sh) / 400 : max(sample_hyp.cf{1, 1}.period);
    elseif prior.Periodic.period.sh > max(sample_hyp.cf{1, 1}.period)
        pxp = min(sample_hyp.cf{1, 1}.period) : 1.2 * (prior.Periodic.period.sh - min(sample_hyp.cf{1, 1}.period)) / 400 : 1.2 * prior.Periodic.period.sh - 0.2 * min(sample_hyp.cf{1, 1}.period);
    else
        pxp = min(sample_hyp.cf{1, 1}.period) : (max(sample_hyp.cf{1, 1}.period) - min(sample_hyp.cf{1, 1}.period)) / 400 : max(sample_hyp.cf{1, 1}.period);
    end
%     pyp = lognpdf(pxp, log(prior.Periodic.period.sh), prior.Periodic.period.sh);
    pyp = gampdf(pxp, prior.Periodic.period.sh, prior.Periodic.period.is);
    figure('rend','painters','pos',[10 10 500 400])
    [xp, np] = histnorm(sample_hyp.cf{1, 1}.period, 200, 'plot');
    hold on
    plot(pxp, max(xp) / max(pyp) * pyp, 'r')
    vline(hyp_data.period, 'g', 'true value', 0.4)
%     vline(prior.Periodic.period.sh, '--r', 'mean of prior', 0.2)
%     vline(mean(sample_hyp.cf{1, 1}.period), '--', 'mean of posterior', 0.3)
%     title(['Prior and posterior distribution of period in ' model_cov{1} ' kernel'])
    
    x_label = 'Period';
    y_label = 'Relative probability';
    xlabel(x_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
    ylabel(y_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
    legend({'posterior', 'prior'}, 'FontSize',15, 'Interpreter', 'latex', 'Location', 'best')
    legend boxoff
    box off
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
    
    % Sigma
    if prior.sigma2.sh < min(sample_hyp.lik.sigma2)
        pxs = 1.2 * prior.sigma2.sh - 0.2 * max(sample_hyp.lik.sigma2) : 1.2 * (max(sample_hyp.lik.sigma2) - prior.sigma2.sh) / 100 : max(sample_hyp.lik.sigma2);
    elseif prior.sigma2.sh > max(sample_hyp.lik.sigma2)
        pxs = min(sample_hyp.lik.sigma2) : 1.2 * (prior.sigma2.sh - min(sample_hyp.lik.sigma2)) / 100 : 1.2 * prior.sigma2.sh - 0.2 * min(sample_hyp.lik.sigma2);
    else
        pxs = min(sample_hyp.lik.sigma2) : (max(sample_hyp.lik.sigma2) - min(sample_hyp.lik.sigma2)) / 100 : max(sample_hyp.lik.sigma2);
    end
    figure('rend','painters','pos',[10 10 500 400])
    [xs, ns] = histnorm(sample_hyp.lik.sigma2, 40, 'plot');
%     pys = lognpdf(pxs, prior.sigma2.sh, prior.sigma2.sh);
    pys = gampdf(pxs, prior.sigma2.sh, prior.sigma2.is);
    hold on
    plot(pxs, max(xs) / max(pys) * pys, 'r')
    vline(hyp_data.sigma2, 'g', 'true value', 0.1)
    vline(prior.sigma2.sh, '--r', 'mean of prior', 0.2)
    vline(mean(sample_hyp.lik.sigma2), '--', 'mean of posterior', 0.3)
%     title(['Prior and posterior distribution of 考 when using ' strjoin(model_cov) ' kernel'])
    
    x_label = 'Noise $\sigma$';
    y_label = 'Relative probability';
    xlabel(x_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
    ylabel(y_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
    legend({'posterior', 'prior'}, 'FontSize',15, 'Interpreter', 'latex', 'Location', 'best')
    legend boxoff
    box off
    set(findall(gca, 'Type', 'Line'),'LineWidth',2);
elseif strcmp(model_cov{1}, 'Matern')
    % length scale
    if prior.Matern.lengthScale.sh < min(sample_hyp.cf{1, 1}.lengthScale)
        pxl = 1.2 * prior.Matern.lengthScale.sh - 0.2 * max(sample_hyp.cf{1, 1}.lengthScale) : 1.2 * (max(sample_hyp.cf{1, 1}.lengthScale) - prior.Matern.lengthScale.sh) / 100 : max(sample_hyp.cf{1, 1}.lengthScale);
    elseif prior.Matern.lengthScale.sh > max(sample_hyp.cf{1, 1}.lengthScale)
        pxl = min(sample_hyp.cf{1, 1}.lengthScale) : 1.2 * (prior.Matern.lengthScale.sh - min(sample_hyp.cf{1, 1}.lengthScale)) / 100 : 1.2 * prior.Matern.lengthScale.sh - 0.2 * min(sample_hyp.cf{1, 1}.lengthScale);
    else
        pxl = min(sample_hyp.cf{1, 1}.lengthScale) : (max(sample_hyp.cf{1, 1}.lengthScale) - min(sample_hyp.cf{1, 1}.lengthScale)) / 100 : max(sample_hyp.cf{1, 1}.lengthScale);
    end
%     pyl = lognpdf(pxl, prior.Matern.lengthScale.sh, prior.Matern.lengthScale.s2);
    pyl = gampdf(pxl, prior.Matern.lengthScale.sh, prior.Matern.lengthScale.is);

    figure('rend','painters','pos',[10 10 500 400])
    [xl, ~] = histnorm(sample_hyp.cf{1, 1}.lengthScale, 40, 'plot');
    hold on
    plot(pxl, max(xl) / max(pyl) * pyl, 'r')
    vline(hyp_data.lengthScale, 'g', 'true value', 0.1)
    vline(prior.Matern.lengthScale.sh, '--r', 'mean of prior', 0.2)
    vline(mean(sample_hyp.cf{1, 1}.lengthScale), '--', 'mean of posterior', 0.3)
    legend('posterior', 'prior')
    title(['Prior and posterior distribution of length scale in ' model_cov{1} ' kernel'])
    ylabel('Probability')
    
    % magnitude
    if prior.Matern.magnSigma2.sh < min(sample_hyp.cf{1, 1}.magnSigma2)
        pxm = 1.2 * prior.Matern.magnSigma2.sh - 0.2 * max(sample_hyp.cf{1, 1}.magnSigma2) : 1.2 * (max(sample_hyp.cf{1, 1}.magnSigma2) - prior.Matern.magnSigma2.sh) / 400 : max(sample_hyp.cf{1, 1}.magnSigma2);
    elseif prior.Matern.magnSigma2.sh > max(sample_hyp.cf{1, 1}.magnSigma2)
        pxm = min(sample_hyp.cf{1, 1}.magnSigma2) : 1.2 * (prior.Matern.magnSigma2.sh - min(sample_hyp.cf{1, 1}.magnSigma2)) / 400 : 1.2 * prior.Matern.magnSigma2.sh - 0.2 * min(sample_hyp.cf{1, 1}.magnSigma2);
    else
        pxm = min(sample_hyp.cf{1, 1}.magnSigma2) : (max(sample_hyp.cf{1, 1}.magnSigma2) - min(sample_hyp.cf{1, 1}.magnSigma2)) / 400 : max(sample_hyp.cf{1, 1}.magnSigma2);
    end
%     pym = lognpdf(pxm, prior.Matern.magnSigma2.sh, prior.Matern.magnSigma2.sh);
    pym = gampdf(pxm, prior.Matern.magnSigma2.sh, prior.Matern.magnSigma2.is);
    figure('rend','painters','pos',[10 10 500 400])
    [xm, ~] = histnorm(sample_hyp.cf{1, 1}.magnSigma2, 40, 'plot');
    hold on
    plot(pxm, max(xm) / max(pym) * pym, 'r')
    vline(hyp_data.magnSigma2, 'g', 'true value', 0.1)
    vline(prior.Matern.magnSigma2.sh, '--r', 'mean of prior', 0.2)
    vline(mean(sample_hyp.cf{1, 1}.magnSigma2), '--', 'mean of posterior', 0.3)
    legend('posterior', 'prior')
    title(['Prior and posterior distribution of magnitude in ' strjoin(model_cov) ' kernel'])
    ylabel('Probability')
    
    % Sigma
    if prior.sigma2.sh < min(sample_hyp.lik.sigma2)
        pxs = 1.2 * prior.sigma2.sh - 0.2 * max(sample_hyp.lik.sigma2) : 1.2 * (max(sample_hyp.lik.sigma2) - prior.sigma2.sh) / 100 : max(sample_hyp.lik.sigma2);
    elseif prior.sigma2.sh > max(sample_hyp.lik.sigma2)
        pxs = min(sample_hyp.lik.sigma2) : 1.2 * (prior.sigma2.sh - min(sample_hyp.lik.sigma2)) / 100 : 1.2 * prior.sigma2.sh - 0.2 * min(sample_hyp.lik.sigma2);
    else
        pxs = min(sample_hyp.lik.sigma2) : (max(sample_hyp.lik.sigma2) - min(sample_hyp.lik.sigma2)) / 100 : max(sample_hyp.lik.sigma2);
    end
    figure('rend','painters','pos',[10 10 500 400])
    [xs, ~] = histnorm(sample_hyp.lik.sigma2, 40, 'plot');
%     pys = lognpdf(pxs, prior.sigma2.sh, prior.sigma2.sh);
    pys = gampdf(pxs, prior.sigma2.sh, prior.sigma2.is);
    hold on
    plot(pxs, max(xs) / max(pys) * pys, 'r')
    vline(hyp_data.sigma2, 'g', 'true value', 0.1)
    vline(prior.sigma2.sh, '--r', 'mean of prior', 0.2)
    vline(mean(sample_hyp.lik.sigma2), '--', 'mean of posterior', 0.3)
    legend('posterior', 'prior')
    title(['Prior and posterior distribution of 考 when using ' strjoin(model_cov) ' kernel'])
    ylabel('Probability')
elseif strcmp(model_cov, 'Exp')
    % length scale
    if prior.Exp.lengthScale.sh < min(sample_hyp.cf{1, 1}.lengthScale)
        pxl = 1.2 * prior.Exp.lengthScale.sh - 0.2 * max(sample_hyp.cf{1, 1}.lengthScale) : 1.2 * (max(sample_hyp.cf{1, 1}.lengthScale) - prior.Exp.lengthScale.sh) / 100 : max(sample_hyp.cf{1, 1}.lengthScale);
    elseif prior.Exp.lengthScale.sh > max(sample_hyp.cf{1, 1}.lengthScale)
        pxl = min(sample_hyp.cf{1, 1}.lengthScale) : 1.2 * (prior.Exp.lengthScale.sh - min(sample_hyp.cf{1, 1}.lengthScale)) / 100 : 1.2 * prior.Exp.lengthScale.sh - 0.2 * min(sample_hyp.cf{1, 1}.lengthScale);
    else
        pxl = min(sample_hyp.cf{1, 1}.lengthScale) : (max(sample_hyp.cf{1, 1}.lengthScale) - min(sample_hyp.cf{1, 1}.lengthScale)) / 100 : max(sample_hyp.cf{1, 1}.lengthScale);
    end
%     pyl = lognpdf(pxl, prior.Matern.lengthScale.sh, prior.Matern.lengthScale.s2);
    pyl = gampdf(pxl, prior.Exp.lengthScale.sh, prior.Exp.lengthScale.is);

    figure('rend','painters','pos',[10 10 500 400])
    [xl, ~] = histnorm(sample_hyp.cf{1, 1}.lengthScale, 40, 'plot');
    hold on
    plot(pxl, max(xl) / max(pyl) * pyl, 'r')
    vline(hyp_data.lengthScale, 'g', 'true value', 0.1)
    vline(prior.Exp.lengthScale.sh, '--r', 'mean of prior', 0.2)
    vline(mean(sample_hyp.cf{1, 1}.lengthScale), '--', 'mean of posterior', 0.3)
    legend('posterior', 'prior')
    title(['Prior and posterior distribution of length scale in ' model_cov{1} ' kernel'])
    ylabel('Probability')
    
    % magnitude
    if prior.Exp.magnSigma2.sh < min(sample_hyp.cf{1, 1}.magnSigma2)
        pxm = 1.2 * prior.Matern.magnSigma2.sh - 0.2 * max(sample_hyp.cf{1, 1}.magnSigma2) : 1.2 * (max(sample_hyp.cf{1, 1}.magnSigma2) - prior.Exp.magnSigma2.sh) / 400 : max(sample_hyp.cf{1, 1}.magnSigma2);
    elseif prior.Exp.magnSigma2.sh > max(sample_hyp.cf{1, 1}.magnSigma2)
        pxm = min(sample_hyp.cf{1, 1}.magnSigma2) : 1.2 * (prior.Exp.magnSigma2.sh - min(sample_hyp.cf{1, 1}.magnSigma2)) / 400 : 1.2 * prior.Exp.magnSigma2.sh - 0.2 * min(sample_hyp.cf{1, 1}.magnSigma2);
    else
        pxm = min(sample_hyp.cf{1, 1}.magnSigma2) : (max(sample_hyp.cf{1, 1}.magnSigma2) - min(sample_hyp.cf{1, 1}.magnSigma2)) / 400 : max(sample_hyp.cf{1, 1}.magnSigma2);
    end
%     pym = lognpdf(pxm, prior.Exp.magnSigma2.sh, prior.Exp.magnSigma2.sh);
    pym = gampdf(pxm, prior.Exp.magnSigma2.sh, prior.Exp.magnSigma2.is);
    figure('rend','painters','pos',[10 10 500 400])
    [xm, ~] = histnorm(sample_hyp.cf{1, 1}.magnSigma2, 40, 'plot');
    hold on
    plot(pxm, max(xm) / max(pym) * pym, 'r')
    vline(hyp_data.magnSigma2, 'g', 'true value', 0.1)
    vline(prior.Exp.magnSigma2.sh, '--r', 'mean of prior', 0.2)
    vline(mean(sample_hyp.cf{1, 1}.magnSigma2), '--', 'mean of posterior', 0.3)
    legend('posterior', 'prior')
    title(['Prior and posterior distribution of magnitude in ' strjoin(model_cov) ' kernel'])
    ylabel('Probability')
    
    % Sigma
    if prior.sigma2.sh < min(sample_hyp.lik.sigma2)
        pxs = 1.2 * prior.sigma2.sh - 0.2 * max(sample_hyp.lik.sigma2) : 1.2 * (max(sample_hyp.lik.sigma2) - prior.sigma2.sh) / 100 : max(sample_hyp.lik.sigma2);
    elseif prior.sigma2.sh > max(sample_hyp.lik.sigma2)
        pxs = min(sample_hyp.lik.sigma2) : 1.2 * (prior.sigma2.sh - min(sample_hyp.lik.sigma2)) / 100 : 1.2 * prior.sigma2.sh - 0.2 * min(sample_hyp.lik.sigma2);
    else
        pxs = min(sample_hyp.lik.sigma2) : (max(sample_hyp.lik.sigma2) - min(sample_hyp.lik.sigma2)) / 100 : max(sample_hyp.lik.sigma2);
    end
    figure('rend','painters','pos',[10 10 500 400])
    [xs, ~] = histnorm(sample_hyp.lik.sigma2, 40, 'plot');
%     pys = lognpdf(pxs, prior.sigma2.sh, prior.sigma2.sh);
    pys = gampdf(pxs, prior.sigma2.sh, prior.sigma2.is);
    hold on
    plot(pxs, max(xs) / max(pys) * pys, 'r')
    vline(hyp_data.sigma2, 'g', 'true value', 0.1)
    vline(prior.sigma2.sh, '--r', 'mean of prior', 0.2)
    vline(mean(sample_hyp.lik.sigma2), '--', 'mean of posterior', 0.3)
    legend('posterior', 'prior')
    title(['Prior and posterior distribution of 考 when using ' strjoin(model_cov) ' kernel'])
    ylabel('Probability')
elseif strcmp(model_cov, 'RQ')
    % length scale
    if prior.RQ.lengthScale.sh < min(sample_hyp.cf{1, 1}.lengthScale)
        pxl = 1.2 * prior.RQ.lengthScale.sh - 0.2 * max(sample_hyp.cf{1, 1}.lengthScale) : 1.2 * (max(sample_hyp.cf{1, 1}.lengthScale) - prior.RQ.lengthScale.sh) / 100 : max(sample_hyp.cf{1, 1}.lengthScale);
    elseif prior.RQ.lengthScale.sh > max(sample_hyp.cf{1, 1}.lengthScale)
        pxl = min(sample_hyp.cf{1, 1}.lengthScale) : 1.2 * (prior.RQ.lengthScale.sh - min(sample_hyp.cf{1, 1}.lengthScale)) / 100 : 1.2 * prior.RQ.lengthScale.sh - 0.2 * min(sample_hyp.cf{1, 1}.lengthScale);
    else
        pxl = min(sample_hyp.cf{1, 1}.lengthScale) : (max(sample_hyp.cf{1, 1}.lengthScale) - min(sample_hyp.cf{1, 1}.lengthScale)) / 100 : max(sample_hyp.cf{1, 1}.lengthScale);
    end
%     pyl = lognpdf(pxl, prior.Matern.lengthScale.sh, prior.Matern.lengthScale.s2);
    pyl = gampdf(pxl, prior.RQ.lengthScale.sh, prior.RQ.lengthScale.is);

    figure('rend','painters','pos',[10 10 500 400])
    [xl, ~] = histnorm(sample_hyp.cf{1, 1}.lengthScale, 40, 'plot');
    hold on
    plot(pxl, max(xl) / max(pyl) * pyl, 'r')
    vline(hyp_data.lengthScale, 'g', 'true value', 0.1)
    vline(prior.RQ.lengthScale.sh, '--r', 'mean of prior', 0.2)
    vline(mean(sample_hyp.cf{1, 1}.lengthScale), '--', 'mean of posterior', 0.3)
    legend('posterior', 'prior')
    title(['Prior and posterior distribution of length scale in ' model_cov{1} ' kernel'])
    ylabel('Probability')
    
    % magnitude
    if prior.RQ.magnSigma2.sh < min(sample_hyp.cf{1, 1}.magnSigma2)
        pxm = 1.2 * prior.RQ.magnSigma2.sh - 0.2 * max(sample_hyp.cf{1, 1}.magnSigma2) : 1.2 * (max(sample_hyp.cf{1, 1}.magnSigma2) - prior.RQ.magnSigma2.sh) / 400 : max(sample_hyp.cf{1, 1}.magnSigma2);
    elseif prior.RQ.magnSigma2.sh > max(sample_hyp.cf{1, 1}.magnSigma2)
        pxm = min(sample_hyp.cf{1, 1}.magnSigma2) : 1.2 * (prior.RQ.magnSigma2.sh - min(sample_hyp.cf{1, 1}.magnSigma2)) / 400 : 1.2 * prior.RQ.magnSigma2.sh - 0.2 * min(sample_hyp.cf{1, 1}.magnSigma2);
    else
        pxm = min(sample_hyp.cf{1, 1}.magnSigma2) : (max(sample_hyp.cf{1, 1}.magnSigma2) - min(sample_hyp.cf{1, 1}.magnSigma2)) / 400 : max(sample_hyp.cf{1, 1}.magnSigma2);
    end
%     pym = lognpdf(pxm, prior.Matern.magnSigma2.sh, prior.Matern.magnSigma2.sh);
    pym = gampdf(pxm, prior.RQ.magnSigma2.sh, prior.RQ.magnSigma2.is);
    figure('rend','painters','pos',[10 10 500 400])
    [xm, ~] = histnorm(sample_hyp.cf{1, 1}.magnSigma2, 40, 'plot');
    hold on
    plot(pxm, max(xm) / max(pym) * pym, 'r')
    vline(hyp_data.magnSigma2, 'g', 'true value', 0.1)
    vline(prior.RQ.magnSigma2.sh, '--r', 'mean of prior', 0.2)
    vline(mean(sample_hyp.cf{1, 1}.magnSigma2), '--', 'mean of posterior', 0.3)
    legend('posterior', 'prior')
    title(['Prior and posterior distribution of magnitude in ' strjoin(model_cov) ' kernel'])
    ylabel('Probability')
    
    % Sigma
    if prior.sigma2.sh < min(sample_hyp.lik.sigma2)
        pxs = 1.2 * prior.sigma2.sh - 0.2 * max(sample_hyp.lik.sigma2) : 1.2 * (max(sample_hyp.lik.sigma2) - prior.sigma2.sh) / 100 : max(sample_hyp.lik.sigma2);
    elseif prior.sigma2.sh > max(sample_hyp.lik.sigma2)
        pxs = min(sample_hyp.lik.sigma2) : 1.2 * (prior.sigma2.sh - min(sample_hyp.lik.sigma2)) / 100 : 1.2 * prior.sigma2.sh - 0.2 * min(sample_hyp.lik.sigma2);
    else
        pxs = min(sample_hyp.lik.sigma2) : (max(sample_hyp.lik.sigma2) - min(sample_hyp.lik.sigma2)) / 100 : max(sample_hyp.lik.sigma2);
    end
    figure('rend','painters','pos',[10 10 500 400])
    [xs, ~] = histnorm(sample_hyp.lik.sigma2, 40, 'plot');
%     pys = lognpdf(pxs, prior.sigma2.sh, prior.sigma2.sh);
    pys = gampdf(pxs, prior.sigma2.sh, prior.sigma2.is);
    hold on
    plot(pxs, max(xs) / max(pys) * pys, 'r')
    vline(hyp_data.sigma2, 'g', 'true value', 0.1)
    vline(prior.sigma2.sh, '--r', 'mean of prior', 0.2)
    vline(mean(sample_hyp.lik.sigma2), '--', 'mean of posterior', 0.3)
    legend('posterior', 'prior')
    title(['Prior and posterior distribution of 考 when using ' strjoin(model_cov) ' kernel'])
    ylabel('Probability')
    
    %alpha
    if prior.RQ.alpha.sh < min(sample_hyp.cf{1, 1}.alpha)
        pxa = 1.2 * prior.RQ.alpha.sh - 0.2 * max(sample_hyp.cf{1, 1}.alpha) : 1.2 * (max(sample_hyp.cf{1, 1}.alpha) - prior.RQ.alpha.sh) / 400 : max(sample_hyp.cf{1, 1}.alpha);
    elseif prior.RQ.alpha.sh > max(sample_hyp.cf{1, 1}.alpha)
        pxa = min(sample_hyp.cf{1, 1}.alpha) : 1.2 * (prior.RQ.alpha.sh - min(sample_hyp.cf{1, 1}.alpha)) / 400 : 1.2 * prior.RQ.alpha.sh - 0.2 * min(sample_hyp.cf{1, 1}.alpha);
    else
        pxa = min(sample_hyp.cf{1, 1}.alpha) : (max(sample_hyp.cf{1, 1}.alpha) - min(sample_hyp.cf{1, 1}.alpha)) / 400 : max(sample_hyp.cf{1, 1}.alpha);
    end
%     pya = lognpdf(pxa, log(prior.RQ.alpha.sh), prior.RQ.alpha.sh);
    pya = gampdf(pxa, prior.RQ.alpha.sh, prior.RQ.alpha.is);
    figure('rend','painters','pos',[10 10 500 400])
    [xa, na] = histnorm(log(sample_hyp.cf{1, 1}.alpha), 200, 'plot');
    hold on
    plot(pxa, max(xa) / max(pya) * pya, 'r')
    vline(hyp_data.alpha, 'g', 'true value', 0.1)
%     vline(prior.period.sh, '--r', 'mean of prior', 0.2)
    vline(mean(sample_hyp.cf{1, 1}.alpha), '--', 'mean of posterior', 0.3)
    legend('posterior', 'prior')
    title(['Prior and posterior distribution of alpha in ' model_cov{1} ' kernel'])
    ylabel('Probability')
end