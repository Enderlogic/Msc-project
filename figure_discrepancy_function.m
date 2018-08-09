%% Number of zero crossing
% clear
close all
rng(990806);
%% Periodic
% Min = -10;
% Max = 10;
% t = Min:0.1:Max;
% sigma = 1;
% lengthscale = 3;
% T = 4;
% Periodic = sigma^2*exp(- 2 * sin(pi * t / T).^2/lengthscale^2);
% figure
% plot(t, Periodic)
% T = 8;
% Periodic = sigma^2*exp(- 2 * sin(pi * t / T).^2/lengthscale^2);
% hold on
% plot(t, Periodic)
% x_label = '$x - y$';
% y_label = '$\kappa\left(x, y\right)$';
% xlabel(x_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% ylabel(y_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% legend({'\boldmath$\sigma = 1, \l = 3, T = 4$', '\boldmath$\sigma = 1, \l = 3, T = 8$'}, 'FontSize',15, 'Interpreter', 'latex', 'Location', 'best')
% legend boxoff
% box off
% xlim([Min,Max]);
% set(gca,'XTick',[Min : (Max - Min) / 2 : Max]);
% ylim([0.7, 1]);
% set(gca, 'YTick', [0.7 : 0.15 : 1]);
% 
% Min = 0;
% Max = 20;
% X = (Min : 0.1 : Max);
% covfunc = {'covPeriodic'}; ell = 3; sf = 1;p = 4; hyp.cov = log([ell;p;sf]);   % periodic
% K = feval(covfunc{:}, hyp.cov, X');
% y1 = mvnrnd(zeros(1, length(X)), K);
% l1 = [X ; y1];
% l2 = [X ; mean(y1) * ones(1, length(y1))];
% inter1 = InterX(l1,l2);
% x_label = '$x$';
% y_label = '$y$';
% figure
% plot(X, y1, 'b')
% hline(mean(y1), ':b');
% hold on
% 
% ell = 1;
% hyp.cov = log([ell;p;sf]);
% K = feval(covfunc{:}, hyp.cov, X');
% y2 = mvnrnd(zeros(1, length(X)), K);
% l3 = [X ; y2];
% l4 = [X ; mean(y2) * ones(1, length(y2))];
% inter2 = InterX(l3,l4);
% plot(X, y2, 'r')
% hline(mean(y2), ':r');
% xlabel(x_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% ylabel(y_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% legend({'\boldmath$l = 3, T_\mathrm{noz} = 10$', '\boldmath$l = 1, T_\mathrm{noz} = 10$'}, 'FontSize',15, 'Interpreter', 'latex', 'Location', 'best')
% plot(inter1(1, :), inter1(2, :), 'rx');
% plot(inter2(1, :), inter2(2, :), 'bx');
% 
% legend boxoff
% box off
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% 
% xlim([Min,Max]);
% set(gca,'XTick',[Min : (Max - Min) / 2 : Max]);
% 
% gradient
% grad1 = gradient(y1, 0.1);
% grad2 = gradient(y2, 0.1);
% figure
% plot(X,grad1, 'b');
% hold on
% plot(X,grad2, 'r');
% y_label = '$\nabla y$';
% xlabel(x_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% ylabel(y_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% legend({'\boldmath$l = 3, T_\mathrm{nog} = 3.98$', '\boldmath$l = 1, T_\mathrm{nog} = 29.52$'}, 'FontSize',15, 'Interpreter', 'latex', 'Location', 'best')
% legend boxoff
% box off
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% xlim([Min,Max]);
% set(gca,'XTick',[Min : (Max - Min) / 2 : Max]);
% 
% ylim([-6, 3]);
% set(gca, 'YTick', [-6 : 1.5 : 3]);

%% SE
% Min = -10;
% Max = 10;
% t = Min:0.1:Max;
% sigma = 1;
% lengthscale = 1;
% SE = sigma^2*exp(- t .^ 2 / 2 /lengthscale^2);
% figure
% plot(t, SE)
% lengthscale = 2;
% SE = sigma^2*exp(- t .^ 2 / 2 /lengthscale^2);
% hold on
% plot(t, SE)
% x_label = '$x - y$';
% y_label = '$\kappa\left(x, y\right)$';
% xlabel(x_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% ylabel(y_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% legend({'\boldmath$\sigma = 1, l = 1 $', '\boldmath$\sigma = 1, l = 2 $'}, 'FontSize',15, 'Interpreter', 'latex', 'Location', 'best')
% legend boxoff
% box off
% xlim([Min,Max]);
% set(gca,'XTick',[Min : (Max - Min) / 2 : Max]);
% ylim([0, 1]);
% set(gca, 'YTick', [0 : 0.5 : 1]);

% Min = 0;
% Max = 15;
% X = (Min : 0.1 : Max);
% covfunc = {'covSEiso'}; ell = 1; sf = 1; hyp.cov = log([ell;sf]);   
% K = feval(covfunc{:}, hyp.cov, X');
% y1 = mvnrnd(zeros(1, length(X)), K + 0*eye(length(X)));
% l1 = [X ; y1];
% l2 = [X ; mean(y1) * ones(1, length(y1))];
% inter1 = InterX(l1,l2);
% x_label = '$x$';
% y_label = '$y$';
% figure
% plot(X, y1, 'b')
% hline(mean(y1), ':b');
% hold on
% 
% ell = 2;
% hyp.cov = log([ell;sf]);
% K = feval(covfunc{:}, hyp.cov, X');
% y2 = mvnrnd(zeros(1, length(X)), K + 0*eye(length(X)));
% l3 = [X ; y2];
% l4 = [X ; mean(y2) * ones(1, length(y2))];
% inter2 = InterX(l3,l4);
% plot(X, y2, 'r')
% hline(mean(y2), ':r');
% 
% xlabel(x_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% ylabel(y_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% legend({'\boldmath$l = 1, T_\mathrm{noz} = 4$', '\boldmath$l = 2, T_\mathrm{noz} = 1$'}, 'FontSize',15, 'Interpreter', 'latex', 'Location', 'best')
% plot(inter1(1, :), inter1(2, :), 'rx');
% plot(inter2(1, :), inter2(2, :), 'bx');
% legend boxoff
% box off
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% 
% xlim([Min,Max]);
% set(gca,'XTick',[Min : (Max - Min) / 2 : Max]);
% 
% % gradient
% grad1 = gradient(y1, 0.1);
% grad2 = gradient(y2, 0.1);
% figure
% plot(X,grad1, 'b');
% hold on
% plot(X,grad2, 'r');
% y_label = '$\nabla y$';
% xlabel(x_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% ylabel(y_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% legend({'\boldmath$l = 3, T_\mathrm{nog} = 3.98$', '\boldmath$l = 1, T_\mathrm{nog} = 29.52$'}, 'FontSize',15, 'Interpreter', 'latex', 'Location', 'best')
% legend boxoff
% box off
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% xlim([Min,Max]);
% set(gca,'XTick',[Min : (Max - Min) / 2 : Max]);
% 
% % ylim([-6, 3]);
% % set(gca, 'YTick', [-6 : 1.5 : 3]);
%% Matern
% Min = -15;
% Max = 15;
% t = Min:0.1:Max;
% sigma = 1;
% lengthscale = 5;
% Matern = sigma ^ 2 * exp(- abs(t) / lengthscale);
% figure
% plot(t, Matern)
% hold on
% Matern = sigma ^ 2 * (1 + sqrt(3) * abs(t) / lengthscale) .* exp(-sqrt(3) * abs(t) / lengthscale);
% plot(t, Matern)
% Matern = sigma ^ 2 * (1 + sqrt(5) * abs(t) / lengthscale + 5 * t .^ 2 / 3 / lengthscale ^ 2) .* exp(-sqrt(5) * abs(t) / lengthscale);
% plot(t, Matern)
% SE = sigma^2*exp(-t.^2/2/lengthscale^2);
% plot(t, SE)
% 
% x_label = '$x - y$';
% y_label = '$\kappa\left(x, y\right)$';
% xlabel(x_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% ylabel(y_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% legend({'\boldmath$\nu = 1/2$', '\boldmath$\nu = 3/2$', '\boldmath$\nu = 5/2$', '\boldmath$\nu\rightarrow\infty$'}, 'FontSize',15, 'Interpreter', 'latex', 'Location', 'best')
% legend boxoff
% box off
% xlim([Min,Max]);
% set(gca,'XTick',[Min : (Max - Min) / 2 : Max]);
% ylim([0, 1]);
% set(gca, 'YTick', [0 : 0.5 : 1]);
% 
% Min = 0;
% Max = 20;
% X = (Min : 0.1 : Max)';
% covfunc = {@covMaterniso, 1}; ell = 3; sf = 1; hyp.cov = log([ell; sf]);
% K = feval(covfunc{:}, hyp.cov, X);
% y1 = mvnrnd(zeros(1, length(X)), K);
% x_label = '$x$';
% y_label = '$y$';
% figure
% plot(X, y1)
% covfunc = {@covMaterniso, 3};
% hyp.cov = log([ell; sf]);
% K = feval(covfunc{:}, hyp.cov, X);
% y2 = mvnrnd(zeros(1, length(X)), K);
% hold on
% plot(X, y2)
% covfunc = {@covMaterniso, 5};
% hyp.cov = log([ell; sf]);
% K = feval(covfunc{:}, hyp.cov, X);
% y3 = mvnrnd(zeros(1, length(X)), K);
% plot(X, y3)
% covfunc = {@covSEiso};
% hyp.cov = log([ell; sf]);
% K = feval(covfunc{:}, hyp.cov, X);
% y4 = mvnrnd(zeros(1, length(X)), K);
% plot(X, y4)
% 
% xlabel(x_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% ylabel(y_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% legend({'\boldmath$\nu = 1/2$', '\boldmath$\nu = 3/2$', '\boldmath$\nu = 5/2$', '\boldmath$\nu\rightarrow\infty$'}, 'FontSize',15, 'Interpreter', 'latex', 'Location', 'best')
% legend boxoff
% box off
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% xlim([Min,Max]);
% set(gca,'XTick',[Min : (Max - Min) / 2 : Max]);

%% noise
% Min = 0;
% Max = 30;
% X = (Min : 0.1 : Max);
% covfunc = {'covSEiso'}; ell = 6; sf = 4; hyp.cov = log([ell;sf]);   
% K = feval(covfunc{:}, hyp.cov, X');
% y1 = mvnrnd(zeros(1, length(X)), K);
% y1n = mvnrnd(y1, 0.1*eye(length(X)));
% x_label = '$x$';
% y_label = '$y$';
% figure
% plot(X, y1, 'b')
% hold on
% 
% ell = 4;
% hyp.cov = log([ell;sf]);
% K = feval(covfunc{:}, hyp.cov, X');
% y2 = mvnrnd(zeros(1, length(X)), K);
% y2n = mvnrnd(y2, 0.1*eye(length(X)));
% plot(X, y2, 'r')
% 
% xlabel(x_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% ylabel(y_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% legend({'\boldmath$l = 6$', '\boldmath$l = 4$'}, 'FontSize',15, 'Interpreter', 'latex', 'Location', 'best')
% legend boxoff
% box off
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% 
% xlim([Min,Max]);
% set(gca,'XTick',[Min : (Max - Min) / 2 : Max]);
% 
% figure
% plot(X, y1n, 'b')
% hold on
% plot(X, y2n, 'r')
% 
% xlabel(x_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% ylabel(y_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% legend({'\boldmath$l = 6$', '\boldmath$l = 4$'}, 'FontSize',15, 'Interpreter', 'latex', 'Location', 'best')
% legend boxoff
% box off
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% 
% xlim([Min,Max]);
% set(gca,'XTick',[Min : (Max - Min) / 2 : Max]);
% 
% % gradient
% grad1 = gradient(y1, 0.1);
% grad2 = gradient(y2, 0.1);
% norm1 = norm(grad1);
% norm2 = norm(grad2);
% % figure
% % plot(X,grad1, 'b');
% % hold on
% % plot(X,grad2, 'r');
% % y_label = '$\nabla y$';
% % xlabel(x_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% % ylabel(y_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% % legend({'\boldmath$l = 6, T_\mathrm{nog} = 3.98$', '\boldmath$l = 4, T_\mathrm{nog} = 29.52$'}, 'FontSize',15, 'Interpreter', 'latex', 'Location', 'best')
% % legend boxoff
% % box off
% % set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% % xlim([Min,Max]);
% % set(gca,'XTick',[Min : (Max - Min) / 2 : Max]);
% 
% % ylim([-6, 3]);
% % set(gca, 'YTick', [-6 : 1.5 : 3]);
% 
% 
% grad1n = gradient(y1n, 0.1);
% grad2n = gradient(y2n, 0.1);
% norm1n = norm(grad1n);
% norm2n = norm(grad2n);
% figure
% plot(X,grad1n, 'b');
% hold on
% plot(X,grad2n, 'r');
% y_label = '$\nabla y$';
% xlabel(x_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% ylabel(y_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% legend({'\boldmath$l = 6, T_\mathrm{nog} = 40.63$', '\boldmath$l = 4, T_\mathrm{nog} = 39.29$'}, 'FontSize',15, 'Interpreter', 'latex', 'Location', 'best')
% legend boxoff
% box off
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% xlim([Min,Max]);
% set(gca,'XTick',[Min : (Max - Min) / 2 : Max]);
% 
% % ylim([-6, 3]);
% % set(gca, 'YTick', [-6 : 1.5 : 3]);
% beta = 0.9;
% par1 = (1 - beta) * fliplr(beta .^ (0 : length(y1n) - 1));
% par2 = beta .^ (1 - length(y1n) : 0);
% par3 = 1 - beta .^ (1 : length(y1n));
% 
% y1r = par1 .* y1n;
% y1r = cumsum(y1r);
% y1r = y1r .* par2;
% y1r = y1r ./ par3;
% 
% y2r = par1 .* y2n;
% y2r = cumsum(y2r);
% y2r = y2r .* par2;
% y2r = y2r ./ par3;
%                 
% grad1r = gradient(y1r, 0.1);
% grad2r = gradient(y2r, 0.1);
% norm1r = norm(grad1r);
% norm2r = norm(grad2r);
% figure
% plot(X,grad1r, 'b');
% hold on
% plot(X,grad1, '--b');
% hold on
% plot(X,grad2r, 'r');
% hold on
% plot(X,grad2, '--r');
% y_label = '$\nabla^\mathrm{mom} y$';
% xlabel(x_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% ylabel(y_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% legend({'\boldmath$l = 6, T_\mathrm{nog} = 11.02$', '\boldmath$l = 6, T_\mathrm{nog} = 10.77$', '\boldmath$l = 4, T_\mathrm{nog} = 15.65$', '\boldmath$l = 4, T_\mathrm{nog} = 15.85$'}, 'FontSize',15, 'Interpreter', 'latex', 'Location', 'best')
% legend boxoff
% box off
% set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
% xlim([Min,Max]);
% set(gca,'XTick',[Min : (Max - Min) / 2 : Max]);
% 
% % ylim([-6, 3]);
% % set(gca, 'YTick', [-6 : 1.5 : 3]);

%% chi square

Min = 12;
Max = 20;
% X = (Min : (Max-Min)/399 : Max);
X = x_test';
covfunc = {'covPeriodic'}; ell = 3; p=1; sf = 10; hyp.cov = log([ell;p;sqrt(sf)]);   
K = feval(covfunc{:}, hyp.cov, X');
K = (K + K') / 2;
% y1 = mvnrnd(zeros(1, length(X)), K + 0.01*eye(length(X)));
y1 = data_test';
% y1n = mvnrnd(y1, 0.01*eye(length(X)));
x_label = '$x$';
y_label = '$y$';
figure
plot(X, y1, 'b')
hold on

covfunc = {'covMaterniso',1};
ell = [18.1200137625339]; sf = [5.75511317283578];
hyp.cov = log([ell;sqrt(sf)]);
K = feval(covfunc{:}, hyp.cov, X');
K = (K + K') / 2;
y2 = mvnrnd(zeros(1, length(X)), K + ([1.83600810556985e-05])*eye(length(X)));
% y2n = mvnrnd(y2, 0.01*eye(length(X)));
plot(X, y2, 'r')

xlabel(x_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
ylabel(y_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% legend({'Observed data', 'Replicated data'}, 'FontSize',15, 'Interpreter', 'latex', 'Location', 'best')
legend boxoff
box off
set(findall(gca, 'Type', 'Line'),'LineWidth',2);

xlim([Min,Max]);
set(gca,'XTick',[Min : (Max - Min) / 2 : Max]);

% chi square
chi = y1 / (K+ ([1.83600810556985e-05])*eye(length(X))) * y1';
figure
plot(0:1:2*length(y1)-2, chi2pdf(0:1:2*length(y1)-2,length(y1)), 'b');
vline(chi, 'r')
1 - chi2cdf(chi, length(y1))
chi_f = floor(chi);
plotshaded([1:chi_f], [zeros(1, length([1:chi_f])); chi2pdf(1:chi_f,length(y1))], 'r')
x_label = '$D_{\chi^2}$';
y_label = '$Pr\left(D_{\chi^2}\right)$';
xlabel(x_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
ylabel(y_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% legend({'\boldmath$l = 6$', '\boldmath$l = 4$'}, 'FontSize',15, 'Interpreter', 'latex', 'Location', 'best')
legend boxoff
box off
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
