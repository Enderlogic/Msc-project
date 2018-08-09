clear
close all
%% SE
Min = -10;
Max = 10;
t = Min:0.1:Max;
sigma = 1;
lengthscale = 1;
SE = sigma^2*exp(-t.^2/2/lengthscale^2);
plot(t, SE)
lengthscale = 2;
SE = sigma^2*exp(-t.^2/2/lengthscale^2);
hold on
plot(t, SE)
x_label = '$\mathrm{x - y}$';
y_label = '$\mathrm{\kappa\left(x, y\right)}$';
xlabel(x_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
ylabel(y_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
legend({'\boldmath$\sigma = 1, \l = 1$', '\boldmath$\sigma = 1, \l = 2$'}, 'FontSize',15, 'Interpreter', 'latex', 'Location', 'best')
legend boxoff
box off
xlim([Min,Max]);
set(gca,'XTick',[Min : (Max - Min) / 2 : Max]);
ylim([0, 1]);
set(gca, 'YTick', [0 : 0.5 : 1]);

Min = 0;
Max = 15;
X = (Min : 0.1 : Max);
covfunc = {@covSEiso}; ell = 1; sf = 1; hyp.cov = log([ell; sf]);
K = feval(covfunc{:}, hyp.cov, X');
y1 = mvnrnd(zeros(1, length(X)), K);
x_label = '$\mathrm{x}$';
y_label = '$\mathrm{y}$';
figure
plot(X, y1)
ell = 2;
hyp.cov = log([ell; sf]);
K = feval(covfunc{:}, hyp.cov, X');
y2 = mvnrnd(zeros(1, length(X)), K);
hold on
plot(X, y2)
xlabel(x_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
ylabel(y_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
legend({'\boldmath$\sigma = 1, \l = 1$', '\boldmath$\sigma = 1, \l = 2$'}, 'FontSize',15, 'Interpreter', 'latex', 'Location', 'best')
legend boxoff
box off
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
xlim([Min,Max]);
set(gca,'XTick',[Min : (Max - Min) / 2 : Max]);
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
% X = (Min : 0.1 : Max)';
% covfunc = {'covPeriodic'}; ell = 3; sf = 1;p = 4; hyp.cov = log([ell;p;sf]);   % periodic
% K = feval(covfunc{:}, hyp.cov, X);
% y1 = mvnrnd(zeros(1, length(X)), K);
% x_label = '$x$';
% y_label = '$y$';
% figure
% plot(X, y1)
% p = 8;
% hyp.cov = log([ell;p;sf]);
% K = feval(covfunc{:}, hyp.cov, X);
% y2 = mvnrnd(zeros(1, length(X)), K);
% hold on
% plot(X, y2)
% xlabel(x_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% ylabel(y_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% legend({'\boldmath$\sigma = 1, \l = 3, T = 4$', '\boldmath$\sigma = 1, \l = 3, T = 8$'}, 'FontSize',15, 'Interpreter', 'latex', 'Location', 'best')
% legend boxoff
% box off
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% xlim([Min,Max]);
% set(gca,'XTick',[Min : (Max - Min) / 2 : Max]);

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
% covfunc = {@covMaterniso, 1}; ell = 5; sf = 1; hyp.cov = log([ell; sf]);
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
%% Composite

% Min = -40;
% Max = 40;
% t = Min:0.1:Max;
% sigma.SE = 1;
% lengthscale.SE = 8;
% SE = sigma.SE ^ 2 * exp(- t .^2 / 2 / lengthscale.SE ^ 2);
% sigma.Periodic = 1;
% lengthscale.Periodic = 3;
% T = 4;
% Periodic = sigma.Periodic ^ 2 * exp(- 2 * sin(pi * t / T) .^ 2/lengthscale.Periodic ^ 2);
% 
% plot(t, SE + Periodic)
% hold on
% plot(t, SE .* Periodic)
% x_label = '$x - y$';
% y_label = '$\kappa\left(x, y\right)$';
% xlabel(x_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% ylabel(y_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% legend({'$\mathrm{SE} + \mathrm{Periodic}$', '$\mathrm{SE} * \mathrm{Periodic}$'}, 'FontSize',15, 'Interpreter', 'latex', 'Location', 'best')
% legend boxoff
% box off
% xlim([Min,Max]);
% set(gca,'XTick',[Min : (Max - Min) / 2 : Max]);
% ylim([0, 2]);
% set(gca, 'YTick', [0 : 1 : 2]);
% 
% Min = 0;
% Max = 40;
% X = (Min : 0.1 : Max)';
% 
% covfunc = {'covSum',{@covSEiso, @covPeriodic}}; ell.SE = 8; sf.SE = 1;
% sf.Periodic = 1; ell.Periodic = 3; p = 4; cov_SE = log([ell.SE; sf.SE]); 
% cov_Periodic = log([ell.Periodic;p;sf.Periodic]);
% hyp.cov = [cov_SE; cov_Periodic];
% K = feval(covfunc{:}, hyp.cov, X);
% y1 = mvnrnd(zeros(1, length(X)), K);
% x_label = '$x$';
% y_label = '$y$';
% figure
% plot(X, y1)
% ell = 2;
% covfunc = {'covProd',{@covSEiso, @covPeriodic}};
% K = feval(covfunc{:}, hyp.cov, X);
% y2 = mvnrnd(zeros(1, length(X)), K);
% hold on
% plot(X, y2)
% xlabel(x_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% ylabel(y_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
% legend({'\boldmath$\sigma = 1, \l = 0.5$', '\boldmath$\sigma = 1, \l = 2$'}, 'FontSize',15, 'Interpreter', 'latex', 'Location', 'best')
% legend boxoff
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% xlim([Min,Max]);
% set(gca,'XTick',[Min : (Max - Min) / 2 : Max]);
% legend({'$\mathrm{SE} + \mathrm{Periodic}$', '$\mathrm{SE} * \mathrm{Periodic}$'}, 'FontSize',15, 'Interpreter', 'latex', 'Location', 'best')
% legend boxoff
% box off
% xlim([Min,Max]);
% set(gca,'XTick',[Min : (Max - Min) / 2 : Max]);
