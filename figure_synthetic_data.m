data_name = {'Matern', '5'};
hyp_data.magnSigma2 = 10; % squared magnitude
hyp_data.lengthScale = 3;
hyp_data.period = 2;
hyp_data.sigma2 = 0.01;
hyp_data.alpha = 5;
hyp_data.length = 20; % choose the length of data
hyp_data.random_seed = true; % choose random seed for data, if seed is true, the data will be the same; if seed is false, the data will be different every time
hyp_data.data_name = data_name; % choose a kernel for generation model ('SE'; 'Periodic'; 'Matern', '3' (nu = 3/2); 'Matern', '5' (nu = 5/2); 'RQ')
num_data = 1000; % choose sampling number

[x, data] = data_generation(hyp_data, num_data);

num_data_train = round(num_data * ratio_train);
data_train = data(1 : num_data_train);
x_train = x(1 : num_data_train);
data_test = data(num_data_train + 1 : num_data);
x_test = x(num_data_train + 1 : num_data);
disp('Data generation complete!')

figure
plot(x_train, data_train, 'b')
hold on
plot(x_test, data_test, 'r')
x_label = '$x$';
y_label = '$y$';
xlabel(x_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
ylabel(y_label, 'FontSize',25, 'FontWeight','bold', 'Interpreter', 'latex')
legend({'train', 'test'}, 'FontSize',15, 'Interpreter', 'latex', 'Location', 'best')
legend boxoff
box off
set(findall(gca, 'Type', 'Line'),'LineWidth',2);

% ylim([-7.5, -3]);
% set(gca, 'YTick', [-7.5 : 1.5 : -3]);