Files = dir('model\');
criteria = 'chi_square'; % choose one criteria for PPC or MMD ('number_of_zero'; 'norm_of_gradient'; 'chi_square'; 'mmd')
for k = 3 : length(Files)
    FileNames = Files(k).name;
    load(['model\' FileNames])
    p_value = p_value_calculator(data_test, data_test_rep, criteria, sample_type, model_cov, sample_hyp, x_train, data_train, x_test);
    figure('Position', [0,0,500,500])
    plot(x_test, data_test)
    hold on
    plot(x_test, data_test_rep(:,1:5))
    legend('observation', 'replication1', 'replication2', 'replication3', 'replication4', 'replication5', 'Location', 'best')
    if strcmp(sample_type, 'para')
        pd = 'p(y^{rep}|para)';
    elseif strcmp(sample_type, 'para&obs')
        pd =  'p(y^{rep}|para, y^{obs})';
    else
        error('The type of sampling is invalid')
    end
    title({['Generation model: ', strcat(strjoin(hyp_data.data_cov)), ' Fitted model: ', strcat(strjoin(model_cov))] ;...
        ['posterior distribution: ', pd] ; ['p-value: ', num2str(p_value)]});
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    fig = gcf;
    fig.PaperPositionMode = 'auto'
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    filename = strcat('model/figure/',strjoin(hyp_data.data_cov), '_', strjoin(model_cov), '_', sample_type, '.eps');
    print(fig, filename, '-dpdf')
    close all;
end