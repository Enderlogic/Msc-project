Files = dir('model\');
% criteria = 'chi_square'; % choose one criteria for PPC or MMD ('number_of_zero'; 'norm_of_gradient'; 'chi_square'; 'mmd')
for k = 4 : length(Files)
    FileNames = Files(k).name;
    load(['model\' FileNames])
    criteria = 'chi_square'; % choose one criteria for PPC or MMD ('number_of_zero'; 'norm_of_gradient'; 'chi_square'; 'mmd')
    p_value_chi = p_value_calculator(data_test, data_rep, criteria, model_cov, sample_location, sample_hyp, x_train, data_train, x_test);
    criteria = 'norm_of_gradient'; % choose one criteria for PPC or MMD ('number_of_zero'; 'norm_of_gradient'; 'chi_square'; 'mmd')
    p_value_gra = p_value_calculator(data_test, data_rep, criteria, model_cov, sample_location, sample_hyp, x_train, data_train, x_test);
    if strcmp(sample_type, 'para') || strcmp(sample_type, 'both')
        figure('Position', [0,0,500,500])
        plot(x_test, data_test)
        hold on
        plot(x_test, data_rep.para(:,1:5), '.')
        pd = 'p(y^{rep}|para)';
        legend('observation', 'replication1', 'replication2', 'replication3', 'replication4', 'replication5', 'Location', 'best')
        title({['Generation model: ', strcat(strjoin(data_name)), ' Fitted model: ', strcat(strjoin(model_cov))] ;...
            ['posterior distribution: ', pd] ; ['p-value: chi square: ', num2str(p_value_chi.para), '    norm of gradient: ', num2str(p_value_gra.para)]});
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
        filename = strcat('model/figure/',strjoin(data_name), '_', strjoin(model_cov), '_', sample_location, '_para.eps');
        print(fig, filename, '-dpdf')
        hold off
        close all;
    end
    if strcmp(sample_type, 'para&obs') || strcmp(sample_type, 'both')
        figure('Position', [0,0,500,500])
        plot(x_test, data_test)
        hold on
        plot(x_test, data_rep.para_obs(:,1:5), '.')
        pd =  'p(y^{rep}|para, y^{obs})';
        legend('observation', 'replication1', 'replication2', 'replication3', 'replication4', 'replication5', 'Location', 'best')
        title({['Generation model: ', strcat(strjoin(data_name)), ' Fitted model: ', strcat(strjoin(model_cov))] ;...
            ['posterior distribution: ', pd] ; ['p-value: chi square: ', num2str(p_value_chi.para_obs), '    norm of gradient: ', num2str(p_value_gra.para_obs)]});
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
        filename = strcat('model/figure/',strjoin(data_name), '_', strjoin(model_cov), '_', sample_location, '_para&obs.eps');
        print(fig, filename, '-dpdf')
        hold off
        close all;
    end
end