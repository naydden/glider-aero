function createfigure(X1, Y1, llegenda, eixX, eixY, titol,nomArxiu)
    %CREATEFIGURE2(X1, Y1)
    %  X1:  vector of x data
    %  Y1:  vector of y data

    % Create figure
    figure1 = figure;

    % Create axes
    axes1 = axes('Parent',figure1,...
        'ZColor',[0.313725501298904 0.313725501298904 0.313725501298904],...
        'YMinorTick','on',...
        'YColor',[0.313725501298904 0.313725501298904 0.313725501298904],...
        'XMinorTick','on',...
        'XColor',[0.313725501298904 0.313725501298904 0.313725501298904]);
    box(axes1,'on');
    grid(axes1,'on');
    hold(axes1,'all');

    % Create plot
    plot(X1,Y1,'Parent',axes1,'LineWidth',2,'DisplayName',llegenda);

    % Create xlabel
    xlabel(eixX,'Interpreter','latex','FontSize',18);

    % Create ylabel
    ylabel(eixY,'Interpreter','latex','FontSize',18);

    % Create title
    title(titol,'Interpreter','latex','FontSize',24);

    % Create legend
    legend1 = legend(axes1,'show');
    set(legend1,'Interpreter','latex','Location','Best','LineWidth',1,...
        'FontSize',18);
    arxiu = strcat(nomArxiu,'.eps');
    path = './Report/plots/';
    print(strcat(path,arxiu),'-depsc2')
end
