function [] = solutionplot(nPeriods,supportX,action,CCP_total,CCP_total_sophi,CCP_total_naive,hh)
    % plot many periods
    % action is the column to plot
    nSuppX = length(supportX);
    nXT = size(CCP_total,1);
    Ttemp = reshape(1:nXT,nSuppX,nPeriods)';

    for graphi = 1:nSuppX
        data_plot = [(1:nPeriods)', CCP_total(Ttemp(:,graphi),action),CCP_total_sophi(Ttemp(:,graphi),action),CCP_total_naive(Ttemp(:,graphi),action)];
        subplot(2,ceil(nSuppX/2),graphi)
        plot(data_plot(:,1),data_plot(:,2),'g',...  % g is consistent agent
            data_plot(:,1),data_plot(:,3),'b',...   % b is sophisticated agent
            data_plot(:,1),data_plot(:,4),'r');     % r is naive agent
        title("X=" + supportX(graphi));
        axis([-inf inf 0 1]);
    end
    lgd = legend(hh);
    set(lgd,...
    'Position',[0.257417727332612 0.00369974614018029 0.494004666846347 0.0479166663187917],...
    'Orientation','horizontal');
end