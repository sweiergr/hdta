function [choices,iX] = simulateData_CCP(CCP,capPi2,nPeriods,nFirms)
    nSuppX = size(capPi2,1);
    nA = size(CCP,2);
    % Get the stationary distribution of the state variable x.
    % Since there are two transition matrices, there are two stationary distributions
    oneMinPi2 = eye(nSuppX)-capPi2';
    p2Inf = [oneMinPi2(1:nSuppX-1,:);ones(1,nSuppX)]\[zeros(nSuppX-1,1);1];
    % From the stationary distribution, get the sample of X1 (X for the
    % first period),here we internally assume that the period 0 all firms
    % choose action 2.
    iX = randomdiscrete(p2Inf*ones(1,nFirms));
    % Get the choices for the first period from CCP
    CCP_t = CCP(1:nSuppX,:);
    uniformDraws = ones(nA-1,1)*random('unif',zeros(1,nFirms),ones(1,nFirms));
    cumulativeCCP_t = cumsum(CCP_t,2);
    choices = sum([ones(1,nFirms);cumulativeCCP_t(iX,1:nA-1)'<=uniformDraws]);
    % Find the terminating firms and put them into the last columns.
    terminalcol = find(choices(end,:)==1);
    remainingcol = setdiff((1:nFirms),terminalcol);
    choices = choices(:,[remainingcol,terminalcol]);
    iX = iX(:,[remainingcol,terminalcol]);
    % Generate data for the following periods
    for t = 2:nPeriods 
        nFirms_remaining = length(remainingcol);
        transition_prob = [capPi2(iX(end,1:nFirms_remaining),:)]';
        temp_newstate = zeros(nSuppX,nFirms_remaining);
        for f = 1:nFirms_remaining
            row_index = choices(end,f);
            temp_newstate(:,f)=transition_prob(nSuppX*(row_index-2)+1:nSuppX*(row_index-2)+nSuppX,f);
        end
        iX_temp = randomdiscrete(temp_newstate);
        iX_temp = [iX_temp,zeros(1,nFirms-nFirms_remaining)];
        iX = [iX;iX_temp];
        
        CCP_t = CCP(t*nSuppX-nSuppX+1:t*nSuppX,:);
        temp_actionprob = zeros(nA,nFirms_remaining);
        for f = 1:nFirms_remaining
            temp_actionprob(:,f)=CCP_t(iX(end,f),:);
        end
        uniformDraws = ones(nA-1,1)*random('unif',zeros(1,nFirms_remaining),ones(1,nFirms_remaining));
        cumulativeP = cumsum(temp_actionprob);
        choice_temp = sum([ones(1,nFirms_remaining);cumulativeP(1:nA-1,:)<=uniformDraws]);
        choice_temp = [choice_temp,zeros(1,nFirms-nFirms_remaining)];
        choices = [choices;choice_temp];
        % Find the terminating firms and put them into the last columns.
        terminalcol = find(choices(end,:)==1);
        remainingcol = setdiff((1:nFirms_remaining),terminalcol);
        choices = choices(:,[remainingcol,terminalcol,nFirms_remaining+1:nFirms]);
        iX = iX(:,[remainingcol,terminalcol,nFirms_remaining+1:nFirms]);
    end
end