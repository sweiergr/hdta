%{
    Function to simulate data from a set of CCPs.
    Inputs:
    - CCP: Conditional choice probabilities based on model parameters.
    - capPi2: Transition matrix for state variable
    - nPeriods: Number of time periods for which to simulate the model.
    - nFirms: Number of cross-section units to simulate ("sample size")
    - stream: Random stream for Matlab

    This function builds on code from "Dynamic Discrete Choice Models: Methods, 
    Matlab Code, and Exercises" by Jaap H. Abbring and Tobias J. Klein.
    See http://ddc.abbring.org for the original code and extensive documentation.

%}

function [choices,iX] = simulateData_CCP(CCP,capPi2,nPeriods,nFirms,stream)

    RandStream.setGlobalStream(stream);    
    nSuppX = size(capPi2,1);
    nA = size(CCP,2);
    % Set the initial state variable with equal probability to each value.
    p2Inf=ones(nSuppX,1)*1/nSuppX;
    % From the stationary distribution, get the sample of X1 (X for the
    % first period). Here, we assume that in period 0 all firms
    % choose action 2.
    iX = randomdiscrete(p2Inf*ones(1,nFirms));

    % Calculate choices in the first period from CCP.
    CCP_t = CCP(1:nSuppX,:);
    uniformDraws = ones(nA-1,1)*random('unif',zeros(1,nFirms),ones(1,nFirms));
    cumulativeCCP_t = cumsum(CCP_t,2);
    choices = sum([ones(1,nFirms);cumulativeCCP_t(iX,1:nA-1)'<=uniformDraws]);
    
    % Find the terminating "firms" and put them into the last columns.
    terminalcol = find(choices(end,:)==1);
    remainingcol = setdiff((1:nFirms),terminalcol);
    choices = choices(:,[remainingcol,terminalcol]);
    iX = iX(:,[remainingcol,terminalcol]);
    
    % Generate data for the following periods.
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