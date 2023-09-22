function [nll] = ...
         nLogLik(choices,iX,supportX,nA,capPi2,theta,delta,beta,betatilde,gamma,flowpayoffs,CCP_generator,backInduction)
    nSuppX = size(supportX,1);
    nPeriods = size(iX,1);
    nFirms = size(iX,2);
    [u1,u2] = flowpayoffs(supportX,theta); 
            
    % terminal period choice-specific value function:
    W_terminal = [u1,u2];
    % generating CCP by extreme value distribution
    CCP_terminal = CCP_generator(W_terminal);  
    % terminal period value function (not choice specific, i.e. ex ante). 
    V_terminal = W_terminal(:,end) + gamma-log(CCP_terminal(:,end));
    % backward induction to derive CCP for each period
    % for exponential case, beta = 1
    [CCP_total,~,~] = backInduction(CCP_terminal,V_terminal,u1,u2,gamma,capPi2,beta,betatilde,delta,nPeriods);
    
    CCP_combined = CCP_total;
    p = ones(nPeriods,nFirms);
    % calculate the negative log likelihood
    for t = 1:nPeriods
        CCP_t = CCP_combined(t*nSuppX-nSuppX+1:t*nSuppX,:);
        remaining_firms = sum(iX(t,:)~=0);
        for f = 1:remaining_firms
            p(t,f) = CCP_t(iX(t,f),choices(t,f));
        end
    end
    nll = -sum(sum(log(p)));
end