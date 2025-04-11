%{
    Function to evaluate the negative log-likelihood as a function of the model parameters.
    Inputs:
    - choices: choices data
    - iX: data for index of state.
    - supportX: support of state vector.
    - capPi2: transition matrix for state variable
    - gamma: mean of the random shock
    - theta: parameters in the flow utility
    - delta: exponential discount factor
    - beta: present-bias discount factor 
    - betatilde: present-bias discount factor for naive agent
    - flowpayoffs: flowpayoffs function
    - CCP_generator: CCP_generator function 
    - backInduction: backInduction function

    Outputs:
    - nll: negative log-likelihood
%}
function [nll] = ...
         nLogLik(choices,iX,supportX,capPi2,theta,delta,beta,betatilde,gamma,flowpayoffs,CCP_generator,backInduction)
    nSuppX = size(supportX,1);
    nPeriods = size(iX,1);
    nFirms = size(iX,2);
    [u1,u2] = flowpayoffs(supportX,theta); 
            
    % Terminal period choice-specific value function.
    W_terminal = [u1,u2];
    % Generating CCP from logit model (type-1 extreme value).
    CCP_terminal = CCP_generator(W_terminal);  
    % Terminal period value function (NOT choice-specific, i.e., ex ante). 
    V_terminal = W_terminal(:,end) + gamma-log(CCP_terminal(:,end));
    
    % Backward induction to derive CCP for each period.
    % For exponential discounting case: beta = 1.
    [CCP_total,~,~] = backInduction(CCP_terminal,V_terminal,u1,u2,gamma,capPi2,beta,betatilde,delta,nPeriods);
    CCP_combined = CCP_total;
    p = ones(nPeriods,nFirms);
    
    % Calculate the negative log likelihood
    for t = 1:nPeriods
        CCP_t = CCP_combined(t*nSuppX-nSuppX+1:t*nSuppX,:);
        remaining_firms = sum(iX(t,:)~=0);
        for f = 1:remaining_firms
            p(t,f) = CCP_t(iX(t,f),choices(t,f));
        end
    end
    nll = -sum(sum(log(p)));
end