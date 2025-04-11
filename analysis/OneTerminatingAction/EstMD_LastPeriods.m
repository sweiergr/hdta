%{
    Function to calculate objective function for minimum distance estimation using final periods of data.
    Inputs:
    - nPeriods: number of time periods in the data.
    - supportX: support of state vector.
    - capPi2: transition matrix for state variable
    - beta: present-bias discount factor
    - delta: exponential discount factor
    - gamma: mean of the random shock
    - CCP: CCPs as estimated from data

%}
function [z, u1, SSres]=EstMD_LastPeriods(nPeriods,supportX,capPi2,beta,delta,gamma,CCP)
    
    % Estimation without incorporating exclusion restrictions.
    nSuppX = length(supportX);
    Q_t_pb = repmat(CCP(:,2),1,nSuppX).*repmat(capPi2,nPeriods,1);
    betadelta = beta*delta;
    phi_21 = log(CCP(:,2)) - log(CCP(:,1));
    m1 = gamma - log(CCP(:,1));
    phi_wide = reshape(phi_21,[nSuppX,nPeriods]);
    m1_wide = reshape(m1,[nSuppX,nPeriods]);
    logCCP_wide=reshape(log(CCP(:,1)),[nSuppX,nPeriods]);
    T = nPeriods;
    A = capPi2*(m1_wide(:,T-1) - m1_wide(:,T));
    B = capPi2*Q_t_pb(nSuppX*(T-2)+1:nSuppX*(T-1),:)/(capPi2)*(phi_wide(:,T-1)-phi_wide(:,T));
    C = phi_wide(:,T-2) - phi_wide(:,T-1);
    zerocond =  betadelta.*(A-B) + delta.*B -C;
    u1 = 1/(betadelta)*inv(capPi2)*(phi_wide(:,T-1)-phi_wide(:,T))-m1_wide(:,T);
    z = sum(zerocond.*zerocond);

    % Compute R-squared
    fitted = [A-B B]*([betadelta delta]');
    SSres = sum((C - fitted).*(C - fitted));
    SStot = sum((C - mean(C)).*(C - mean(C)));
    Rsqr = 1 - SSres/SStot;
end
