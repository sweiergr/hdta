%{
    Function to calculate objective function for minimum distance estimation using final periods of data and additional exclusion restriction.
    Inputs:
    - nPeriods: number of time periods in the data.
    - supportX: support of state vector.
    - capPi2: transition matrix for state variable
    - beta: present-bias discount factor
    - delta: exponential discount factor
    - gamma: mean of the random shock
    - CCP: CCPs as estimated from data

%}
function [z, uK, SSres]=EstMD_LastPeriods_restriction(nPeriods,supportX,capPi2,beta,delta,gamma,CCP)
    
    nSuppX = length(supportX);
    Q_t_pb = repmat(CCP(:,2),1,nSuppX).*repmat(capPi2,nPeriods,1);
    betadelta = beta*delta;
    phi_21 = log(CCP(:,2)) - log(CCP(:,1));
    m1 = gamma - log(CCP(:,1));
    phi_wide = reshape(phi_21,[nSuppX,nPeriods]);
    m1_wide = reshape(m1,[nSuppX,nPeriods]);
    logCCP_wide=reshape(log(CCP(:,1)),[nSuppX,nPeriods]);
    T = nPeriods;
    % Estimation incorporating the exclusion restriction.
    b=[1 0;1 0;0 1;0 1];
    a=capPi2*b;
    coef=(a'*a)\a';
    uK_tilde=coef*(capPi2*logCCP_wide(:,T)-1/betadelta*(phi_wide(:,T) - phi_wide(:,T-1)));
    uK=b*uK_tilde;
    A = capPi2*(m1_wide(:,T-1) - m1_wide(:,T));
    C = phi_wide(:,T-2) - phi_wide(:,T-1);
    B = capPi2*Q_t_pb(nSuppX*(T-2)+1:nSuppX*(T-1),:)*( m1_wide(:,T)+uK);
    zerocond =  betadelta*A+ betadelta*(delta-betadelta)*B-C;
    z = sum(zerocond.*zerocond);
    
    % Compute R-squared
    fitted = betadelta*A+ betadelta*(delta-betadelta)*B;
    SSres = sum((C - fitted).*(C - fitted));
    SStot = sum((C - mean(C)).*(C - mean(C)));
    Rsqr = 1 - SSres/SStot;
end
