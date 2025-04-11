%{
    Function to implement OLS-estimator using the final three periods of data.
    Inputs:
    - nPeriods: number of time periods in data.
    - supportX: support of state vector.
    - capPi2: transition matrix for state variable
    - gamma: mean of the random shock
    - CCP: CCPs as estimated from the data.
    Output is estimates of discount factors and the reference utility.
%}

function [EstResult, SSres]=EstOLS_LastPeriods(nPeriods,supportX,gamma,capPi2,CCP)
    nSuppX = length(supportX);
    % Action 2 is the continuation action (k), and action 1 is the terminating
    % action(K)
    Q_t_bar = repmat(CCP(:,2),1,nSuppX).*repmat(capPi2,nPeriods,1);
    phi_21 = log(CCP(:,2)) - log(CCP(:,1));
    phi_wide = reshape(phi_21,[nSuppX,nPeriods]);
    CCP_K_wide = reshape(CCP(:,1),[nSuppX,nPeriods]);
    T = nPeriods;
    Q_T_1_bar = Q_t_bar(nSuppX*(T-2)+1:nSuppX*(T-1),:);
    MatrixA = capPi2*(log(CCP_K_wide(:,T)) - log(CCP_K_wide(:,T-1)) - ...
                     Q_T_1_bar /(capPi2) * (phi_wide(:,T-1) - phi_wide(:,T)));
    MatrixB = capPi2* Q_T_1_bar /(capPi2) * (phi_wide(:,T-1) - phi_wide(:,T));
    Omega = [MatrixA MatrixB];
    z = (Omega'*Omega)\Omega'*(phi_wide(:,T-2) - phi_wide(:,T-1));
    
    % Compute R-squared
    fitted = Omega*z;
    y = phi_wide(:,T-2) - phi_wide(:,T-1);
    SSres = sum((y - fitted).*(y - fitted));
    SStot = sum((y - mean(y)).*(y - mean(y)));
    Rsqr = 1 - SSres/SStot;
    
    betadelta = z(1);
    delta = z(2);
    beta=z(1)/z(2);
    discount=[beta delta];
    u1 = -gamma + log(CCP_K_wide(:,T)) + 1/betadelta.*inv(capPi2)*(phi_wide(:,T-1) - phi_wide(:,T));
    EstResult = [discount'; u1];
end