% This file implements our OLS estimation, returning discount factors and
% the referene utility.
% Using the last three periods data.

function [EstResult]=EstOLS_LastPeriods(nPeriods,supportX,gamma,capPi2,CCP)

    nSuppX = length(supportX);
    % action 2 is the continuation action(k), and action 1 is the terminating
    % action(K)
    Q_t_bar = repmat(CCP(:,2),1,nSuppX).*repmat(capPi2,nPeriods,1);
    
    phi_21 = log(CCP(:,2)) - log(CCP(:,1));
    phi_wide = reshape(phi_21,[nSuppX,nPeriods]);
    CCP_K_wide = reshape(CCP(:,1),[nSuppX,nPeriods]);
    
    T = nPeriods;
    Q_T_1_bar = Q_t_bar(nSuppX*(T-2)+1:nSuppX*(T-1),:);
    MatrixA = capPi2*(log(CCP_K_wide(:,T)) - log(CCP_K_wide(:,T-1)) - ...
                     Q_T_1_bar * inv(capPi2) * (phi_wide(:,T-1) - phi_wide(:,T)));
    MatrixB = capPi2* Q_T_1_bar * inv(capPi2) * (phi_wide(:,T-1) - phi_wide(:,T));
    Omega = [MatrixA MatrixB];
    z = (Omega'*Omega)\Omega'*(phi_wide(:,T-2) - phi_wide(:,T-1));
    betadelta = z(1);
    delta = z(2);
    beta=z(1)/z(2);
    discount=[beta delta];
    u1 = -gamma + log(CCP_K_wide(:,T)) + 1/betadelta.*inv(capPi2)*(phi_wide(:,T-1) - phi_wide(:,T));
    EstResult = [discount'; u1];
end