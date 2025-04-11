%{
    Function to calculate and verify moment conditions for naive agent.

%}
function [z,Rsqr] = verification_naive(beta,delta,supportX,capPi2,gamma,theta,CCP,CCP_tilde,V,flowpayoffs)
    nSuppX = length(supportX);
    % Check whether the moment conditions can solve the parameters
    % Since there is no inversion of any matrix, we can directly check the moment
    % condition from the proposition. Using the last three periods of data is enough.
    phi_21 = log(CCP(1:end,2)) - log(CCP(1:end,1));
    % Differencing the T-1 and T-2 phi_21 (27)
    dphi_21 = phi_21(end-nSuppX*3+1:end-nSuppX*2) - phi_21(end-nSuppX*2+1:end-nSuppX);
    zerocheckMoment = dphi_21 - beta*delta.*capPi2*(log(exp((1/beta.* phi_21(end-nSuppX*2+1:end-nSuppX)+(1-1/beta).*phi_21(end-nSuppX+1:end)))+ones(nSuppX,1)) +log(CCP(end-nSuppX+1:end,1)));
    z = sum(zerocheckMoment.*zerocheckMoment);
    % Compute R-squared.
    fitted = beta*delta.*capPi2*(log(exp((1/beta.* phi_21(end-nSuppX*2+1:end-nSuppX)+(1-1/beta).*phi_21(end-nSuppX+1:end)))+ones(nSuppX,1)) +log(CCP(end-nSuppX+1:end,1)));
    y = dphi_21;
    SSres = sum((y - fitted).*(y - fitted));
    SStot = sum((y - mean(y)).*(y - mean(y)));
    Rsqr = 1 - SSres/SStot;
end