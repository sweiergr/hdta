function z = verification_naive(beta,delta,supportX,capPi2,gamma,theta,CCP,CCP_tilde,V,flowpayoffs)

nSuppX = length(supportX);
[u1,u2] = flowpayoffs(supportX,theta);
% verification part
%% check whether the moment conditions can recover the parameters. Since there is no inversion to any matrix, directly check the moment condition.
% Using the last three periods of data should be enough, see proposition in paper.
phi_21 = log(CCP(1:end,2)) - log(CCP(1:end,1));
% differencing the T-1 and T-2 phi_21 (Eq 27 in paper)
dphi_21 = phi_21(end-nSuppX*3+1:end-nSuppX*2) - phi_21(end-nSuppX*2+1:end-nSuppX);
% % check (32 in paper)
zerocheckMoment = dphi_21 - beta*delta.*capPi2*(log(exp((1/beta.* phi_21(end-nSuppX*2+1:end-nSuppX)+(1-1/beta).*phi_21(end-nSuppX+1:end)))+ones(nSuppX,1)) +log(CCP(end-nSuppX+1:end,1)));
z = sum(zerocheckMoment.*zerocheckMoment);
end