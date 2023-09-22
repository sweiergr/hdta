% criterion function:
% minimum distance between u1 solved from data in period 1-3 and u1
% solved from data in period 2-4
function [z] = criterianNoinv(beta,delta,theta,capPi2,CCP_total_sophi,supportX,nPeriods,gamma,flowpayoffs) 

   nSuppX = length(supportX);
   [u1,~] = flowpayoffs(supportX,theta);
    dphi_21 = log(CCP_total_sophi(nSuppX+1:end,2)) - log(CCP_total_sophi(nSuppX+1:end,1)) - ...
    ( log(CCP_total_sophi(1:nSuppX*(nPeriods-1),2)) - log(CCP_total_sophi(1:nSuppX*(nPeriods-1),1)) );
    m1 = gamma - log(CCP_total_sophi(:,1));
    t = 1;
    Q_2_pb = repmat(CCP_total_sophi( nSuppX*t+1 : nSuppX*t+nSuppX,2),1,nSuppX).*capPi2;
    Q_3_pb = repmat(CCP_total_sophi( nSuppX*(t+1)+1 : nSuppX*(t+1)+nSuppX,2),1,nSuppX).*capPi2;
    Q_4_pb = repmat(CCP_total_sophi( nSuppX*(t+2)+1 : nSuppX*(t+2)+nSuppX,2),1,nSuppX).*capPi2;
    h_2 = inv(eye(nSuppX) - delta*(1-beta).*Q_2_pb)*(1/(beta*delta).*inv(capPi2)*dphi_21(nSuppX*(t-1)+1:nSuppX*(t-1)+nSuppX)+m1(nSuppX*t+1:nSuppX*t+nSuppX) + u1);
    h_3 = inv(eye(nSuppX) - delta*(1-beta).*Q_3_pb)*(1/(beta*delta).*inv(capPi2)*dphi_21(nSuppX*(t+0)+1:nSuppX*(t+0)+nSuppX)+m1(nSuppX*(t+1)+1:nSuppX*(t+1)+nSuppX) + u1);
    h_4 = inv(eye(nSuppX) - delta*(1-beta).*Q_4_pb)*(1/(beta*delta).*inv(capPi2)*dphi_21(nSuppX*(t+1)+1:nSuppX*(t+1)+nSuppX)+m1(nSuppX*(t+2)+1:nSuppX*(t+2)+nSuppX) + u1);
    tempz_123= h_2 - m1(nSuppX*(t+1)+1:nSuppX*(t+1)+nSuppX) - u1 - delta*(1-beta).*Q_3_pb*h_3;
    tempz_234 = h_3 - m1(nSuppX*(t+2)+1:nSuppX*(t+2)+nSuppX) - u1 - delta*(1-beta).*Q_4_pb*h_4;
    z = sum(tempz_123.*tempz_123 + tempz_234.*tempz_234);
end