% criterion function
% minimum distance between u1 solved from data in period 1-3 and u1
% solved from data in period 2-4
function [c,ceq] = myconstraint(beta,delta,theta,supportX,capPi2,CCP_total_sophi,nSuppX,nPeriods,gamma,flowpayoffs,CCP_generator,backInduction) 

   c = [];
   
   % calculate CCP based on the parameters
   [u1,u2] = flowpayoffs(supportX,theta); 
            
    % terminal period choice-specific value function:
    W_terminal = [u1,u2];
    % generating CCP by extreme value distribution
    CCP_terminal = CCP_generator(W_terminal);  
    % terminal period value function (not choice specific, i.e. ex ante). 
    V_terminal = W_terminal(:,end) + gamma-log(CCP_terminal(:,end));
    
    % backward induction to derive CCP for each period
    % for exponential case, beta = 1
    [CCP_total_sophi,~,~] = backInduction(CCP_terminal,V_terminal,u1,u2,gamma,capPi2,beta,beta,delta,nPeriods);
       
   
    dphi_21 = log(CCP_total_sophi(nSuppX+1:end,2)) - log(CCP_total_sophi(nSuppX+1:end,1)) - ...
    ( log(CCP_total_sophi(1:nSuppX*(nPeriods-1),2)) - log(CCP_total_sophi(1:nSuppX*(nPeriods-1),1)) );
    m1 = gamma - log(CCP_total_sophi(:,1));
    t = 1;
    Q_2_pb = repmat(CCP_total_sophi( nSuppX*t+1 : nSuppX*t+nSuppX,2),1,nSuppX).*capPi2;
    Q_3_pb = repmat(CCP_total_sophi( nSuppX*(t+1)+1 : nSuppX*(t+1)+nSuppX,2),1,nSuppX).*capPi2;
    Q_4_pb = repmat(CCP_total_sophi( nSuppX*(t+2)+1 : nSuppX*(t+2)+nSuppX,2),1,nSuppX).*capPi2;
        
     Left2 = delta*(1-beta).*Q_2_pb;
     Left3 = delta*(1-beta).*Q_3_pb;
     Left4 = delta*(1-beta).*Q_4_pb;
     Right2 = 1/(beta*delta).*inv(capPi2)*dphi_21(nSuppX*1-nSuppX+1:nSuppX*1)+m1(nSuppX*2-nSuppX+1:nSuppX*2);
     Right3 = 1/(beta*delta).*inv(capPi2)*dphi_21(nSuppX*2-nSuppX+1:nSuppX*2)+m1(nSuppX*3-nSuppX+1:nSuppX*3);
     Right4 = 1/(beta*delta).*inv(capPi2)*dphi_21(nSuppX*3-nSuppX+1:nSuppX*3)+m1(nSuppX*4-nSuppX+1:nSuppX*4);
    
     % express u1 as function of beta and delta
     LeftMatrix_123 = inv(eye(nSuppX)-Left2) - eye(nSuppX)- Left3/(eye(nSuppX)-Left3);
     LeftMatrix_234 = inv(eye(nSuppX)-Left3) - eye(nSuppX)- Left4/(eye(nSuppX)-Left4);
     RightVector_123 = - (eye(nSuppX)-Left2)\Right2 + m1(nSuppX*3-nSuppX+1:nSuppX*3) + Left3*((eye(nSuppX)-Left3)\Right3);
     RightVector_234 = - (eye(nSuppX)-Left3)\Right3 + m1(nSuppX*4-nSuppX+1:nSuppX*4) + Left4*((eye(nSuppX)-Left4)\Right4);
     u1_123= LeftMatrix_123\RightVector_123;
     u1_234= LeftMatrix_234\RightVector_234;
     ceq = sum((u1_123-u1_234).*(u1_123-u1_234));
end