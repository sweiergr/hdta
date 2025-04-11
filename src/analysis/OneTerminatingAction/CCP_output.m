%{
    Function to calculate CCPs given the estimated parameters.
    Inputs:
    - theta: parameters in the flow utility
    - beta: present-bias discount factor 
    - betatilde: present-bias discount factor for naive agent
    - delta: exponential discount factor
    - nPeriods: number of time periods in the data.
    - gamma: mean of the random shock
    - supportX: support of state vector.
    - capPi2: transition matrix for state variable
    - backInduction: backInduction function
    - CCP_generator: CCP_generator function 
    - flowpayoffs: flowpayoffs function
    
    Output: CCPs computed by backward induction.

%}



function [CCPs] = CCP_output(theta,beta,betatilde,delta,nPeriods,gamma,supportX,capPi2,...
                    backInduction,CCP_generator,flowpayoffs)
    
   [u1,u2] = flowpayoffs(supportX,theta);
    W_terminal = [u1,u2];
    CCP_terminal = CCP_generator(W_terminal);    
    V_terminal = W_terminal(:,1) + gamma-log(CCP_terminal(:,1)); % action 1 is the terminating action      
    [CCPs,~,~,~] = backInduction(CCP_terminal,V_terminal,... 
                             u1,u2,gamma,capPi2,beta,betatilde,delta,nPeriods);
end