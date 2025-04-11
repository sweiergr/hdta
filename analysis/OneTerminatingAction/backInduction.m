%{
    Function to derive CCPs for each period using backward induction given the terminal
    period. Inputs to the function are:
    - CCP_T: CCPs for the last period T
    - V_T: ex-ante value function for the last period T
    - u1: flow utility for taking action 1
    - u2: flow utility for taking action 0
    - gamma: mean of the random shock
    - capPi2: transition matrix for state variable
    - beta: present-bias discount factor
    - betatilde: present-bias discount factor for naive agent
    - delta: exponential discount factor
    - nPeriods: number of time periods in the data.

    Outputs:
    - CCP_total: CCPs for all types of agents based on the parameter beta and delta
    - V_total: perceived long-run value function.
    - W_total: choice-specific value function.
    - CCP_total_tilde: perceived CCPs for naive agents based on parameter betatilde and delta



%}
function [CCP_total,V_total,W_total,CCP_total_tilde] = backInduction(CCP_T,V_T,u1,u2,gamma,capPi2,beta,betatilde,delta,nPeriods)
    
    [nSuppX,nA]=size(CCP_T);
    V_last = V_T;
    CCP_total = CCP_T;
    CCP_total_tilde = CCP_T;
    V_total = V_last;
    
    % Note that the choice specific value function for the final period is
    % just flow utility, and here flow utility is invariant for different
    % previous action.
    W_total = [u1,u2];
    
    % Loop over time periods.
    for t = 1:nPeriods-1
    % First, compute the choice specific value function for period t-1.
    % Note that V_last is the ex-ante value function, not choice-specific.
    % Each column corresponds to a specific choice in the last period.
        if (betatilde ==beta) % Combine consistent agent and sophisticated agent
            % Sophisticated case: beta != 1 and betatilde = beta.
            [W] = bellman(u1,u2,beta,delta,V_last,capPi2);
            % Calculate CCP for t-1 period
            CCP = CCP_generator(W);
            % Calculate the ex ante value function.
            V_last = W(:,1)+gamma-log(CCP(:,1)) + ...
                delta*(1-beta)*CCP(:,2).*(capPi2*V_total(1:nSuppX));
        else            
            % Naive agent case: beta!=1 and betatilde != beta
            [W] = bellman(u1,u2,beta,delta,V_last,capPi2);
            [Z] = bellman(u1,u2,betatilde,delta,V_last,capPi2);
            % Calculate CCP for period t-1.
            CCP = CCP_generator(W);
            CCP_tilde = CCP_generator(Z);
            % Calculate the ex-ante value function.
            V_last = Z(:,1)+gamma-log(CCP_tilde(:,1)) + ...
                delta*(1-betatilde)*CCP_tilde(:,2).*(capPi2*V_total(1:nSuppX));
            CCP_total_tilde = [CCP_tilde;CCP_total_tilde];
        end % end distinction between sophisticated and naive agents case.

        CCP_total = [CCP;CCP_total];
        V_total = [V_last;V_total];
        W_total = [W;W_total];
    end % end loop over time periods.
    % Drop the zero period long-run continuation value part since we start
    % from the first period. This long-run continuation value is not
    % meaningful for the zero period.
    V_total = V_total(nSuppX+1:nSuppX*nPeriods,1);
end