% backward induction to derive the CCP for each period given the terminal
% period
function [CCP_total,V_total,W_total,CCP_total_tilde] = backInduction(CCP_T,V_T,u1,u2,gamma,capPi2,beta,betatilde,delta,nPeriods)
    [nSuppX,nA]=size(CCP_T);
    V_last = V_T;
    CCP_total = CCP_T;
    CCP_total_tilde = CCP_T;
    V_total = V_last;
    
    % Note that the choice specific value function for the final period is
    % just the flow utility. Here, the flow utility is invariant for different
    % previous actions.
    W_total = [u1,u2];
    
    for t = 1:nPeriods-1
        if (betatilde ==beta) % combine consistent agent and sophisticated agent
            % sophisticated case: beta != 1 and betatilde = beta
            [W] = bellman(u1,u2,beta,delta,V_last,capPi2);
            % then calculate CCP for t-1 period
            CCP = CCP_generator(W);
            % calculate the ex ante value function (not choice specific)
            V_last = W(:,1)+gamma-log(CCP(:,1)) + ...
                delta*(1-beta)*CCP(:,2).*(capPi2*V_total(1:nSuppX));
        else            
            % naive case: beta!=1 and betatilde != beta
            [W] = bellman(u1,u2,beta,delta,V_last,capPi2);
            [Z] = bellman(u1,u2,betatilde,delta,V_last,capPi2);
            % then calculate CCP for t-1 period
            CCP = CCP_generator(W);
            CCP_tilde = CCP_generator(Z);
            % calculate the value function (not choice specific)
            V_last = Z(:,1)+gamma-log(CCP_tilde(:,1)) + ...
                delta*(1-betatilde)*CCP_tilde(:,2).*(capPi2*V_total(1:nSuppX));
            CCP_total_tilde = [CCP_tilde;CCP_total_tilde];
        end
        CCP_total = [CCP;CCP_total];
        V_total = [V_last;V_total];
        W_total = [W;W_total];
    end
    V_total = V_total(nSuppX+1:nSuppX*nPeriods,1);
end






