%{
    Function to calculate CCPs for a standard logit model as function of choice-specific value function.
%}

function [CCP]=CCP_generator(W)
    nSuppX = size(W,1);
    nA = size(W,2);
    
    CCP = zeros(nSuppX,nA);
    for xi = 1:nSuppX
        for ai = 1:nA
            max_W=max(W(xi,:));
            CCP(xi,ai) = exp(W(xi,ai)-max_W)/sum(exp(W(xi,:)-max_W));
        end
    end
end