% Generate CCP as function of deltaV 
% Written by Chao Wang 
% Created on Dec 26th, 2019

function [CCP]=CCP_generator(W)
    nSuppX = size(W,1);
    nA = size(W,2);
    CCP = zeros(nSuppX,nA);
    for xi = 1:nSuppX
        for ai = 1:nA
            CCP(xi,ai) = exp(W(xi,ai))/sum(exp(W(xi,:)));
        end
    end
end