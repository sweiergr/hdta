%{
    Function to compute flow utility from actions 1 and 2 as function of state variables *supportX* and parameter vector *theta.

%}
function [u1,u2] = flowpayoffs(supportX,theta)
% Utility function used for DDCOneTerminate.m
       nSuppX = length(supportX);
       u = reshape(theta,nSuppX,2);
       u1 = u(:,1);
       u2 = u(:,2);
end
