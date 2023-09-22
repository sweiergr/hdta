function [u1,u2] = flowpayoffs(supportX,theta)
    u1 = theta(1) + theta(2) *supportX;
    u2 = theta(3) + theta(4)*supportX;
end
