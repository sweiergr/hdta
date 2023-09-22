function [W] = bellman(u1,u2,beta,delta,V,capPi2)
    % This is the choice-specific value function
    % continuation value
    r2 = (beta*delta).* (capPi2*V);
    % value function
    W = [u1, u2+r2];
end 