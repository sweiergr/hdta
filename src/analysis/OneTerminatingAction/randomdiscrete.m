%{		
	The function |randomDiscrete| returns a random draw from the distribution of $(Y_1,\ldots,Y_n)$; with $Y_1,\ldots,Y_n$ independently discretely distributed with (not necessarily identical) distributions on $\{1,2,\ldots,k\}$ for some $2\leq k<\infty$.

    This function follows "Dynamic Discrete Choice Models: Methods, Matlab Code,
    and Exercises" by Jaap H. Abbring and Tobias J. Klein.
    See http://ddc.abbring.org for the original code and extensive documentation.
%}

function y = randomdiscrete(p)
    nSupp = size(p,1);
    nVar = size(p,2);
    
    % iid sampling from uniform distribution [0,1]
    uniformDraws = ones(nSupp-1,1)*random('unif',zeros(1,nVar),ones(1,nVar));
    % get the cumulative prob of p for each column (i.e. for each firm)
    cumulativeP = cumsum(p);
    y = sum([ones(1,nVar);cumulativeP(1:nSupp-1,:) <= uniformDraws]);
end