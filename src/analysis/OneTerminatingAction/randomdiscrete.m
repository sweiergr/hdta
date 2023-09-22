function y = randomdiscrete(p)
    nSupp = size(p,1);
    nVar = size(p,2);
    % iid sampling from uniform distribution [0,1]
    uniformDraws = ones(nSupp-1,1)*random('unif',zeros(1,nVar),ones(1,nVar));
    % get the cumulative prob of p for each column (i.e. for each firm)
    cumulativeP = cumsum(p);
    y = sum([ones(1,nVar);cumulativeP(1:nSupp-1,:) <= uniformDraws]);
end