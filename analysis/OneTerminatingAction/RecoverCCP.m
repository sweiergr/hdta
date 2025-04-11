%{
    Function to recover CCps from  observed choices and index matrices.
    Inputs:
    - choices: Observed choices in the data.
    - iX: Index of the states
    - nPeriods: Number of periods of data.
    - nSuppX: Number of states
    - nA: Number of actions
%}
function CCP_recover = RecoverCCP(choices,iX,nPeriods,nSuppX,nA)
    nFirms = size(choices,2);
    CCP_recover = zeros(nSuppX*nPeriods,nA);
    CCP_temp = zeros(nSuppX,nA);
    for ti = 1:nPeriods
        for xi = 1:nSuppX
            for k = 1:nA
                CCP_temp(xi,k) = sum(choices(ti,:)==k & iX(ti,:)==xi)/...
                    sum(iX(ti,:)==xi);
            end
        end
        CCP_recover(nSuppX*(ti-1)+1: nSuppX*(ti-1)+nSuppX,:)= CCP_temp;
    end
    CCP_recover = logitfit(CCP_recover,nFirms);
end