% logit-fit of CCPs
function CCP_adjusted = logitfit(CCP_original,nFirms)
     nRow = size(CCP_original,1);
     nA = size(CCP_original,2);
     CCP_adjusted = zeros(nRow,nA);
     for i = 1:nRow
         CCP_adjusted(i,:) = CCP_original(i,:);
         for j = 1: nA
             if CCP_original(i,j) ==0
                 CCP_adjusted(i,:) = CCP_original(i,:)-1/nFirms/(nA-1);
                 CCP_adjusted(i,j) = CCP_adjusted(i,j) + 1/nFirms *nA/(nA-1);
             end
         end
     end
end