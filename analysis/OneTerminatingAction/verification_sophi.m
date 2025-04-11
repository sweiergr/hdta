%{
    Function to verify moment conditions for sophisticated agent.

%}
function [z, u1_t, Matrix2_t,Matrix5_t,dphi_21,m1]=verification_sophi(supportX,nPeriods,capPi2,beta,delta,gamma,theta,CCP,V,verification,flowpayoffs)

nSuppX = length(supportX);
[u1,~] = flowpayoffs(supportX,theta);
I = eye(nSuppX);
residpb = delta*(1-beta);

% Part 1: Check the full rank condition required by Theorem.
% 1. Check whether the composite transition matrix is singular or not by
% checking the singular value.
    % Calculate Q_t_bar.
    Q_t_bar = repmat(CCP(:,2),1,nSuppX).*repmat(capPi2,nPeriods,1);
    Matrix2_t = capPi2;
    Matrix5_t = [];
    LeftMatrix_t = [];
    singularvalue_5 = zeros(nSuppX,nPeriods-1); 
    condnum2 = [];
    condnum5 = [];
   for t = 1:nPeriods-1
    Matrix5 = Q_t_bar(nSuppX*(t-1)+1:nSuppX*t,:)-Q_t_bar(nSuppX*(t)+1:nSuppX*(t+1),:);
    Matrix5_t = [Matrix5_t;Matrix5];
    LeftMatrix = residpb*inv(I-residpb.*Q_t_bar(nSuppX*(t-1)+1:nSuppX*t,:)) *Matrix5 *inv(I-residpb.*Q_t_bar(nSuppX*(t)+1:nSuppX*(t+1),:));
    LeftMatrix_t = [LeftMatrix_t;LeftMatrix];
        if verification ==1
            singularvalue_2 = svd(Matrix2_t); % calculate singular value
            condnum2 = [cond(Matrix2_t)]; % calculate condition number
            singularvalue_5(:,t) = svd(Matrix5); % calculate singular value
            condnum5 = [condnum5;cond(Matrix5)]; % calculate condition number
        end
   end
   
   % Part 2: Check whether the moment condition solves the parameters.
   % Moment condition satisfies rank condition.
   % Step 1: Given beta,delta and u1, check whether the equations hold.
    dphi_21 = log(CCP(nSuppX+1:end,2)) - log(CCP(nSuppX+1:end,1)) - ...
        ( log(CCP(1:nSuppX*(nPeriods-1),2)) - log(CCP(1:nSuppX*(nPeriods-1),1)) );
    m1 = gamma - log(CCP(:,1));
    dphi_wide = reshape(dphi_21,[nSuppX,nPeriods-1]);
    m1_wide = reshape(m1,[nSuppX,nPeriods]);
    V_wide = reshape(V,[nSuppX,nPeriods-1]);
    CCP_K_wide = reshape(CCP(:,1),[nSuppX,nPeriods]);
    CCP_k_wide = reshape(CCP(:,2),[nSuppX,nPeriods]);
   % Calculate the matrix Omega.
   T = nPeriods;
   Q_T_1_bar = Q_t_bar(nSuppX*(T-2)+1:nSuppX*(T-1),:);
   MatrixA = capPi2*(log(CCP_K_wide(:,T-1)) - log(CCP_K_wide(:,T)) - ...
                     Q_T_1_bar * inv(capPi2) * (dphi_wide(:,end)));
   MatrixB = capPi2* Q_T_1_bar * inv(capPi2) * (dphi_wide(:,end));
   Omega = [MatrixA MatrixB];
   singular_Omega = svd(Omega);
   svd(Omega'*Omega);
   % Verify discount factors based on Omega.
   DiscountFactor = inv(Omega'*Omega)*Omega'*(dphi_wide(:,end-1));
    
   % Part 3: Verifying non-stationarity of CCPs.
   % Calculate the differences of consecutive periods CCP_k_wide.
   dCCP_k_wide = CCP_k_wide(:,1:end-1) - CCP_k_wide(:,2:end);

if verification ==1
    disp(["condnum2",num2str(condnum2')]);
    disp(["condnum5",num2str(condnum5')]);
    minsingular2 = min(min(singularvalue_2));
    minsingular5 = min(min(singularvalue_5));
    disp(["min singular:",num2str(min([minsingular2,minsingular5]))]);
end
% Export the singularity tables and condition numbers to LaTeX.
if verification == 1
    % First, put the singular values and conditional numbers into the same
    % table
    singular2.data = [singularvalue_2'];
    for i = 1:nSuppX
        singular2.tableColLabels{i} = ['SV',num2str(i)];
    end
    % singular2.tableColLabels{nSuppX+1} = 'Condition number';
    singular2.dataFormatMode = 'col';
    singular2.tablePlacement = 'h!';
    singular2.dataFormat = {'%.4f'};
    singular2.tableBorders = 3; % three line table; 1 is for all border table
    singular2.tableColumnAlignment = 'c';
    %     singular2.tableCaption = 'Singular values and condition numbers for $\bm{Q}_2(I - \delta (1-\beta)\bm{Q}^{pb}_{t+1})$ ';
    singular2.tableCaption = 'Singular values for $\bm{Q}_0$ ';
    singular2.tableLabel = 'singular value Q_pb';
    singular2.makeCompleteLatexDocument = 0;
    latex = latexTable(singular2);
    FID = fopen(project_paths('OUT_TABLES','singularityQ.tex'), 'w');
    for k=1:size(latex, 1)
        fprintf(FID, '%s \n', latex{k,:});
    end
    fclose(FID);
    clear FID

    % Export Omega and singular_Omega
    singularOmega.data = [Omega;singular_Omega'];
    for i = 1:2
        singularOmega.tableColLabels{i} = ['col',num2str(i)];
    end
    for i = 1:nSuppX
        singularOmega.tableRowLabels{i} = ['$\Omega_{',num2str(i),',\cdot}$'];
    end
    singularOmega.tableRowLabels{nSuppX+1} = 'SV($\Omega$)';
    singularOmega.dataFormatMode = 'row';
    singularOmega.tablePlacement = 'h!';
    singularOmega.dataFormat = {'%.4f'};
    singularOmega.tableBorders = 3; % three line table; 1 is for all border table
    singularOmega.tableColumnAlignment = 'c';
    %     singularOmega.tableCaption = 'Singular values and condition numbers for $\bm{Q}_2(I - \delta (1-\beta)\bm{Q}^{pb}_{t+1})$ ';
    singularOmega.tableCaption = '$\Omega$ and its singular values ';
    singularOmega.tableLabel = 'Omega and its singular values';
    singularOmega.makeCompleteLatexDocument = 0;
    latex = latexTable(singularOmega);
    FID = fopen(project_paths('OUT_TABLES','singularityOmega.tex'), 'w');
    for k=1:size(latex, 1)
        fprintf(FID, '%s \n', latex{k,:});
    end
    fclose(FID);
    clear FID
    
    % Export DiscountFactor verified.
    DF.data = [DiscountFactor'];
    DF.tableColLabels = {'$\delta\beta$','$\delta$'};
    DF.dataFormatMode = 'col';
    DF.tablePlacement = 'h!';
    DF.dataFormat = {'%.4f'};
    DF.tableBorders = 3; % three line table; 1 is for all border table
    DF.tableColumnAlignment = 'c';
    %       DF.tableCaption = 'Singular values and condition numbers for $\bm{Q}_2(I - \delta (1-\beta)\bm{Q}^{pb}_{t+1})$ ';
    DF.tableCaption = 'Verified discount factors ';
    DF.tableLabel = 'verified DF';
    DF.makeCompleteLatexDocument = 0;
    latex = latexTable(DF);
    FID = fopen(project_paths('OUT_TABLES','verifiedDF.tex'), 'w');
    for k=1:size(latex, 1)
        fprintf(FID, '%s \n', latex{k,:});
    end
    fclose(FID);
    clear FID
    
    % Export dCCP_k_wide
    CCPk.data = dCCP_k_wide;
    for t = 1:nPeriods-2
        CCPk.tableColLabels{t} = ['$p_{0,T-',num2str(nPeriods-t),'}-p_{0,T-',num2str(nPeriods-t-1),'}$'];
    end
    CCPk.tableColLabels{nPeriods-1} = ['$p_{0,T-1}-p_{0,T}$'];
    for i = 1:nSuppX
        CCPk.tableRowLabels{i} = [char(supportX(i))];
    end
    CCPk.dataFormatMode = 'row';
    CCPk.tablePlacement = 'h!';
    CCPk.dataFormat = {'%.4f'};
    CCPk.tableBorders = 3; % three line table; 1 is for all border table
    CCPk.tableColumnAlignment = 'c';
%     CCPk.tableCaption = 'Singular values and condition number for $ (I - \delta (1-\beta)\bm{Q}^{pb}_{t})^{-1} -I -\delta(1-\beta)\bm{Q}^{pb}_t(I - \delta (1-\beta)\bm{Q}^{pb}_{t+1})^{-1} $ ';
    CCPk.tableCaption = '$p_{k,t}-p_{k,t+1}$';
    CCPk.tableLabel = 'CCP difference';
    CCPk.makeCompleteLatexDocument = 0;
    latex = latexTable(CCPk);
    FID = fopen(project_paths('OUT_TABLES','CCPdifference.tex'), 'w');
    for k=1:size(latex, 1)
        fprintf(FID, '%s \n', latex{k,:});
    end
    fclose(FID);
    clear FID
end

Q_beta = delta*(1-beta).*Q_t_bar; 
Right =  1/(beta*delta).*inv(capPi2)*dphi_wide+m1_wide(:,2:end);
RightVector = [];
u1_t = [];
zerocheck1_t = [];
zerocheck2_t = [];
zerocheck3_t = [];
zerocheck4_t = [];
zerocheck_u_t = [];

for t = 1:nPeriods-2
    Q_t = Q_beta(nSuppX*t+1:nSuppX*(t+1),:);
    Q_tplus1 = Q_beta(nSuppX*(t+1)+1:nSuppX*(t+2),:);
    RV =  -(I-Q_t)\Right(:,t) + m1_wide(:,t+2)+Q_tplus1*((I - Q_tplus1)\Right(:,t+1));
    RightVector = [RightVector,RV];
    % u1_t are the u1 calculated with different three periods of data,
    % For example, the first column is from the period 1-3, and the second
    % column is from the period 2-4.
    u1_t = [u1_t,LeftMatrix_t(nSuppX*t+1 : nSuppX*t+nSuppX,:)\RightVector(:,t)];
    % Zero check 
    if verification == 1        
        % check equation (18) of V1 expression
        zerocheck1 = dphi_wide(:,t) - beta*delta.*capPi2*(V_wide(:,t+1)-V_wide(:,t));
        zerocheck1_t = [zerocheck1_t,zerocheck1];
        Q_tplus1_pb = Q_t_bar(nSuppX*t+1:nSuppX*t+nSuppX,:);
        zerocheck2 = dphi_wide(:,t) - beta*delta.*capPi2*(I - delta*(1-beta).*Q_tplus1_pb)*V_wide(:,t+1)...
                                            + beta*delta.*capPi2*(m1_wide(:,t+2)+u1);
        zerocheck2_t = [zerocheck2_t,zerocheck2];                          

        % check (19) (v expression)
        zerocheck3 = V_wide(:,t) - m1_wide(:,t+1)-u1-delta*(1-beta).*Q_tplus1_pb*V_wide(:,t+1);
        zerocheck3_t = [zerocheck3_t,zerocheck3];

        % check (18) last equality 
        zerocheck4 = V_wide(:,t+1) - inv(I - delta*(1-beta).*Q_tplus1_pb)*(1/(beta*delta).*inv(capPi2)*dphi_wide(:,t)+m1_wide(:,t+1) + u1);
        zerocheck4_t = [zerocheck4_t,zerocheck4];

        zerocheck_u= LeftMatrix_t(nSuppX*t+1 : nSuppX*t+nSuppX,:)*u1 ...
     -(- (I-Q_t)\Right(:,t) + m1_wide(:,t+2)+Q_tplus1 *((I - Q_tplus1)\Right(:,t+1)));
        zerocheck_u_t = [zerocheck_u_t,zerocheck_u];
    end
end
    if verification == 1
        % export the true u1, and the computed utility u1_t
        u_verify.data = [u1,u1_t];
        for x = 1:nSuppX
            % u_verify.tableRowLabels{x} = ['x=',num2str(x)];
            u_verify.tableRowLabels{x} = [char(supportX(x))];
        end
		u_verify.tableColLabels{1} = 'True $u_1$';
        u_verify.tableColLabels{2} = '$u^{123}_1$';
        u_verify.tableColLabels{3} = '$u^{234}_1$';
        u_verify.dataFormatMode = 'row';
        u_verify.tablePlacement = 'h!';
        u_verify.dataFormat = {'%.4f'};
        u_verify.tableBorders = 3; % three line table; 1 is for all border table
        u_verify.tableColumnAlignment = 'c';
        u_verify.tableCaption = 'true $u_1$ and computed $u_1$ from three consecutive periods of data';
        u_verify.tableLabel = 'u_verify';
        u_verify.makeCompleteLatexDocument = 0;
        latex = latexTable(u_verify);
        FID = fopen(project_paths('OUT_TABLES','Uverify.tex'), 'w');
        for k=1:size(latex, 1)
            fprintf(FID, '%s \n', latex{k,:});
        end
        fclose(FID);
        clear FID
        
        zerocheck = max(max([zerocheck1_t,zerocheck2_t,zerocheck4_t]));
        disp(["verify: eqn(18)",num2str(zerocheck)]);
        disp(["verify: V",num2str(max(max(zerocheck3_t)))]);
        disp(["verify u express: ",num2str(max(max(zerocheck_u_t)))]);
    end
    t = nPeriods -3;
    z = sum((u1_t(:,t)-u1_t(:,t+1)).*(u1_t(:,t)-u1_t(:,t+1)));
end