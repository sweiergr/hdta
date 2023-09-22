
%{
    This is the main file to replicate the simulations in the paper.
    Written by Chao Wang, Stefan Weiergraeber, Ruli Xiao in June 2022.

%}

clear
clc

maxj = 100; % maximum times of sampling
epsilon = [0.01];

nPeriods = 4;
nFirms_alternatives =[10000,50000,100000,1000000,20000,5000];
nA = 2;

thetaY = 1; % this equals to one means theta is estimated
thetaNormalizationT = 2; % this equals to 2 means normalized theta (normalizing the terminating action) is estimated
thetaNormalization = 3; % this equals to 3 means normalized theta (normalizing the other action) is estimated

nSuppX =4;
supportX = [2,3,7,9]';
% since the transitional matrix is different with action 2 and action 3,
% set capPi2 and capPi3 separately.
capPi2 = 1./(1+abs(ones(nSuppX,1)*(1:nSuppX)-(1:nSuppX)'*ones(1,nSuppX)));
capPi2 = capPi2./(sum(capPi2,2)*ones(1,nSuppX));

parameter.change = zeros(8,3);
parameter.change(1,:) = [2.5,0.8,3.8]; % change theta1
parameter.change(2,:) = [0.7,1.2,2.2];  % change theta2
parameter.change(3,:) = [0,1.7,3.7]; % change theta3
parameter.change(4,:) = [1,1.15,2.15]; % change theta4
parameter.change(5,:) = [0.8, 0.55, 0.15]; % change delta
parameter.change(6,:) = [0.4, 0.5, 0.1]; % change beta
parameter.change(7,:) = [1, 0.6, 0.3]; % change betatilde
parameter.change(8,:) = [0.57721, 0, 1]; % change gamma

% parameter name list:
% NamesPara ={'theta1','theta2','theta3','theta4','delta','beta','betatilde','gamma'};
parameterlist = parameter.change(:,1)';

rng(430)
for r = 1:1
    % initialize the parameter list:
    parameterlist = parameter.change(:,1)';
    for j = 1:1

        parameterlist(r) = parameter.change(r,j);
        
        theta = parameterlist(1:4);
        delta = parameterlist(5);
        beta = parameterlist(6);
        betatilde = parameterlist(7);
        gamma = parameterlist(8);

        % just for flow utility, the subscript 1,2,3 are the actions in the current period.
        % action 1 is the terminating action
        [u1,u2] = flowpayoffs(supportX,theta); 
            
        % terminal period choice-specific value function: 
        W_terminal = [u1,u2];
        % generating CCP by extreme value distribution: 
        CCP_terminal = CCP_generator(W_terminal);    
        % terminal period value function (not choice specific, i.e. ex ante)
        V_terminal = W_terminal(:,1) + gamma-log(CCP_terminal(:,1)); % action 1 is the terminating action
        
        % backward induction to derive CCP for each period
        % for exponential case, beta = 1
        [CCP_total,V_total,W_total,~] = backInduction(CCP_terminal,V_terminal, ...
                                                u1,u2,gamma,capPi2,1,1,delta,nPeriods);
        % for sophisticated case, beta is set to be 0.9
        [CCP_total_sophi,V_total_sophi,W_total_sophi,~] = backInduction(CCP_terminal,V_terminal, ...
                                                u1,u2,gamma,capPi2,beta,beta,delta,nPeriods);
        % for naive case, beta is set to be 0.9 and betatilde is set to be 0.97
        [CCP_total_naive,V_total_naive,W_total_naive,CCP_total_tilde] = backInduction(CCP_terminal,V_terminal, ...
                                                 u1,u2,gamma,capPi2,beta,betatilde,delta,nPeriods);
         
        figure(1)
        % Export data for creating graph in Python.
        CCP_data = [CCP_total, CCP_total_sophi, CCP_total_naive];
        save(project_paths('OUT_ANALYSIS','fig1Data.mat'),'CCP_data');

    end
end

%% verification part - fully naive agent
disp("naive agent verification")
% calculate the moment condition at the true beta and delta
z_momentvalue = verification_naive(beta,delta,supportX,capPi2,gamma,theta,CCP_total_naive,CCP_total_tilde,V_total_naive,@flowpayoffs);
epsil = randn(size(CCP_total_naive,1),1).*1E-16;
CCP_total_Natemp(:,1) = CCP_total_naive(:,1) + epsil;
CCP_total_Natemp(:,2) = CCP_total_naive(:,2) - epsil;
temp_Natemp = abs(CCP_total_Natemp-CCP_total_naive)./CCP_total_naive;

objFun_naive_moment = @(bd)verification_naive(bd(1),bd(2),supportX,capPi2,gamma,theta,CCP_total_Natemp,CCP_total_tilde,V_total_naive,@flowpayoffs);
initialv = [0.1,0.01];
[argGMMEst,fv,exitflag,~] = fmincon(objFun_naive_moment,initialv,[],[],[],[],[0,0],[1,1],[]);
disp("disturbed true CCP result")
disp(argGMMEst)
disp("------------------------------");

%% verification part - sophisticated agent
disp("sophi agent verification")
verification = 1;

[z_true,u1_verify,M2_t,M5_t,dphi_true,m1_true] = verification_sophi(supportX,nPeriods,capPi2,beta,delta,gamma,theta,CCP_total_sophi,V_total_sophi,verification, @flowpayoffs);

verification = 0;
% Small disturbance to CCP, check whether the minimum distance estimator is consistent. 
epsil = randn(size(CCP_total_sophi,1),1).*1E-3;
CCP_total_sTtemp(:,1) = CCP_total_sophi(:,1) + epsil;
CCP_total_sTtemp(:,2) = CCP_total_sophi(:,2) - epsil;
temp_sTtemp = abs(CCP_total_sTtemp-CCP_total_sophi)./CCP_total_sophi;

% Check whether matrix 2 and matrix 5 in assumptions 4 & 5 induced by
% disturbed CCP is close to those from true CCP.
[z_disturb, u1_disturb,M2_tDisturbed,M5_tDisturbed,dphi_disturbed,m1_disturbed] = verification_sophi(supportX,nPeriods,capPi2,beta,delta,gamma,theta,CCP_total_sTtemp,V_total_sophi,verification,@flowpayoffs);
M2check = M2_t - M2_tDisturbed;
M5check = M5_t - M5_tDisturbed;
dphi_check = dphi_true-dphi_disturbed;
m1_check = m1_true - m1_disturbed;

% Given different combinations of beta and delta, check whether u1 from different periods
% of data are the same. If they are the same any way, then beta and delta are not identified.
bband = 0.2;
dband = 0.2;
ngrid = 101;
beta0 = linspace(beta-bband,beta+bband,ngrid);
delta0 = linspace(delta-dband,delta+dband,ngrid);
u1_123 = zeros(nSuppX,length(beta0)*length(delta0));
u1_234 = zeros(nSuppX,length(beta0)*length(delta0));
z = ones(length(beta0),length(delta0));
u_index = (find(abs(beta0-beta)<0.001)-1)*length(delta0)+find(abs(delta0-delta)<0.001);
disp(["truevalue is",num2str(u_index),"th column"]);
for i = 1:length(beta0)
    for j = 1:length(delta0)     
    % 1. Criterion distance for different beta and delta (sophisticated
    % agent with 4 consecutive periods)
    [z_P4(i,j),u1_bd,~,~,~,~] = verification_sophi(supportX,nPeriods,capPi2,beta0(i),delta0(j),gamma,theta,CCP_total_sophi,V_total_sophi,verification, @flowpayoffs); 
    % 2. Criterion distance for different beta and delta (sophisticated
    % agent with last three peirods)
    [z_L3(i,j),~] = EstMD_LastPeriods(nPeriods,supportX,capPi2,beta0(i)*delta0(j),delta0(j),gamma,CCP_total_sophi);
    % 3. Criterion distance for different beta and delta (naive agent with last three periods)
    [z_naive_L3(i,j)] = verification_naive(beta0(i),delta0(j),supportX,capPi2,gamma,theta,CCP_total_naive,CCP_total_tilde,V_total_naive,@flowpayoffs);
    end
end
z_P4 = sqrt(z_P4);
z_L3 = sqrt(z_L3);
z_naive_L3 = sqrt(z_naive_L3);

% Output data for graph formatting in Python.
[bb,dd] = meshgrid(beta0,delta0);
L3_data = cat(3,bb,dd,z_L3');
save(project_paths('OUT_ANALYSIS','L3_data.mat'),'L3_data');
P4_data = cat(3,bb,dd,z_P4');
save(project_paths('OUT_ANALYSIS','P4_data.mat'),'P4_data');
L3naive_data= cat(3,bb,dd,z_naive_L3');
save(project_paths('OUT_ANALYSIS','L3naive_data.mat'),'L3naive_data');
 
 
% solve only for beta and delta using CCP
objFun = @(bd)verification_sophi(supportX,nPeriods,capPi2,bd(1),bd(2),gamma,theta,CCP_total_sophi,V_total_sophi,verification, @flowpayoffs); % bd(1) is beta, bd(2) is delta
initialv = [beta,delta];
Options = optimset('Display','off','Algorithm','interior-point','AlwaysHonorConstraints','bounds',...
                                   'GradObj','off','TolFun',1E-10,'TolX',1E-10,'DerivativeCheck','off');

 [argminCriterian,fv,exitflag,~] = fminsearch(objFun,initialv);
disp("disturbed true CCP result")
disp(argminCriterian)

% solve only for beta and delta using true CCP but also imposing the
% normalization while the true theta is not zero
theta_normalization = [0,0,theta(3),theta(4)];
objFun = @(bd)verification_sophi(supportX,nPeriods,capPi2,bd(1),bd(2),gamma,theta_normalization,CCP_total_sophi,V_total_sophi,verification, @flowpayoffs); % bd(1) is beta, bd(2) is delta
initialv = [beta+0.2,delta-0.3];
Options = optimset('Display','off','Algorithm','interior-point','AlwaysHonorConstraints','bounds',...
                                   'GradObj','off','TolFun',1E-10,'TolX',1E-10,'DerivativeCheck','off');
[argminCriterian_normalization,fv_norm,exitflag_norm,~] = fminsearch(objFun,initialv);
disp("DF using true CCP and wrong normalization")
disp(argminCriterian_normalization)
% the results are the same


% solve for beta, delta jointly with theta
objFun = @(bd)verification_sophi(supportX,nPeriods,capPi2,bd(1),bd(2),gamma,[bd(3),bd(4),theta(3:4)],...
    CCP_total_sophi,V_total_sophi,verification, @flowpayoffs); % bd(1) is beta, bd(2) is delta, bd(3) is theta1, bd(4) is theta2
initialv = [beta-0.2,delta+0.2,0.1,0.1];
Options = optimset('Display','off','Algorithm','interior-point','AlwaysHonorConstraints','bounds',...
                                   'GradObj','off','TolFun',1E-8,'TolX',1E-10,'DerivativeCheck','off');
[argminCriterian_together,fval,exitflag2,~] = fmincon(objFun,initialv,[],[],[],[],[0,0,-Inf,-Inf],[1,1,Inf,Inf],[]);
disp("disturbed CCP result together")
disp(argminCriterian_together)

%% Estimation
disp('Estimates Starts')
%% Estimate from the last period
% Try to estimate the model using the last three periods 
% minimum distance estimation
[z_true_L3, u1_true_L3]=EstMD_LastPeriods(nPeriods,supportX,capPi2,beta,delta,gamma,CCP_total_sophi);
objFun = @(bd)EstMD_LastPeriods(nPeriods,supportX,capPi2,bd(1),bd(2),gamma,CCP_total_sTtemp); % bd(1) is beta, bd(2) is delta
initialv = [beta-0.001,delta-0.001];
Options = optimset('Display','off','Algorithm','interior-point','AlwaysHonorConstraints','bounds',...
                                   'GradObj','off','TolFun',1E-10,'TolX',1E-10,'DerivativeCheck','off');
[argminCriterian,fv,exitflag,~] = fminsearch(objFun,initialv,Options);
disp("disturbed true CCP Last 3Periods")
disp("beta delta")
argminshow = [argminCriterian(1),argminCriterian(2)];
disp(argminshow);

tic
for N_size = 1:6
    nFirms = nFirms_alternatives(N_size);
parfor j = 1:maxj

    [choices_sophi,iX_sophi] = simulateData_CCP(CCP_total_sophi,capPi2,nPeriods,nFirms);
    [choices_naive,iX_naive] = simulateData_CCP(CCP_total_naive,capPi2,nPeriods,nFirms);
    % recover CCP practice from the data generated
    CCP_recover_sophi = RecoverCCP(choices_sophi,iX_sophi,nPeriods,nSuppX,nA);
    CCP_recover_naive = RecoverCCP(choices_naive,iX_naive,nPeriods,nSuppX,nA);
    % OLS estimator from last three periods data
    EstimatesOLS_L3(:,j) = EstOLS_LastPeriods(nPeriods,supportX,gamma,capPi2,CCP_recover_sophi);

    %% Minimum Distance estimation for complete naive agent case
    objFun_naive_moment = @(bd)verification_naive(exp(bd(1))/(1+exp(bd(1))),exp(bd(2))/(1+exp(bd(2))),supportX,capPi2,gamma,theta,CCP_recover_naive,CCP_total_tilde,V_total_naive,@flowpayoffs);
    initialv = [log(0.1/0.9),log(0.6/0.4)];
    Options = optimset('Display','off','Algorithm','interior-point','AlwaysHonorConstraints','bounds',...
                                    'GradObj','off','TolFun',1E-8,'TolX',1E-10,'DerivativeCheck','off');
    [EstimatesMoment_naive(:,j),fv,exitflag_MD_naive(:,j)] = fminsearch(objFun_naive_moment,initialv,Options);
    EstimatesMoment_naive(:,j) = exp(EstimatesMoment_naive(:,j))./(1+exp(EstimatesMoment_naive(:,j)));


    %% Maximum Likelihood Estimation: This is probably the most relevant approach for applied work.
    if nFirms <= 20001
        [Est_sDsM(:,j),~,exitflag_mle(:,j)] = MLEstimation(choices_sophi,iX_sophi,supportX,nA,capPi2,gamma,2,theta,thetaY,...
                [],@nLogLik,@flowpayoffs,@CCP_generator,@backInduction);
        [Est_sDsM_normalization_terminating(:,j),~,exitflag_mle_norm_t(:,j)] = MLEstimation(choices_sophi,iX_sophi,supportX,nA,capPi2,gamma,2,theta,thetaNormalizationT,...
                [],@nLogLik,@flowpayoffs,@CCP_generator,@backInduction);
        [Est_sDsM_normalization(:,j),~,exitflag_mle_norm(:,j)] = MLEstimation(choices_sophi,iX_sophi,supportX,nA,capPi2,gamma,2,theta,thetaNormalization,...
                [],@nLogLik,@flowpayoffs,@CCP_generator,@backInduction);
        [Est_nDnM(:,j),~,exitflag_mle_naive(:,j)] = MLEstimation(choices_naive,iX_naive,supportX,nA,capPi2,gamma,3,theta,thetaY,...
                [],@nLogLik,@flowpayoffs,@CCP_generator,@backInduction);
    end
    disp(['iter:' num2str(j)])
end
clear I_eDeM I_eDsM I_sDeM I_sDsM
clear MLE_eDeM MLE_eDsM MLE_sDeM MLE_sDsM

% Run the counterfactuals
if nFirms == 10000
    Est_counter = Est_sDsM;
    Est_counter_norm = Est_sDsM_normalization;
    Est_counter_norm_terminating = Est_sDsM_normalization_terminating;
end

Est_naiveMD = mean(EstimatesMoment_naive,2);
Est_sophiOLS_L3 = mean(EstimatesOLS_L3,2);
% compute the mean and std deviation, and then do ttest
para_sDsM = mean(Est_sDsM,2);
para_sDsM_normalization = mean(Est_sDsM_normalization,2);
para_sDsM_normalization_terminating = mean(Est_sDsM_normalization_terminating,2);
para_nDnM = mean(Est_nDnM,2);
SD_sDsM = std(Est_sDsM,0,2);
SD_sDsM_normalization = std(Est_sDsM_normalization,0,2);
SD_sDsM_normalization_terminating = std(Est_sDsM_normalization_terminating,0,2);
SD_nDnM = std(Est_nDnM,0,2);
averageTime = toc;
disp(['total time = ', num2str(toc/60),'m']);
SD_sophiOLS_L3 = std(EstimatesOLS_L3,0,2);
SD_naiveMD = std(EstimatesMoment_naive,0,2);
% save the parameters estimated in file
% M means model specified, D means data used
% sMsD means estimating parameters of sophisticated agents with
% sophisticated agent data
% eMeD means estimating parameters of exponential agents with
% actually sophisticated agent data
% save("parameterlist_counter1.mat",'para_mean');
% save('estimatesSE.mat','SD_sDsM')

% Here we basically test two cases: case 1, the underlying model is present
bdtrue = [beta;delta];
if exist(project_paths('OUT_TABLES','bdEst.mat'),'file')~=2
    save(project_paths('OUT_TABLES','bdEst.mat'),"bdtrue")
end
if nFirms == 100000
     Est100000_naive = Est_naiveMD;
     SD100000_naive = SD_naiveMD;
     Est100000_OLS_L3 = Est_sophiOLS_L3;
     SD100000_OLS_L3 = SD_sophiOLS_L3;
elseif nFirms == 10000
     Est10000_naive = Est_naiveMD;
     SD10000_naive = SD_naiveMD;
     Est10000_mle = para_sDsM;
     SD10000_mle = SD_sDsM;
     Est10000_mle_norm = para_sDsM_normalization;
     SD10000_mle_norm = SD_sDsM_normalization;
     Est10000_mle_norm_terminating = para_sDsM_normalization_terminating;
     SD10000_mle_norm_terminating = SD_sDsM_normalization_terminating;
     Est10000_mle_naive = para_nDnM;
     SD10000_mle_naive = SD_nDnM;
     Est10000_OLS_L3 = Est_sophiOLS_L3;
     SD10000_OLS_L3 = SD_sophiOLS_L3;
elseif nFirms == 50000
     Est50000_naive = Est_naiveMD;
     SD50000_naive = SD_naiveMD;
     Est50000_OLS_L3 = Est_sophiOLS_L3;
     SD50000_OLS_L3 = SD_sophiOLS_L3;
elseif nFirms == 1000000
     Est1000000_naive = Est_naiveMD;
     SD1000000_naive = SD_naiveMD;
     Est1000000_OLS_L3 = Est_sophiOLS_L3;
     SD1000000_OLS_L3 = SD_sophiOLS_L3;
elseif nFirms == 5000
     Est5000_naive = Est_naiveMD;
     SD5000_naive = SD_naiveMD;
     Est5000_mle = para_sDsM;
     SD5000_mle = SD_sDsM;
     Est5000_mle_norm = para_sDsM_normalization;
     SD5000_mle_norm = SD_sDsM_normalization;
     Est5000_mle_norm_terminating = para_sDsM_normalization_terminating;
     SD5000_mle_norm_terminating = SD_sDsM_normalization_terminating;
     Est5000_mle_naive = para_nDnM;
     SD5000_mle_naive = SD_nDnM;
     Est5000_OLS_L3 = Est_sophiOLS_L3;
     SD5000_OLS_L3 = SD_sophiOLS_L3;
elseif nFirms == 20000
    Est20000_mle = para_sDsM;
    SD20000_mle = SD_sDsM;
    Est20000_mle_norm = para_sDsM_normalization;
    SD20000_mle_norm = SD_sDsM_normalization;
    Est20000_mle_norm_terminating = para_sDsM_normalization_terminating;
    SD20000_mle_norm_terminating = SD_sDsM_normalization_terminating;
    Est20000_mle_naive = para_nDnM;
    SD20000_mle_naive = SD_nDnM;
end



end

%% counterfactual analysis
% comparing the case with and without normalization under the tax on not
% adopting
nPeriods_counter = 10;
CCP_total_sophi_longer = backInduction(CCP_terminal,V_terminal,... 
                             u1,u2,gamma,capPi2,beta,beta,delta,nPeriods_counter);

theta_counter = Est10000_mle(3:end) - [0;0;0.5;0]; % full estimation with counterfactual
beta_counter = Est10000_mle(1);
delta_counter = Est10000_mle(2);
[u1_counter,u2_counter] = flowpayoffs(supportX,theta_counter);
W_terminal_non = [u1_counter,u2_counter];
CCP_terminal_non = CCP_generator(W_terminal_non);    
V_terminal_non = W_terminal_non(:,1) + gamma-log(CCP_terminal_non(:,1)); % action 1 is the terminating action      
[CCP_counter,~,~,~] = backInduction(CCP_terminal_non,V_terminal_non,... 
                             u1_counter,u2_counter,gamma,capPi2,beta_counter,beta_counter,delta_counter,nPeriods_counter);
                         
theta_counter_normTerm = [0;0;Est10000_mle_norm_terminating(3:end)] - [0;0;0.5;0]; 
beta_counter_normTerm = Est10000_mle_norm_terminating(1);
delta_counter_normTerm = Est10000_mle_norm_terminating(2);
[u1_counter_normTerm,u2_counter_normTerm] = flowpayoffs(supportX,theta_counter_normTerm);
W_terminal_normTerm = [u1_counter_normTerm,u2_counter_normTerm];
CCP_terminal_normTerm = CCP_generator(W_terminal_normTerm);    
V_terminal_normTerm = W_terminal_normTerm(:,1) + gamma-log(CCP_terminal_normTerm(:,1)); % action 1 is the terminating action      
[CCP_counter_normTerm,~,~,~] = backInduction(CCP_terminal_normTerm,V_terminal_normTerm,... 
                             u1_counter_normTerm,u2_counter_normTerm,gamma,capPi2,beta_counter_normTerm,beta_counter_normTerm,delta_counter_normTerm,nPeriods_counter);
legendc = {'true model','cf w/o norm','cf w/ norm'};
% Export data for creating graph in Python.
CCP_counter_normterm_data = [CCP_total_sophi_longer,CCP_counter,CCP_counter_normTerm];
save(project_paths('OUT_ANALYSIS','CCPCounterNormTermData.mat'),'CCP_counter_normterm_data');


theta_counter = Est10000_mle(3:end) + [0.5;0;0;0]; % full estimation with counterfactual on terminating action
beta_counter = Est10000_mle(1);
delta_counter = Est10000_mle(2);
[u1_counter,u2_counter] = flowpayoffs(supportX,theta_counter);
W_terminal_non = [u1_counter,u2_counter];
CCP_terminal_non = CCP_generator(W_terminal_non);    
V_terminal_non = W_terminal_non(:,1) + gamma-log(CCP_terminal_non(:,1)); % action 1 is the terminating action      
[CCP_counter,~,~,~] = backInduction(CCP_terminal_non,V_terminal_non,... 
                             u1_counter,u2_counter,gamma,capPi2,beta_counter,beta_counter,delta_counter,nPeriods_counter);

theta_counter_norm = [Est10000_mle_norm(3:end);0;0] + [0.5;0;0;0];
beta_counter_norm = Est10000_mle_norm(1);
delta_counter_norm = Est10000_mle_norm(2);
[u1_counter_norm,u2_counter_norm] = flowpayoffs(supportX,theta_counter_norm);
W_terminal_norm = [u1_counter_norm,u2_counter_norm];
CCP_terminal_norm = CCP_generator(W_terminal_norm);    
V_terminal_norm = W_terminal_norm(:,1) + gamma-log(CCP_terminal_norm(:,1)); % action 1 is the terminating action      
[CCP_counter_norm,~,~,~] = backInduction(CCP_terminal_norm,V_terminal_norm,... 
                             u1_counter_norm,u2_counter_norm,gamma,capPi2,beta_counter_norm,beta_counter_norm,delta_counter_norm,nPeriods_counter);
legendc = {'true model','cf w/o norm','cf w/ norm'};
CCP_total_sophi_longer = backInduction(CCP_terminal,V_terminal,... 
                             u1,u2,gamma,capPi2,beta,beta,delta,nPeriods_counter);
% Export data for creating graph in Python.
CCP_counter_data = [CCP_total_sophi_longer,CCP_counter,CCP_counter_norm];
save(project_paths('OUT_ANALYSIS','CCPCounterData.mat'),'CCP_counter_data');

% creating the confidence interval (95%), just choose one action to print: here I choose adoption decision                       
CCP_sophi_counter_CI = zeros(nPeriods_counter*nSuppX,maxj);
CCP_sophi_counter_CI_u2 = zeros(nPeriods_counter*nSuppX,maxj);
CCP_sophi_counter_norm_terminating_CI = zeros(nPeriods_counter*nSuppX,maxj);
CCP_sophi_counter_norm_CI = zeros(nPeriods_counter*nSuppX,maxj);
for i = 1:maxj
    % counterfactual for the correct model (change utility on non-terminating action)
    theta_counter_i = Est_counter(3:end,i) - [0;0;0.5;0]; % full estimation with counterfactual
    beta_counter_i = Est_counter(1);
    delta_counter_i = Est_counter(2);
    [u1_counter_i,u2_counter_i] = flowpayoffs(supportX,theta_counter_i);
    W_terminal_non_i = [u1_counter_i,u2_counter_i];
    CCP_terminal_non_i = CCP_generator(W_terminal_non_i);    
    V_terminal_non_i = W_terminal_non_i(:,1) + gamma-log(CCP_terminal_non_i(:,1)); % action 1 is the terminating action      
    [CCP_sophi_counter,~,~,~] = backInduction(CCP_terminal_non_i,V_terminal_non_i,... 
                                 u1_counter_i,u2_counter_i,gamma,capPi2,beta_counter_i,beta_counter_i,delta_counter_i,nPeriods_counter);
    CCP_sophi_counter_CI(:,i) = CCP_sophi_counter(:,1); % just choose action 1 to execute cf

    % counterfactual for the model with normalization on terminating action
    % then the counterfactual is done on non-terminating action utility
    theta_counter_norm_terminating_i = [0;0;Est_counter_norm_terminating(3:end,i)] - [0;0;0.5;0]; % full estimation with counterfactual
    beta_counter_norm_terminating_i = Est_counter_norm_terminating(1);
    delta_counter_norm_terminating_i = Est_counter_norm_terminating(2);
    [u1_counter_norm_terminating_i,u2_counter_norm_terminating_i] = flowpayoffs(supportX,theta_counter_norm_terminating_i);
    W_terminal_non_i = [u1_counter_norm_terminating_i,u2_counter_norm_terminating_i];
    CCP_terminal_non_i = CCP_generator(W_terminal_non_i);    
    V_terminal_non_i = W_terminal_non_i(:,1) + gamma-log(CCP_terminal_non_i(:,1)); % action 1 is the terminating action      
    [CCP_sophi_counter,~,~,~] = backInduction(CCP_terminal_non_i,V_terminal_non_i,... 
                                 u1_counter_norm_terminating_i,u2_counter_norm_terminating_i,gamma,capPi2,beta_counter_norm_terminating_i,beta_counter_norm_terminating_i,delta_counter_norm_terminating_i,nPeriods_counter);
    CCP_sophi_counter_norm_terminating_CI(:,i) = CCP_sophi_counter(:,1); % just choose action 1 to execute cf

    
    % counterfactual for the correct model (change utility on terminating action)
    theta_counter_i = Est_counter(3:end,i) + [0.5;0;0;0]; % full estimation with counterfactual
    beta_counter_i = Est_counter(1);
    delta_counter_i = Est_counter(2);
    [u1_counter_i,u2_counter_i] = flowpayoffs(supportX,theta_counter_i);
    W_terminal_non_i = [u1_counter_i,u2_counter_i];
    CCP_terminal_non_i = CCP_generator(W_terminal_non_i);    
    V_terminal_non_i = W_terminal_non_i(:,1) + gamma-log(CCP_terminal_non_i(:,1)); % action 1 is the terminating action      
    [CCP_sophi_counter,~,~,~] = backInduction(CCP_terminal_non_i,V_terminal_non_i,... 
                                 u1_counter_i,u2_counter_i,gamma,capPi2,beta_counter_i,beta_counter_i,delta_counter_i,nPeriods_counter);
    CCP_sophi_counter_CI_u2(:,i) = CCP_sophi_counter(:,1); % just choose action 1 to execute cf

    % counterfactual for the model with normalization on non-terminating action
    % then the counterfactual is done on terminating action utility
    theta_counter_norm_i = [Est_counter_norm(3:end,i);0;0] + [0.5;0;0;0]; % full estimation with counterfactual
    beta_counter_norm_i = Est_counter_norm(1);
    delta_counter_norm_i = Est_counter_norm(2);
    [u1_counter_norm_i,u2_counter_norm_i] = flowpayoffs(supportX,theta_counter_norm_i);
    W_terminal_non_i = [u1_counter_norm_i,u2_counter_norm_i];
    CCP_terminal_non_i = CCP_generator(W_terminal_non_i);    
    V_terminal_non_i = W_terminal_non_i(:,1) + gamma-log(CCP_terminal_non_i(:,1)); % action 1 is the terminating action      
    [CCP_sophi_counter,~,~,~] = backInduction(CCP_terminal_non_i,V_terminal_non_i,... 
                                 u1_counter_norm_i,u2_counter_norm_i,gamma,capPi2,beta_counter_norm_i,beta_counter_norm_i,delta_counter_norm_i,nPeriods_counter);
    CCP_sophi_counter_norm_CI(:,i) = CCP_sophi_counter(:,1); % just choose action 1 to execute cf

end

% get the confidence band
% for each band, save the mean (first column), lower bound (second column), and upper bound (thrid column)
CCP_sophi_counter_band = zeros(nPeriods_counter*nSuppX,3);
CCP_sophi_counter_band_u2 = zeros(nPeriods_counter*nSuppX,3);
CCP_sophi_counter_norm_terminating_band = zeros(nPeriods_counter*nSuppX,3);
CCP_sophi_counter_norm_band = zeros(nPeriods_counter*nSuppX,3);
conf = 0.95;
alpha = 1 - conf;
pLo = alpha/2;
pUp = 1 - alpha/2;
for xt=1:nPeriods_counter*nSuppX
    SEM = std(CCP_sophi_counter_CI(xt,:));
    ts = tinv([pLo,pUp],length(CCP_sophi_counter_CI(xt,:))-1);
    CCP_sophi_counter_band(xt,1) = mean(CCP_sophi_counter_CI(xt,:));
    CCP_sophi_counter_band(xt,2:3) = mean(CCP_sophi_counter_CI(xt,:)) + ts*SEM; 
    SEM = std(CCP_sophi_counter_norm_terminating_CI(xt,:));
    ts = tinv([pLo,pUp],length(CCP_sophi_counter_norm_terminating_CI(xt,:))-1);
    CCP_sophi_counter_norm_terminating_band(xt,1) = mean(CCP_sophi_counter_norm_terminating_CI(xt,:));
    CCP_sophi_counter_norm_terminating_band(xt,2:3) = mean(CCP_sophi_counter_norm_terminating_CI(xt,:)) + ts*SEM; 
    SEM = std(CCP_sophi_counter_CI_u2(xt,:));
    ts = tinv([pLo,pUp],length(CCP_sophi_counter_CI_u2(xt,:))-1);
    CCP_sophi_counter_band_u2(xt,1) = mean(CCP_sophi_counter_CI_u2(xt,:));
    CCP_sophi_counter_band_u2(xt,2:3) = mean(CCP_sophi_counter_CI_u2(xt,:)) + ts*SEM; 
    SEM = std(CCP_sophi_counter_norm_CI(xt,:));
    ts = tinv([pLo,pUp],length(CCP_sophi_counter_norm_CI(xt,:))-1);
    CCP_sophi_counter_norm_band(xt,1) = mean(CCP_sophi_counter_norm_CI(xt,:));
    CCP_sophi_counter_norm_band(xt,2:3) = mean(CCP_sophi_counter_norm_CI(xt,:)) + ts*SEM; 
end
% Export data for creating graph in Python.
CCP_counter_normTerm_banddata = [CCP_total_sophi_longer(:,1),CCP_sophi_counter_band,CCP_sophi_counter_norm_terminating_band];
save(project_paths('OUT_ANALYSIS','CCPCounterNormTermBandData.mat'),'CCP_counter_normTerm_banddata');
CCP_counter_banddata = [CCP_total_sophi_longer(:,1),CCP_sophi_counter_band_u2,CCP_sophi_counter_norm_band];
save(project_paths('OUT_ANALYSIS','CCPCounterBandData.mat'),'CCP_counter_banddata');

% Export relevant data for table formatting in Python.
MD_table_labels = {'true','beta_hat_50000','sd_50000','beta_hat_100000','sd_100000','beta_hat_1000000','sd_1000000'};
estMDnaive_table = table(bdtrue,Est50000_naive,SD50000_naive,Est100000_naive,SD100000_naive,Est1000000_naive,SD1000000_naive, 'RowNames',{'$\beta$','$\delta$'},'VariableNames',MD_table_labels);
writetable(estMDnaive_table, project_paths('OUT_ANALYSIS','md_table_data_naive.csv'),'writeRowNames',true);
OLS_table_labels = {'true','beta_hat_50000','sd_50000','beta_hat_100000','sd_100000','beta_hat_1000000','sd_1000000'};
estOLS_table = table([bdtrue;u1],Est50000_OLS_L3,SD50000_OLS_L3,Est100000_OLS_L3,SD100000_OLS_L3,Est1000000_OLS_L3,SD1000000_OLS_L3, 'RowNames',{'$\beta$','$\delta$','$u_{11}$','$u_{12}$','$u_{13}$','$u_{14}$'},'VariableNames',OLS_table_labels);
writetable(estOLS_table, project_paths('OUT_ANALYSIS','ols_table_data.csv'),'writeRowNames',true);
MLE_table_labels = {'true','beta_hat_5000','sd_5000','beta_hat_10000','sd_10000','beta_hat_20000','sd_20000'};
estMLEsophi_table = table([bdtrue;theta'],Est5000_mle,SD5000_mle,Est10000_mle,SD10000_mle,Est20000_mle,SD20000_mle, 'RowNames',{'$\beta$','$\delta$','$\theta_{1}$','$\theta_{2}$','$\theta_{3}$','$\theta_{4}$'},'VariableNames',MLE_table_labels);
writetable(estMLEsophi_table, project_paths('OUT_ANALYSIS','mle_table_data.csv'),'writeRowNames',true);
estMLEsophi_norm_table = table([bdtrue;theta(1:2)'],Est5000_mle_norm,SD5000_mle_norm,Est10000_mle_norm,SD10000_mle_norm,Est20000_mle_norm,SD20000_mle_norm, 'RowNames',{'$\beta$','$\delta$','$\theta_{1}$','$\theta_{2}$'},'VariableNames',MLE_table_labels);
writetable(estMLEsophi_norm_table, project_paths('OUT_ANALYSIS','mle_norm_table_data.csv'),'writeRowNames',true);
estMLEsophi_normTerm_table = table([bdtrue;theta(3:4)'],Est5000_mle_norm_terminating,SD5000_mle_norm_terminating,Est10000_mle_norm_terminating,SD10000_mle_norm_terminating,Est20000_mle_norm_terminating,SD20000_mle_norm_terminating, 'RowNames',{'$\beta$','$\delta$','$\theta_{3}$','$\theta_{4}$'},'VariableNames',MLE_table_labels);
writetable(estMLEsophi_normTerm_table, project_paths('OUT_ANALYSIS','mle_normTerm_table_data.csv'),'writeRowNames',true);
estMLEnaive_table = table([bdtrue;theta'],Est5000_mle_naive,SD5000_mle_naive,Est10000_mle_naive,SD10000_mle_naive,Est20000_mle_naive,SD20000_mle_naive, 'RowNames',{'$\beta$','$\delta$','$\theta_{1}$','$\theta_{2}$','$\theta_{3}$','$\theta_{4}$'},'VariableNames',MLE_table_labels);
writetable(estMLEnaive_table, project_paths('OUT_ANALYSIS','mle_naive_table_data.csv'),'writeRowNames',true);
% Save post-simulation workspace.
save(project_paths('OUT_ANALYSIS','postestimation_DDC_oneTA.mat'));
