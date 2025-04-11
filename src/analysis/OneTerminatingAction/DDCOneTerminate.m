% Main script to run all the simulations discussed in the paper.
clear
clc

%% Set parameters of simulations.
maxj = 100; % number of maximum simulation runs.
nPeriods = 4; % number of periods.
nFirms_alternatives =[10000,50000,100000,1000000,20000,5000,30000]; % sample size.
nA = 2; % number of actions.
  
% Define parameter vector.
thetaY = 1; % if this  equals one, full theta is estimated (no normalization).
thetaNormalizationT = 2; % if this equals to 2, normalized theta (normalizing the terminating action) is estimated.
thetaNormalization = 3; % if this equals to 3, normalized theta (normalizing the other action) is estimated.
theta_rn = 4; % if this equals 4, theta (making utility on x1 and x2 the same) is estimated with additional exclusion restrictions.
theta_rn_norm = 5; % if this equals 5, normalized theta (normalizing terminating action) is estimated  using exclusion restrictions.
theta_rn_norm1 = 6; % if this equals 6, normalized theta on only one value of x2 (normalizing on one state value) is estimated.

% Set parameters of utility function.
utilityx2=[2.1 3.5]';  
utilityx1=[4.2 3.8]';
supportX = ["$x_1 = x_{1L}, x_2 = x_{2L}$","$x_1 = x_{1L}, x_2 = x_{2H}$","$x_1 = x_{1H}, x_2 = x_{2L}$","$x_1 = x_{1H}, x_2 = x_{2H}$"]';

% Transition matrix for price quality adjustment measure.
Transit_x1=[0.6 0.4
            0.1 0.9];
% Transition matrix for "electriciyy price".
Transit_x2=[0.8 0.2
            0.2 0.8];
% Construct the combinations of (x1 x2) as (L L;L H; H L; H H)
capPi2=kron(Transit_x1,Transit_x2);

% Define the parameters.
parameter.change = zeros(8,3);
parameter.change(1,:) = [2.5,0.8,3.8]; % change theta1
parameter.change(2,:) = [0.7,1.2,2.2];  % change theta2
parameter.change(3,:) = [0,1.7,3.7]; % change theta3
parameter.change(4,:) = [1,1.15,2.15]; % change theta4
parameter.change(5,:) = [0.8, 0.55, 0.15]; % change delta
parameter.change(6,:) = [0.6, 0.6, 0.5]; % change beta
parameter.change(7,:) = [1, 0.6, 0.3]; % change betatilde
parameter.change(8,:) = [0,0.57721,1]; % change gamma
% Parameter name list.
parameterlist = parameter.change(:,1)';
% Set random seed to ensure simulation is reproducible.
rng(430)
% Initialize the parameter list:
parameterlist = parameter.change(:,1)';
% Extract parameter components.
% parameterlist(1) = parameter.change(1,1);
theta = [utilityx1;utilityx2]';
delta = parameterlist(5);
beta = parameterlist(6);
betatilde = parameterlist(7);
gamma = parameterlist(8);
% For the flow utility, the subscript 1, 2, 3 indicates the action in the current period. Action 1 is the terminating action
u1_ext=kron(utilityx1,ones(2,1));
u2_ext=repmat(utilityx2,2,1);
% Recycle flow payoffs function for use in MLE estimation.
[u1,u2] = flowpayoffs(supportX,[u1_ext;u2_ext]);
nSuppX =size(u1,1);
% Terminal period choice-specific value function. 
W_terminal = [u1,u2];
% Calculate CCPs for model with extreme value distributed error term.
CCP_terminal = CCP_generator(W_terminal);    
% Terminal period ex-ante value function (not choice-specific)
V_terminal = W_terminal(:,1) + gamma-log(CCP_terminal(:,1)); % action 1 is the terminating action.
% Backward induction to derive CCPs for each period.
% For exponential discounting, beta = 1.
[CCP_total,V_total,W_total,~] = backInduction(CCP_terminal,V_terminal, ...
                                        u1,u2,gamma,capPi2,1,1,delta,nPeriods);
% For sophisticated agent, beta is set to 0.9.
[CCP_total_sophi,V_total_sophi,W_total_sophi,~] = backInduction(CCP_terminal,V_terminal, ...
                                        u1,u2,gamma,capPi2,beta,beta,delta,nPeriods);
% For naive case, beta is set to 0.9 and betatilde is set to 0.97.
[CCP_total_naive,V_total_naive,W_total_naive,CCP_total_tilde] = backInduction(CCP_terminal,V_terminal, ...
                                            u1,u2,gamma,capPi2,beta,betatilde,delta,nPeriods);
% Export data for creating graph in Python.
CCP_data = [CCP_total, CCP_total_sophi, CCP_total_naive];
save(project_paths('OUT_ANALYSIS','fig1Data.mat'),'CCP_data');

%% Verification  - naive agent
disp("naive agent verification")
% Calculate the moment conditions at the true beta and delta.
[z_momentvalue,~] = verification_naive(beta,delta,supportX,capPi2,gamma,theta,CCP_total_naive,CCP_total_tilde,V_total_naive,@flowpayoffs);
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

%% Verification - sophisticated agent
disp("sophi agent verification")
verification = 1;
theta_verification = [kron(theta(1:2)',ones(2,1));repmat(theta(3:4)',2,1)]; % construct it to be 
[z_true,u1_verify,M2_t,M5_t,dphi_true,m1_true] = verification_sophi(supportX,nPeriods,capPi2,beta,delta,gamma,...
    theta_verification,CCP_total_sophi,V_total_sophi,verification, @flowpayoffs);
verification = 0;
% Introduce a small disturbance to CCP, check whether the minimum distance
% estimator is valid.
epsil = randn(size(CCP_total_sophi,1),1).*1E-3;
CCP_total_sTtemp(:,1) = CCP_total_sophi(:,1) + epsil;
CCP_total_sTtemp(:,2) = CCP_total_sophi(:,2) - epsil;
temp_sTtemp = abs(CCP_total_sTtemp-CCP_total_sophi)./CCP_total_sophi;

% Check whether Matrix 2 and Matrix 5 in assumption 4 & 5 induced by
% disturbed CCP is closed to those from the true CCP.
[z_disturb, u1_disturb,M2_tDisturbed,M5_tDisturbed,dphi_disturbed,m1_disturbed] = verification_sophi(supportX,nPeriods,capPi2,beta,delta,gamma,...
    theta_verification,CCP_total_sTtemp,V_total_sophi,verification,@flowpayoffs);
M2check = M2_t - M2_tDisturbed;
M5check = M5_t - M5_tDisturbed;
dphi_check = dphi_true-dphi_disturbed;
m1_check = m1_true - m1_disturbed;

% Given different combinations of beta and delta, check whether u1 from different periods
% of data are the same. If they are the same anyway, then beta and delta
% are not identified.
bband = 0.2;
dband = 0.2;
ngrid = 101;
beta0 = linspace(beta-bband,beta+bband,ngrid);
delta0 = linspace(delta-dband,delta+dband,ngrid);
u1_123 = zeros(nSuppX,length(beta0)*length(delta0));
u1_234 = zeros(nSuppX,length(beta0)*length(delta0));
z = ones(length(beta0),length(delta0));
u_index = (find(abs(beta0-beta)<0.001)-1)*length(delta0)+find(abs(delta0-delta)<0.001);
disp(["true value is",num2str(u_index),"th column"]);
for i = 1:length(beta0)
    for j = 1:length(delta0)     
        % 1. Criterion distance for different beta and delta (sophisticated
        % agent with 4 consecutive periods)
        [z_P4(i,j),u1_bd,~,~,~,~] = verification_sophi(supportX,nPeriods,capPi2,beta0(i),delta0(j),gamma,theta_verification,CCP_total_sophi,V_total_sophi,verification, @flowpayoffs); 
        % 2. Criterion distance for different beta and delta (sophisticated
        % agent with last three periods)
        [z_L3(i,j),~] = EstMD_LastPeriods(nPeriods,supportX,capPi2,beta0(i),delta0(j),gamma,CCP_total_sophi);
        % 3. Criterion distance for different beta and delta (naive
        % agent with last three periods)
        [z_naive_L3(i,j)] = verification_naive(beta0(i),delta0(j),supportX,capPi2,gamma,theta,CCP_total_naive,CCP_total_tilde,V_total_naive,@flowpayoffs);
     end
 end
 z_P4 = sqrt(z_P4);
 z_L3 = sqrt(z_L3);
 z_naive_L3 = sqrt(z_naive_L3);
 disp(["true u1",u1']);
 disp(["u1_123 at true bd", u1_123(:,u_index)']);
 disp(["u1_234 at true bd", u1_234(:,u_index)']);
 figure(10)
[bb,dd] = meshgrid(beta0,delta0);
 contourf(bb,dd,z_P4',100);
 xlabel("beta");ylabel("delta");
 close(figure(10));
 
 figure(11)
 [bb,dd] = meshgrid(beta0,delta0);
 contourf(bb,dd,z_L3',100);
 xlabel("beta");ylabel("delta");
 close(figure(11));
 
 figure(12)
[bb,dd] = meshgrid(beta0,delta0);
 contourf(bb,dd,z_naive_L3',100);
 xlabel("beta");ylabel("delta");
% Export graph data for formatting in Python.
% This is a 3-d data structure that contains the beta grid in the first
% level, delta grid in the second level, and the actual data in the third
% level.
P4_data = cat(3,bb,dd,z_P4');
save(project_paths('OUT_ANALYSIS','P4_data.mat'),'P4_data');
L3naive_data= cat(3,bb,dd,z_naive_L3');
save(project_paths('OUT_ANALYSIS','L3naive_data.mat'),'L3naive_data');
close(figure(12));

%% Estimation Part
disp('Simulation and estimation starts.')
% Estimation using the last period.
sc = parallel.pool.Constant(RandStream('Threefry','Seed', 161866));
tic
% Loop over sample sizes.
for N_size = 1:length(nFirms_alternatives)
    nFirms = nFirms_alternatives(N_size);
   % Loop over simulation runs.
   parfor j = 1:maxj
         stream = sc.Value;
         stream.Substream = j;
         [choices_expon,iX_expon] = simulateData_CCP(CCP_total,capPi2,nPeriods,nFirms,stream);
         [choices_sophi,iX_sophi] = simulateData_CCP(CCP_total_sophi,capPi2,nPeriods,nFirms,stream);
         [choices_naive,iX_naive] = simulateData_CCP(CCP_total_naive,capPi2,nPeriods,nFirms,stream);
     
    %% recover CCP practice from the data generated
        CCP_recover_expon = RecoverCCP(choices_expon,iX_expon,nPeriods,nSuppX,nA);
        CCP_recover_sophi = RecoverCCP(choices_sophi,iX_sophi,nPeriods,nSuppX,nA);
        CCP_recover_naive = RecoverCCP(choices_naive,iX_naive,nPeriods,nSuppX,nA);
    %% OLS estimator using final three periods data.
    [EstimatesOLS_L3(:,j),SSres_OLS(:,j)] = EstOLS_LastPeriods(nPeriods,supportX,gamma,capPi2,CCP_recover_sophi);
    
    %% Minimum distance estimation     
    % Set initial values for beta and delta.
    initialv=0.9*[log(0.98*beta/(1-0.98*beta)) log(0.98*delta/(1-0.98*delta))];
    % Minimum distance estimation using final three periods of data without exclusion restrictions. 
    objFun = @(bd)EstMD_LastPeriods(nPeriods,supportX,capPi2,exp(bd(1))/(1+exp(bd(1))),exp(bd(2))/(1+exp(bd(2))),gamma,CCP_recover_sophi); % bd(1) is beta, bd(2) is delta
    Options = optimset('Display','off','Algorithm','interior-point','AlwaysHonorConstraints','bounds',...
                                       'GradObj','off','TolFun',1E-10,'TolX',1E-10,'DerivativeCheck','off');
    [EstimatesMoment(:,j),fval,exitflag_ROLS(:,j)] = fminsearch(objFun,initialv,Options);
    EstimatesMoment(:,j) = exp(EstimatesMoment(:,j))./(1+exp(EstimatesMoment(:,j)));
    bdROLS = EstimatesMoment(:,j);
    [~,~,SSres_ROLS(:,j)] = EstMD_LastPeriods(nPeriods,supportX,capPi2,bdROLS(1),bdROLS(2),gamma,CCP_recover_sophi);
    % Minimum distance estimation using final three periods of data with additional exclusion restrictions. 
    objFun = @(bd)EstMD_LastPeriods_restriction(nPeriods,supportX,capPi2,exp(bd(1))/(1+exp(bd(1))),exp(bd(2))/(1+exp(bd(2))),gamma,CCP_recover_sophi); % bd(1) is beta, bd(2) is delta
    Options = optimset('Display','off','Algorithm','interior-point','AlwaysHonorConstraints','bounds',...
                                       'GradObj','off','TolFun',1E-10,'TolX',1E-10,'DerivativeCheck','off');
    [EstimatesMoment_rn(:,j),fval,exitflag_ROLS_rn(:,j)] = fminsearch(objFun,initialv,Options);
    EstimatesMoment_rn(:,j) = exp(EstimatesMoment_rn(:,j))./(1+exp(EstimatesMoment_rn(:,j)));
    bdROLS_rn = EstimatesMoment_rn(:,j);
    [~,~,SSres_ROLS_rn(:,j)] = EstMD_LastPeriods_restriction(nPeriods,supportX,capPi2,bdROLS_rn(1),bdROLS_rn(2),gamma,CCP_recover_sophi);
    
    %% Minimum distance estimation for fully naive agent case.
    % Note that theta is used in the MD estimation for the naive agent case only for
    % verification, but not in the actual estimation.
    objFun_naive_moment = @(bd)verification_naive(exp(bd(1))/(1+exp(bd(1))),exp(bd(2))/(1+exp(bd(2))),supportX,capPi2,gamma,[u1;u2],CCP_recover_naive,CCP_total_tilde,V_total_naive,@flowpayoffs);
    Options = optimset('Display','off','Algorithm','interior-point','AlwaysHonorConstraints','bounds',...
                                       'GradObj','off','TolFun',1E-8,'TolX',1E-10,'DerivativeCheck','off');
    [EstimatesMoment_naive(:,j),fv,exitflag_MD_naive(:,j)] = fminsearch(objFun_naive_moment,initialv,Options);
    EstimatesMoment_naive(:,j) = exp(EstimatesMoment_naive(:,j))./(1+exp(EstimatesMoment_naive(:,j)));
    bdMD_naive = EstimatesMoment_naive(:,j);
    [SSres_MD_naive(:,j), RsqrMD_naive(:,j)] = verification_naive(bdMD_naive(1),bdMD_naive(2),supportX,capPi2,gamma,[u1;u2],CCP_recover_naive,CCP_total_tilde,V_total_naive,@flowpayoffs);
    
    %% Maximum Likelihood Estimation     
    if N_size == 1 || N_size == 7 || N_size == 5 || N_size== 6
        theta = [u1;u2]';
        % eDsM = exponential data and sophisticated agent model
        [Est_eDeM_rn(:,j),Likhood_eDeM_rn(:,j),exitflag_eDeM_rn(:,j)] = MLEstimation(choices_expon,iX_expon,supportX,capPi2,gamma,1,theta,theta_rn,...
                @nLogLik,@flowpayoffs,@CCP_generator,@backInduction,1,delta);
        [Est_eDsM_rn(:,j),Likhood_eDsM_rn(:,j),exitflag_eDsM_rn(:,j)] = MLEstimation(choices_expon,iX_expon,supportX,capPi2,gamma,2,theta,theta_rn,...
                @nLogLik,@flowpayoffs,@CCP_generator,@backInduction,1,delta);
        [Est_eDnM_rn(:,j),Likhood_eDnM_rn(:,j),exitflag_eDnM_rn(:,j)] = MLEstimation(choices_expon,iX_expon,supportX,capPi2,gamma,3,theta,theta_rn,...
                @nLogLik,@flowpayoffs,@CCP_generator,@backInduction,1,delta);
        [Est_sDsM(:,j),Likhood_sDsM(:,j),exitflag_sDsM(:,j)] = MLEstimation(choices_sophi,iX_sophi,supportX,capPi2,gamma,2,theta,thetaY,...
                @nLogLik,@flowpayoffs,@CCP_generator,@backInduction,beta,delta);
        [Est_sDsM_normalization_terminating(:,j),Likhood_sDsM_norm_term(:,j),exitflag_sDsM_norm_term(:,j)] = MLEstimation(choices_sophi,iX_sophi,supportX,capPi2,gamma,2,theta,thetaNormalizationT,...
                @nLogLik,@flowpayoffs,@CCP_generator,@backInduction,beta,delta);
        [Est_sDsM_normalization(:,j),Likhood_sDsM_norm(:,j),exitflag_sDsM_norm(:,j)] = MLEstimation(choices_sophi,iX_sophi,supportX,capPi2,gamma,2,theta,thetaNormalization,...
                @nLogLik,@flowpayoffs,@CCP_generator,@backInduction,beta,delta);
        % Adding exclusion restrictions to estimat.
        [Est_sDsM_rn(:,j),Likhood_sDsM_rn(:,j),exitflag_sDsM_rn(:,j)] = MLEstimation(choices_sophi,iX_sophi,supportX,capPi2,gamma,2,theta,theta_rn,...
                @nLogLik,@flowpayoffs,@CCP_generator,@backInduction,beta,delta);
        [Est_sDsM_rn_norm(:,j),Likhood_sDsM_rn_norm(:,j),exitflag_sDsM_rn_norm(:,j)] = MLEstimation(choices_sophi,iX_sophi,supportX,capPi2,gamma,2,theta,theta_rn_norm,...
                @nLogLik,@flowpayoffs,@CCP_generator,@backInduction,beta,delta);
        [Est_sDsM_rn_norm1(:,j),Likhood_sDsM_rn_norm1(:,j),exitflag_sDsM_rn_norm1(:,j)] = MLEstimation(choices_sophi,iX_sophi,supportX,capPi2,gamma,2,theta,theta_rn_norm1,...
                @nLogLik,@flowpayoffs,@CCP_generator,@backInduction,beta,delta);
        % Naive agent ML-estimation
        [Est_nDnM(:,j),Likhood_nDnM(:,j),exitflag_nDnM(:,j)] = MLEstimation(choices_naive,iX_naive,supportX,capPi2,gamma,3,theta,thetaY,...
                @nLogLik,@flowpayoffs,@CCP_generator,@backInduction,beta,delta);
    end % end if-condition
end % end parfor-loop over simulation runs.

clear I_eDeM I_eDsM I_sDeM I_sDsM
clear MLE_eDeM MLE_eDsM MLE_sDeM MLE_sDsM

% Store counterfactual results.
if nFirms == nFirms_alternatives(1)
    Est_counter = Est_sDsM;
    Est_counter_norm = Est_sDsM_normalization;
    Est_counter_norm_terminating = Est_sDsM_normalization_terminating;
end

Est_sophiMD = mean(EstimatesMoment,2);
Est_sophiMD_rn = mean(EstimatesMoment_rn,2);
Est_naiveMD = mean(EstimatesMoment_naive,2);
Est_sophiOLS_L3 = mean(EstimatesOLS_L3,2);
SD_sophiMD = std(EstimatesMoment,0,2);
SD_sophiMD_rn = std(EstimatesMoment_rn,0,2);
SD_sophiOLS_L3 = std(EstimatesOLS_L3,0,2);
SD_naiveMD = std(EstimatesMoment_naive,0,2);
% Compute mean and std deviation.
para_eDeM_rn = mean(Est_eDeM_rn,2);
para_eDsM_rn = mean(Est_eDsM_rn,2);
para_eDnM_rn = mean(Est_eDnM_rn,2);
para_sDsM = mean(Est_sDsM,2);
para_sDsM_normalization = mean(Est_sDsM_normalization,2);
para_sDsM_normalization_terminating = mean(Est_sDsM_normalization_terminating,2);
para_sDsM_rn = mean(Est_sDsM_rn,2);
para_sDsM_rn_norm = mean(Est_sDsM_rn_norm,2);
para_sDsM_rn_norm1 = mean(Est_sDsM_rn_norm1,2);
para_nDnM = mean(Est_nDnM,2);

SD_eDeM_rn = std(Est_eDeM_rn,0,2);
SD_eDsM_rn = std(Est_eDsM_rn,0,2);
SD_eDnM_rn = std(Est_eDnM_rn,0,2);
SD_sDsM = std(Est_sDsM,0,2);
SD_sDsM_normalization = std(Est_sDsM_normalization,0,2);
SD_sDsM_normalization_terminating = std(Est_sDsM_normalization_terminating,0,2);
SD_sDsM_rn = std(Est_sDsM_rn,0,2);
SD_sDsM_rn_norm = std(Est_sDsM_rn_norm,0,2);
SD_sDsM_rn_norm1 = std(Est_sDsM_rn_norm1,0,2);
SD_nDnM = std(Est_nDnM,0,2);
averageTime = toc;
disp(['nFirms = ', num2str(nFirms)]);
disp(['total time = ', num2str(toc/60),'m']);

% Save the parameters estimated to file
bdtrue = [beta;delta];
if exist(project_paths('OUT_TABLES','bdEst.mat'),'file')~=2
    save(project_paths('OUT_TABLES','bdEst.mat'),"bdtrue")
end
if N_size == 3
    Est100000 = Est_sophiMD;
    SD100000 = SD_sophiMD;
    Est100000_MD_rn = Est_sophiMD_rn;
    SD100000_MD_rn = SD_sophiMD_rn;
    Est100000_naive = Est_naiveMD;
    SD100000_naive = SD_naiveMD;
    SSres100000_OLS = mean(SSres_OLS);
    SSres100000_ROLS = mean(SSres_ROLS);
    SSres100000_ROLS_rn = mean(SSres_ROLS_rn);
    SSres100000_MD_naive = mean(SSres_MD_naive);
    Est100000_OLS_L3 = Est_sophiOLS_L3;
    SD100000_OLS_L3 = SD_sophiOLS_L3;
elseif N_size == 1
    Est10000 = Est_sophiMD;
    SD10000 = SD_sophiMD;
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
    Est10000_mle_rn = para_sDsM_rn;
    SD10000_mle_rn = SD_sDsM_rn;
    Lik10000_mle = -mean(Likhood_sDsM);
    Lik10000_mle_rn = -mean(Likhood_sDsM_rn);
    Lik10000_mle_rn_eDeM = -mean(Likhood_eDeM_rn);
    Lik10000_mle_rn_eDsM = -mean(Likhood_eDsM_rn);
    Lik10000_mle_rn_eDnM = -mean(Likhood_eDnM_rn);
    Lik10000_mle_sDsM_norm_term = -mean(Likhood_sDsM_norm_term);
    Lik10000_mle_sDsM_norm = -mean(Likhood_sDsM_norm);
    Lik10000_mle_sDsM_rn_norm = -mean(Likhood_sDsM_rn_norm);
    Lik10000_mle_sDsM_rn_norm1 = -mean(Likhood_sDsM_rn_norm1);
    Lik10000_mle_nDnM = -mean(Likhood_nDnM);     
    Est10000_mle_rn_norm = para_sDsM_rn_norm;
    SD10000_mle_rn_norm = SD_sDsM_rn_norm;
    Est10000_mle_rn_norm1 = para_sDsM_rn_norm1;
    SD10000_mle_rn_norm1 = SD_sDsM_rn_norm1;
    Est10000_OLS_L3 = Est_sophiOLS_L3;
    SD10000_OLS_L3 = SD_sophiOLS_L3;
    Est10000_mle_rn_eDeM = [1;para_eDeM_rn];
    SD10000_mle_rn_eDeM = [0;SD_eDeM_rn];
    Est10000_mle_rn_eDsM = para_eDsM_rn;
    SD10000_mle_rn_eDsM = SD_eDsM_rn;
elseif N_size == 2
    Est50000 = Est_sophiMD;
    SD50000 = SD_sophiMD;
    Est50000_MD_rn = Est_sophiMD_rn;
    SD50000_MD_rn = SD_sophiMD_rn;
    Est50000_OLS_L3 = Est_sophiOLS_L3;
    SD50000_OLS_L3 = SD_sophiOLS_L3;
    Est50000_naive = Est_naiveMD;
    SD50000_naive = SD_naiveMD;
    SSres50000_OLS = mean(SSres_OLS);
    SSres50000_ROLS = mean(SSres_ROLS);
    SSres50000_ROLS_rn = mean(SSres_ROLS_rn);
    SSres50000_MD_naive = mean(SSres_MD_naive);
elseif N_size == 4
    Est1000000 = Est_sophiMD;
    SD1000000 = SD_sophiMD;
    Est1000000_MD_rn = Est_sophiMD_rn;
    SD1000000_MD_rn = SD_sophiMD_rn;
    Est1000000_naive = Est_naiveMD;
    SD1000000_naive = SD_naiveMD;
    SSres1000000_OLS = mean(SSres_OLS);
    SSres1000000_ROLS = mean(SSres_ROLS);
    SSres1000000_ROLS_rn = mean(SSres_ROLS_rn);
    SSres1000000_MD_naive = mean(SSres_MD_naive);
    Est1000000_OLS_L3 = Est_sophiOLS_L3;
    SD1000000_OLS_L3 = SD_sophiOLS_L3;
elseif N_size == 6
    Est5000 = Est_sophiMD;
    SD5000 = SD_sophiMD;
    Est5000_naive = Est_naiveMD;
    SD5000_naive = SD_naiveMD;
    Est5000_mle = para_sDsM;
    SD5000_mle = SD_sDsM;
    Lik5000_mle = -mean(Likhood_sDsM);
    Est5000_mle_norm = para_sDsM_normalization;
    SD5000_mle_norm = SD_sDsM_normalization;
    Est5000_mle_norm_terminating = para_sDsM_normalization_terminating;
    SD5000_mle_norm_terminating = SD_sDsM_normalization_terminating;
    Est5000_mle_naive = para_nDnM;
    SD5000_mle_naive = SD_nDnM;
    Est5000_mle_rn = para_sDsM_rn;
    SD5000_mle_rn = SD_sDsM_rn;
    Lik5000_mle_rn = -mean(Likhood_sDsM_rn);
    Est5000_mle_rn_norm = para_sDsM_rn_norm;
    SD5000_mle_rn_norm = SD_sDsM_rn_norm;
    Est5000_mle_rn_norm1 = para_sDsM_rn_norm1;
    SD5000_mle_rn_norm1 = SD_sDsM_rn_norm1;
    Est5000_OLS_L3 = Est_sophiOLS_L3;
    SD5000_OLS_L3 = SD_sophiOLS_L3;
    Est5000_mle_rn_eDeM = [1;para_eDeM_rn];
    SD5000_mle_rn_eDeM = [0;SD_eDeM_rn];
    Est5000_mle_rn_eDsM = para_eDsM_rn;
    SD5000_mle_rn_eDsM = SD_eDsM_rn;
elseif N_size == 5
   Est20000 = Est_sophiMD;
    SD20000 = SD_sophiMD;
    Est20000_mle = para_sDsM;
    SD20000_mle = SD_sDsM;
    Est20000_mle_norm = para_sDsM_normalization;
    SD20000_mle_norm = SD_sDsM_normalization;
    Est20000_mle_norm_terminating = para_sDsM_normalization_terminating;
    SD20000_mle_norm_terminating = SD_sDsM_normalization_terminating;
    Est20000_mle_naive = para_nDnM;
    SD20000_mle_naive = SD_nDnM;
    Est20000_mle_rn = para_sDsM_rn;
    SD20000_mle_rn = SD_sDsM_rn;
    Lik20000_mle = -mean(Likhood_sDsM);
    Lik20000_mle_rn = -mean(Likhood_sDsM_rn);
    Lik20000_mle_rn_eDeM = -mean(Likhood_eDeM_rn);
    Lik20000_mle_rn_eDsM = -mean(Likhood_eDsM_rn);
    Lik20000_mle_rn_eDnM = -mean(Likhood_eDnM_rn);
    Lik20000_mle_sDsM_norm_term = -mean(Likhood_sDsM_norm_term);
    Lik20000_mle_sDsM_norm = -mean(Likhood_sDsM_norm);
    Lik20000_mle_sDsM_rn_norm = -mean(Likhood_sDsM_rn_norm);
    Lik20000_mle_sDsM_rn_norm1 = -mean(Likhood_sDsM_rn_norm1);
    Lik20000_mle_nDnM = -mean(Likhood_nDnM); 
    Est20000_mle_rn_norm = para_sDsM_rn_norm;
    SD20000_mle_rn_norm = SD_sDsM_rn_norm;
    Est20000_mle_rn_norm1 = para_sDsM_rn_norm1;
    SD20000_mle_rn_norm1 = SD_sDsM_rn_norm1;
    Est20000_mle_rn_eDeM = [1;para_eDeM_rn];
    SD20000_mle_rn_eDeM = [0;SD_eDeM_rn];
    Est20000_mle_rn_eDsM = para_eDsM_rn;
    SD20000_mle_rn_eDsM = SD_eDsM_rn;
    Est20000_mle_rn_eDnM = para_eDnM_rn;
    SD20000_mle_rn_eDnM = SD_eDnM_rn;
elseif N_size==7
    Est30000_mle = para_sDsM;
    SD30000_mle = SD_sDsM;
    Est30000_mle_norm = para_sDsM_normalization;
    SD30000_mle_norm = SD_sDsM_normalization;
    Est30000_mle_norm_terminating = para_sDsM_normalization_terminating;
    SD30000_mle_norm_terminating = SD_sDsM_normalization_terminating;
    Est30000_mle_naive = para_nDnM;
    SD30000_mle_naive = SD_nDnM;
    Est30000_mle_rn = para_sDsM_rn;
    SD30000_mle_rn = SD_sDsM_rn;
    Lik30000_mle = -mean(Likhood_sDsM);
    Lik30000_mle_rn = -mean(Likhood_sDsM_rn);
    Lik30000_mle_rn_eDeM = -mean(Likhood_eDeM_rn);
    Lik30000_mle_rn_eDsM = -mean(Likhood_eDsM_rn);
    Lik30000_mle_rn_eDnM = -mean(Likhood_eDnM_rn);
    Lik30000_mle_sDsM_norm_term = -mean(Likhood_sDsM_norm_term);
    Lik30000_mle_sDsM_norm = -mean(Likhood_sDsM_norm);
    Lik30000_mle_sDsM_rn_norm = -mean(Likhood_sDsM_rn_norm);
    Lik30000_mle_sDsM_rn_norm1 = -mean(Likhood_sDsM_rn_norm1);
    Lik30000_mle_nDnM = -mean(Likhood_nDnM); 
    Est30000_mle_rn_norm = para_sDsM_rn_norm;
    SD30000_mle_rn_norm = SD_sDsM_rn_norm;
    Est30000_mle_rn_norm1 = para_sDsM_rn_norm1;
    SD30000_mle_rn_norm1 = SD_sDsM_rn_norm1;
    Est30000_mle_rn_eDeM = [1;para_eDeM_rn];
    SD30000_mle_rn_eDeM = [0;SD_eDeM_rn];
    Est30000_mle_rn_eDsM = para_eDsM_rn;
    SD30000_mle_rn_eDsM = SD_eDsM_rn;
end % end if-condition over sample sizes.
end % end loop over sample sizes.

% SW: Continue here.

%% Counterfactual analysis
% Compare the case "tax on not adopting" with and without normalization.
% Number of simulated time periods.
nPeriods_counter = 10;
% Full CCP with true parameters.
theta_true = [u1;u2];
beta_true = beta;
delta_true = delta;
CCP_true_sophi_longer = CCP_output(theta_true,beta_true,beta_true,delta_true,nPeriods_counter,...
                                    gamma,supportX,capPi2,@backInduction,@CCP_generator,@flowpayoffs);
% Full CCP with counterfactual parameters.
theta_true_counter = [u1;u2] - [0;0;0;0;0.5;0.5;0.5;0.5]; 
beta_true_counter = beta;
delta_true_counter = delta;
CCP_true_counter = CCP_output(theta_true_counter,beta_true_counter,beta_true_counter,delta_true_counter,nPeriods_counter,...
                                    gamma,supportX,capPi2,@backInduction,@CCP_generator,@flowpayoffs);
CCP_true_counter_diff = CCP_true_counter - CCP_true_sophi_longer;
% Full CCP with estimated parameters.
theta_est = Est10000_mle(3:end);
beta_est = Est10000_mle(1);
delta_est = Est10000_mle(2);
CCP_total_sophi_longer = CCP_output(theta_est,beta_est,beta_est,delta_est,nPeriods_counter,...
                                    gamma,supportX,capPi2,@backInduction,@CCP_generator,@flowpayoffs);
% Full CCP with counterfactual parameters.
theta_counter = Est10000_mle(3:end) - [0;0;0;0;0.5;0.5;0.5;0.5]; 
beta_counter = Est10000_mle(1);
delta_counter = Est10000_mle(2);
CCP_counter = CCP_output(theta_counter,beta_counter,beta_counter,delta_counter,nPeriods_counter,...
                                    gamma,supportX,capPi2,@backInduction,@CCP_generator,@flowpayoffs);
CCP_counter_diff = CCP_counter - CCP_total_sophi_longer;
% Full CCP with normalization on terminating action.
theta_normTerm = [0;0;0;0;Est10000_mle_norm_terminating(3:end)]; 
beta_normTerm = Est10000_mle_norm_terminating(1);
delta_normTerm = Est10000_mle_norm_terminating(2);
CCP_normTerm = CCP_output(theta_normTerm,beta_normTerm,beta_normTerm,delta_normTerm,nPeriods_counter,...
                                    gamma,supportX,capPi2,@backInduction,@CCP_generator,@flowpayoffs);
% Full CCP with normalization on terminating action and counterfactual.
theta_counter_normTerm = [0;0;0;0;Est10000_mle_norm_terminating(3:end)] - [0;0;0;0;0.5;0.5;0.5;0.5]; 
beta_counter_normTerm = Est10000_mle_norm_terminating(1);
delta_counter_normTerm = Est10000_mle_norm_terminating(2);
CCP_counter_normTerm = CCP_output(theta_counter_normTerm,beta_counter_normTerm,beta_counter_normTerm,delta_counter_normTerm,nPeriods_counter,...
                                    gamma,supportX,capPi2,@backInduction,@CCP_generator,@flowpayoffs);
CCP_counter_normTerm_diff = CCP_counter_normTerm - CCP_normTerm;                                
% Full CCP with normalization on non-terminating action.                                      
theta_norm = [Est10000_mle_norm(3:end);0;0;0;0];
beta_norm = Est10000_mle_norm(1);
delta_norm = Est10000_mle_norm(2);
CCP_norm = CCP_output(theta_norm,beta_norm,beta_norm,delta_norm,nPeriods_counter,...
                                    gamma,supportX,capPi2,@backInduction,@CCP_generator,@flowpayoffs);
% Full CCP with normalization on non-terminating action and counterfactual (in other direction).
theta_counter_norm = [Est10000_mle_norm(3:end);0;0;0;0] + [0.5;0.5;0.5;0.5;0;0;0;0];
beta_counter_norm = Est10000_mle_norm(1);
delta_counter_norm = Est10000_mle_norm(2);
CCP_counter_norm = CCP_output(theta_counter_norm,beta_counter_norm,beta_counter_norm,delta_counter_norm,nPeriods_counter,...
                                    gamma,supportX,capPi2,@backInduction,@CCP_generator,@flowpayoffs);
CCP_counter_norm_diff = CCP_counter_norm - CCP_norm;
% Full estimation with counterfactual on terminating action.
theta_counter = Est10000_mle(3:end) + [0.5;0.5;0.5;0.5;0;0;0;0]; 
beta_counter = Est10000_mle(1);
delta_counter = Est10000_mle(2);
[u1_counter,u2_counter] = flowpayoffs(supportX,theta_counter);
W_terminal_non = [u1_counter,u2_counter];
CCP_terminal_non = CCP_generator(W_terminal_non);    
V_terminal_non = W_terminal_non(:,1) + gamma-log(CCP_terminal_non(:,1)); % action 1 is the terminating action      
[CCP_counter,~,~,~] = backInduction(CCP_terminal_non,V_terminal_non,... 
                             u1_counter,u2_counter,gamma,capPi2,beta_counter,beta_counter,delta_counter,nPeriods_counter);
legendc = {'true model','cf w/o norm','cf w/ norm'};

%% Export data for creating graph in Python.
CCP_counter_diff_data = [CCP_true_counter_diff,CCP_counter_diff,CCP_counter_normTerm_diff];
save(project_paths('OUT_ANALYSIS','CCPCounterdiffData.mat'),'CCP_counter_diff_data');
legendc = {'true model','cf w/o norm','cf w/ norm'};

% Create confidence interval (95%); here for: adoption decision.                       
CCP_sophi_counter_CI = zeros(nPeriods_counter*nSuppX,maxj);
CCP_sophi_counter_CI_u2 = zeros(nPeriods_counter*nSuppX,maxj);
CCP_sophi_counter_norm_terminating_CI = zeros(nPeriods_counter*nSuppX,maxj);
CCP_sophi_counter_norm_CI = zeros(nPeriods_counter*nSuppX,maxj);
% Loop over states.
for i = 1:maxj
    % Counterfactual for the correct model (change utility on non-terminating action).
    theta_counter_i = Est_counter(3:end,i) - [0;0;0;0;0.5;0.5;0.5;0.5];
    beta_counter_i = Est_counter(1);
    delta_counter_i = Est_counter(2);
    [u1_counter_i,u2_counter_i] = flowpayoffs(supportX,theta_counter_i);
    W_terminal_non_i = [u1_counter_i,u2_counter_i];
    CCP_terminal_non_i = CCP_generator(W_terminal_non_i);    
    V_terminal_non_i = W_terminal_non_i(:,1) + gamma-log(CCP_terminal_non_i(:,1)); % action 1 is the terminating action      
    [CCP_sophi_counter,~,~,~] = backInduction(CCP_terminal_non_i,V_terminal_non_i,... 
                                 u1_counter_i,u2_counter_i,gamma,capPi2,beta_counter_i,beta_counter_i,delta_counter_i,nPeriods_counter);
    CCP_sophi_counter_CI(:,i) = CCP_sophi_counter(:,1); % just choose action 1 to execute cf
    % Counterfactual for the model with normalization on terminating action.
    % In this case, the counterfactual is done on the utility from non-terminating action.
    theta_counter_norm_terminating_i = [0;0;0;0;Est_counter_norm_terminating(3:end,i)] - [0;0;0;0;0.5;0.5;0.5;0.5]; % full estimation with counterfactual
    beta_counter_norm_terminating_i = Est_counter_norm_terminating(1);
    delta_counter_norm_terminating_i = Est_counter_norm_terminating(2);
    [u1_counter_norm_terminating_i,u2_counter_norm_terminating_i] = flowpayoffs(supportX,theta_counter_norm_terminating_i);
    W_terminal_non_i = [u1_counter_norm_terminating_i,u2_counter_norm_terminating_i];
    CCP_terminal_non_i = CCP_generator(W_terminal_non_i);    
    V_terminal_non_i = W_terminal_non_i(:,1) + gamma-log(CCP_terminal_non_i(:,1)); % action 1 is the terminating action      
    [CCP_sophi_counter,~,~,~] = backInduction(CCP_terminal_non_i,V_terminal_non_i,... 
                                 u1_counter_norm_terminating_i,u2_counter_norm_terminating_i,gamma,capPi2,beta_counter_norm_terminating_i,beta_counter_norm_terminating_i,delta_counter_norm_terminating_i,nPeriods_counter);
    CCP_sophi_counter_norm_terminating_CI(:,i) = CCP_sophi_counter(:,1); % just choose action 1 to execute cf
    % Counterfactual for the correct model (change utility on terminating action).
    theta_counter_i = Est_counter(3:end,i) + [0.5;0.5;0.5;0.5;0;0;0;0]; % full estimation with counterfactual
    beta_counter_i = Est_counter(1);
    delta_counter_i = Est_counter(2);
    [u1_counter_i,u2_counter_i] = flowpayoffs(supportX,theta_counter_i);
    W_terminal_non_i = [u1_counter_i,u2_counter_i];
    CCP_terminal_non_i = CCP_generator(W_terminal_non_i);    
    V_terminal_non_i = W_terminal_non_i(:,1) + gamma-log(CCP_terminal_non_i(:,1)); % action 1 is the terminating action      
    [CCP_sophi_counter,~,~,~] = backInduction(CCP_terminal_non_i,V_terminal_non_i,... 
                                 u1_counter_i,u2_counter_i,gamma,capPi2,beta_counter_i,beta_counter_i,delta_counter_i,nPeriods_counter);
    CCP_sophi_counter_CI_u2(:,i) = CCP_sophi_counter(:,1); % just choose action 1 to execute cf
    % Counterfactual for the model with normalization on non-terminating action.
    % In this case, the counterfactual is done on the utility from terminating action.
    theta_counter_norm_i = [Est_counter_norm(3:end,i);0;0;0;0] + [0.5;0.5;0.5;0.5;0;0;0;0]; % full estimation with counterfactual
    beta_counter_norm_i = Est_counter_norm(1);
    delta_counter_norm_i = Est_counter_norm(2);
    [u1_counter_norm_i,u2_counter_norm_i] = flowpayoffs(supportX,theta_counter_norm_i);
    W_terminal_non_i = [u1_counter_norm_i,u2_counter_norm_i];
    CCP_terminal_non_i = CCP_generator(W_terminal_non_i);    
    V_terminal_non_i = W_terminal_non_i(:,1) + gamma-log(CCP_terminal_non_i(:,1)); % action 1 is the terminating action      
    [CCP_sophi_counter,~,~,~] = backInduction(CCP_terminal_non_i,V_terminal_non_i,... 
                                 u1_counter_norm_i,u2_counter_norm_i,gamma,capPi2,beta_counter_norm_i,beta_counter_norm_i,delta_counter_norm_i,nPeriods_counter);
    CCP_sophi_counter_norm_CI(:,i) = CCP_sophi_counter(:,1); % just choose action 1 to execute cf
end % end loop over states for counterfactual calculation.
% Extract the confidence band.
% For each band, save the mean, lower bound, and upper bound in first, second, 
% and third columnb, respectively.
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

%% Export the estimation results to tables for formatting in Python.
MD_table_labels = {'true','beta_hat_50000','sd_50000','beta_hat_100000','sd_100000','beta_hat_1000000','sd_1000000'};
estMD_table = table(bdtrue,Est50000,SD50000,Est100000,SD100000,Est1000000,SD1000000, 'RowNames',{'$\beta$','$\delta$'},'VariableNames',MD_table_labels);
newrow = table(nan,SSres50000_ROLS,nan,SSres100000_ROLS,nan,SSres1000000_ROLS,nan,'RowNames',{'Criteria Fun Value'},'VariableNames',MD_table_labels);
estMD_table = [estMD_table;newrow];
writetable(estMD_table, project_paths('OUT_ANALYSIS','md_table_data.csv'),'writeRowNames',true);
estMD_rn_table = table(bdtrue,Est50000_MD_rn,SD50000_MD_rn,Est100000_MD_rn,SD100000_MD_rn,Est1000000_MD_rn,SD1000000_MD_rn, 'RowNames',{'$\beta$','$\delta$'},'VariableNames',MD_table_labels);
newrow = table(nan,SSres50000_ROLS_rn,nan,SSres100000_ROLS_rn,nan,SSres1000000_ROLS_rn,nan,'RowNames',{'Criteria Fun Value'},'VariableNames',MD_table_labels);
estMD_rn_table = [estMD_rn_table;newrow];
writetable(estMD_rn_table, project_paths('OUT_ANALYSIS','md_rn_table_data.csv'),'writeRowNames',true);
estMDnaive_table = table(bdtrue,Est50000_naive,SD50000_naive,Est100000_naive,SD100000_naive,Est1000000_naive,SD1000000_naive, 'RowNames',{'$\beta$','$\delta$'},'VariableNames',MD_table_labels);
newrow = table(nan,SSres50000_MD_naive,nan,SSres100000_MD_naive,nan,SSres1000000_MD_naive,nan,'RowNames',{'Criteria Fun Value'},'VariableNames',MD_table_labels);
estMDnaive_table = [estMDnaive_table;newrow];
writetable(estMDnaive_table, project_paths('OUT_ANALYSIS','md_table_data_naive.csv'),'writeRowNames',true);
OLS_table_labels = {'true','beta_hat_50000','sd_50000','beta_hat_100000','sd_100000','beta_hat_1000000','sd_1000000'};
estOLS_table = table([bdtrue;u1],Est50000_OLS_L3,SD50000_OLS_L3,Est100000_OLS_L3,SD100000_OLS_L3,Est1000000_OLS_L3,SD1000000_OLS_L3, 'RowNames',{'$\beta$','$\delta$','$u_{11}$','$u_{12}$','$u_{13}$','$u_{14}$'},'VariableNames',OLS_table_labels);
newrow = table(nan,SSres50000_OLS,nan,SSres100000_OLS,nan,SSres1000000_OLS,nan,'RowNames',{'Criteria Fun Value'},'VariableNames',OLS_table_labels);
estOLS_table = [estOLS_table;newrow];
writetable(estOLS_table, project_paths('OUT_ANALYSIS','ols_table_data.csv'),'writeRowNames',true);
MLE_table_labels = {'true','beta_hat_10000','sd_10000','beta_hat_20000','sd_20000','beta_hat_30000','sd_30000'};
estMLEsophi_table = table([bdtrue;u1;u2],Est10000_mle,SD10000_mle,Est20000_mle,SD20000_mle,Est30000_mle,SD30000_mle, 'RowNames',{'$\beta$','$\delta$','$u_1(x_{1L}, x_{2L})$','$u_1(x_{1L}, x_{2H})$','$u_1(x_{1H}, x_{2L})$','$u_1(x_{1H}, x_{2H})$','$u_0(x_{1L}, x_{2H})$','$u_0(x_{1L}, x_{2L})$','$u_0(x_{1H}, x_{2H})$','$u_0(x_{1H}, x_{2L})$'},'VariableNames',MLE_table_labels);
newrow = table(nan,Lik10000_mle,nan,Lik20000_mle,nan,Lik30000_mle,nan,'RowNames',{'Log likelihood'},'VariableNames',MLE_table_labels);
estMLEsophi_table = [estMLEsophi_table;newrow];
writetable(estMLEsophi_table, project_paths('OUT_ANALYSIS','mle_table_data.csv'),'writeRowNames',true);
estMLEsophi_norm_table = table([bdtrue;u1],Est10000_mle_norm,SD10000_mle_norm,Est20000_mle_norm,SD20000_mle_norm,Est30000_mle_norm,SD30000_mle_norm, 'RowNames',{'$\beta$','$\delta$','$u_1(x_{1L}, x_{2L})$','$u_1(x_{1L}, x_{2H})$','$u_1(x_{1H}, x_{2L})$','$u_1(x_{1H}, x_{2H})$'},'VariableNames',MLE_table_labels);
newrow = table(nan,Lik10000_mle_sDsM_norm,nan,Lik20000_mle_sDsM_norm,nan,Lik30000_mle_sDsM_norm,nan,'RowNames',{'Log likelihood'},'VariableNames',MLE_table_labels);
estMLEsophi_norm_table = [estMLEsophi_norm_table;newrow];
writetable(estMLEsophi_norm_table, project_paths('OUT_ANALYSIS','mle_norm_table_data.csv'),'writeRowNames',true);
estMLEsophi_normTerm_table = table([bdtrue;u2],Est10000_mle_norm_terminating,SD10000_mle_norm_terminating,Est20000_mle_norm_terminating,SD20000_mle_norm_terminating, ...
    Est30000_mle_norm_terminating,SD30000_mle_norm_terminating, 'RowNames',{'$\beta$','$\delta$','$u_0(x_{1L}, x_{2H})$','$u_0(x_{1L}, x_{2L})$','$u_0(x_{1H}, x_{2H})$','$u_0(x_{1H}, x_{2L})$'},'VariableNames',MLE_table_labels);
newrow = table(nan,Lik10000_mle_sDsM_norm_term,nan,Lik20000_mle_sDsM_norm_term,nan,Lik30000_mle_sDsM_norm_term,nan,'RowNames',{'Log likelihood'},'VariableNames',MLE_table_labels);
estMLEsophi_normTerm_table = [estMLEsophi_normTerm_table;newrow];
writetable(estMLEsophi_normTerm_table, project_paths('OUT_ANALYSIS','mle_normTerm_table_data.csv'),'writeRowNames',true);
estMLEsophi_rn_table = table([bdtrue;u1([1,3]);u2([1,2])],Est10000_mle_rn,SD10000_mle_rn,Est20000_mle_rn,SD20000_mle_rn,Est30000_mle_rn,SD30000_mle_rn, 'RowNames',{'$\beta$','$\delta$','$u_1(x_{1L}, x_{2L})$','$u_1(x_{1H}, x_{2L})$','$u_0(x_{1L}, x_{2H})$','$u_0(x_{1L}, x_{2L})$'},'VariableNames',MLE_table_labels);
newrow = table(nan,Lik10000_mle_rn,nan,Lik20000_mle_rn,nan,Lik30000_mle_rn,nan,'RowNames',{'Log likelihood'},'VariableNames',MLE_table_labels);
estMLEsophi_rn_table = [estMLEsophi_rn_table;newrow];
writetable(estMLEsophi_rn_table, project_paths('OUT_ANALYSIS','mle_rn_table_data.csv'),'writeRowNames',true);
bdtrue_eD = [1,delta]';
MLE_eD_table_labels = {'true','beta_hat_20000','sd_20000','beta_hat_sophi_20000','sd_sophi_20000','beta_hat_naive_20000','sd_naive_20000'};
estMLEsophi_rn_eD_table = table([bdtrue_eD;u1([1,3]);u2([1,2])],...
    Est20000_mle_rn_eDeM,SD20000_mle_rn_eDeM,...
    Est20000_mle_rn_eDsM,SD20000_mle_rn_eDsM,...
    Est20000_mle_rn_eDnM,SD20000_mle_rn_eDnM,...
    'RowNames',{'$\beta$','$\delta$','$u_1(x_{1L}, x_{2L})$','$u_1(x_{1H}, x_{2L})$','$u_0(x_{1L}, x_{2H})$','$u_0(x_{1L}, x_{2L})$'},'VariableNames',MLE_eD_table_labels);
newrow = table(nan,Lik20000_mle_rn_eDeM,nan,Lik20000_mle_rn_eDsM,nan,Lik20000_mle_rn_eDnM,nan,'RowNames',{'Log likelihood'},'VariableNames',MLE_eD_table_labels);
estMLEsophi_rn_eD_table = [estMLEsophi_rn_eD_table;newrow];
writetable(estMLEsophi_rn_eD_table, project_paths('OUT_ANALYSIS','mle_rn_eD_table_data.csv'),'writeRowNames',true);
estMLEnaive_table = table([bdtrue;u1;u2],Est10000_mle_naive,SD10000_mle_naive,Est20000_mle_naive,SD20000_mle_naive,Est30000_mle_naive,SD30000_mle_naive, 'RowNames',{'$\beta$','$\delta$','$u_1(x_{1L}, x_{2L})$','$u_1(x_{1L}, x_{2H})$','$u_1(x_{1H}, x_{2L})$','$u_1(x_{1H}, x_{2H})$','$u_0(x_{1L}, x_{2H})$','$u_0(x_{1L}, x_{2L})$','$u_0(x_{1H}, x_{2H})$','$u_0(x_{1H}, x_{2L})$'},'VariableNames',MLE_table_labels);
newrow = table(nan,Lik10000_mle_nDnM,nan,Lik20000_mle_nDnM,nan,Lik30000_mle_nDnM,nan,'RowNames',{'Log likelihood'},'VariableNames',MLE_table_labels);
estMLEnaive_table = [estMLEnaive_table;newrow];
writetable(estMLEnaive_table, project_paths('OUT_ANALYSIS','mle_naive_table_data.csv'),'writeRowNames',true);