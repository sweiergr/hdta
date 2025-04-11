%{
    This function implements the MLE estimation.
    Argument "model" indicates which model is estimated:
    - model == 1: estimate the consistent agent model (beta = betatilde=1)
    - model == 2: estimate the sophisticated agent (beta = betatilde <1)
    - model == 3: estimate the naive agent model (beta <betatilde = 1)
    - thetaY == 1 indicates theta is estimated

    Inputs:
    - choices: choices data
    - iX: data for index of state.
    - supportX: support of state vector.
    - capPi2: transition matrix for state variable
    - gamma: mean of the random shock
    - model: which model is estimated
    - theta: parameters in the flow utility
    - thetaY: whether to estimate parameters in flow utility
    - nLogLik: nLogLik function
    - flowpayoffs: flowpayoffs function
    - CCP_generator: CCP_generator function 
    - backInduction: backInduction function
    - beta: present-bias discount factor 
    - delta: exponential discount factor

    Outputs:
    - MLEst: Estimates of Maximum likelihood
    - Likelihood: likelihood value
    - exitflag_mle: exitflag from MLE estimation.
%}
function [MLEst,Likelihood,exitflag_mle] = MLEstimation(choices,iX,supportX,capPi2,gamma,model,theta,thetaY,...
                nLogLik,flowpayoffs,CCP_generator,backInduction,beta,delta)
 
    disturb=0.95;
    initialv=0.9*[log(0.98*beta/(1-0.98*beta)); log(0.98*delta/(1-0.98*delta))];

    if thetaY == 1
        ntheta = length(theta); % ntheta is the number of theta-parameters to be estimated.
    elseif thetaY == 0
        ntheta = 0; % This indicates theta is set to the true value.
    elseif thetaY ==2 || thetaY ==3
        ntheta = 4; % This is the case when a normalized theta is estimated.
        % That is, the utility at one action is normalized to zero
    end
    % Estimate the time-consistent agent model.
    if model ==1 
        if thetaY ==0
            objFunc = @(parameters)nLogLik(choices,iX,supportX,capPi2,theta,...
                        exp(parameters(ntheta+1))/(1+exp(parameters(ntheta+1))),...
                        1,1,gamma,flowpayoffs,CCP_generator,backInduction);  
            Initial = [initialv(2)];        
        elseif thetaY ==1
            objFunc = @(parameters)nLogLik(choices,iX,supportX,capPi2,parameters(2:ntheta+1),...
                        exp(parameters(1))/(1+exp(parameters(1))),...
                        1,1,gamma,flowpayoffs,CCP_generator,backInduction);  
%
            Initial = [initialv(2);disturb*theta']; 
        % Model with exclusion restriction.
        elseif thetaY ==4 
            objFunc = @(parameters)nLogLik(choices,iX,supportX,capPi2,[parameters(2);parameters(2);parameters(3);parameters(3);parameters(4);parameters(5);parameters(4);parameters(5)],...
                        exp(parameters(1))/(1+exp(parameters(1))),...
                        1,1,gamma,flowpayoffs,CCP_generator,backInduction);     
            Initial = [initialv(2);disturb*theta([1,3,5,6])'];             
        end
    % Estimate the sophisticated agent model.
    elseif model ==2 
        if thetaY ==0
            objFunc = @(parameters)nLogLik(choices,iX,supportX,capPi2,theta,...
                        exp(parameters(ntheta+2))/(1+exp(parameters(ntheta+2))),exp(parameters(ntheta+1))/(1+exp(parameters(ntheta+1))),exp(parameters(ntheta+1))/(1+exp(parameters(ntheta+1))),gamma,flowpayoffs,CCP_generator,backInduction);        
        Initial=initialv;
        elseif thetaY ==1
            objFunc = @(parameters)nLogLik(choices,iX,supportX,capPi2,parameters(3:2+ntheta),...
                        exp(parameters(2))/(1+exp(parameters(2))),exp(parameters(1))/(1+exp(parameters(1))),exp(parameters(1))/(1+exp(parameters(1))),gamma,flowpayoffs,CCP_generator,backInduction);        
            Initial=[initialv;disturb*theta'];
        % Estimation with normalization on terminating action.
        elseif thetaY ==2 
            objFunc = @(parameters)nLogLik(choices,iX,supportX,capPi2,[0;0;0;0;parameters(3:2+ntheta)],...
                        exp(parameters(2))/(1+exp(parameters(2))),exp(parameters(1))/(1+exp(parameters(1))),exp(parameters(1))/(1+exp(parameters(1))),gamma,flowpayoffs,CCP_generator,backInduction);        
               Initial = [initialv; disturb*theta(ntheta+1:ntheta*2)'];
        % Estimation with alternative normalization.
        elseif thetaY == 3 
            objFunc = @(parameters)nLogLik(choices,iX,supportX,capPi2,[parameters(3:2+ntheta);0;0;0;0],...
                        exp(parameters(2))/(1+exp(parameters(2))),exp(parameters(1))/(1+exp(parameters(1))),exp(parameters(1))/(1+exp(parameters(1))),gamma,flowpayoffs,CCP_generator,backInduction);        
            Initial = [initialv;disturb*theta(1:ntheta)'];
        % Estimation with exclusion restriction.
        elseif thetaY ==4 
            objFunc = @(parameters)nLogLik(choices,iX,supportX,capPi2,[parameters(3);parameters(3);parameters(4);parameters(4);parameters(5);parameters(6);parameters(5);parameters(6)],...
                        exp(parameters(2))/(1+exp(parameters(2))),exp(parameters(1))/(1+exp(parameters(1))),exp(parameters(1))/(1+exp(parameters(1))),gamma,flowpayoffs,CCP_generator,backInduction);        
             Initial = [initialv;disturb*theta([1,3,5,6])'];
        % Eestimation with exclusion restriction and normalization on terminating action.
        elseif thetaY ==5 
            objFunc = @(parameters)nLogLik(choices,iX,supportX,capPi2,[0;0;0;0;parameters(3);parameters(4);parameters(3);parameters(4)],...
                        exp(parameters(2))/(1+exp(parameters(2))),exp(parameters(1))/(1+exp(parameters(1))),exp(parameters(1))/(1+exp(parameters(1))),gamma,flowpayoffs,CCP_generator,backInduction);        
             Initial = [initialv;disturb*theta([5,6])'];
        % Estimation with exclusion restriction and normalization on one state value.
        elseif thetaY ==6 
            objFunc = @(parameters)nLogLik(choices,iX,supportX,capPi2,[0;0;parameters(3);parameters(3);parameters(4);parameters(5);parameters(4);parameters(5)],...
                        exp(parameters(2))/(1+exp(parameters(2))),exp(parameters(1))/(1+exp(parameters(1))),exp(parameters(1))/(1+exp(parameters(1))),gamma,flowpayoffs,CCP_generator,backInduction);        
              Initial = [initialv;disturb*theta([3,5,6])'];
        end
    % Estimate the naive agent model.
    elseif model == 3 
        if thetaY ==0
            objFunc = @(parameters)nLogLik(choices,iX,supportX,capPi2,theta,...
                        exp(parameters(ntheta+2))/(1+exp(parameters(ntheta+2))),...
                        exp(parameters(ntheta+1))/(1+exp(parameters(ntheta+1))),...
                        1,gamma,flowpayoffs,CCP_generator,backInduction);  
             Initial = initialv;
        elseif thetaY ==1
            objFunc = @(parameters)nLogLik(choices,iX,supportX,capPi2,parameters(3:2+ntheta),...
                        exp(parameters(2))/(1+exp(parameters(2))),...
                        exp(parameters(1))/(1+exp(parameters(1))),...
                        1,gamma,flowpayoffs,CCP_generator,backInduction);  
            Initial = [initialv;disturb*theta'];
        % Estimation with exclusion restriction.
        elseif thetaY ==4 
            objFunc = @(parameters)nLogLik(choices,iX,supportX,capPi2,[parameters(3);parameters(3);parameters(4);parameters(4);parameters(5);parameters(6);parameters(5);parameters(6)],...
                        exp(parameters(2))/(1+exp(parameters(2))),exp(parameters(1))/(1+exp(parameters(1))),1,gamma,flowpayoffs,CCP_generator,backInduction);        
            Initial = [initialv;disturb*theta([1,3,5,6])'];
        end
    end
   options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'off', ...
        'TolFun',1E-10,'TolX',1E-10,'MaxIterations', 10000,'MaxFunctionEvaluations', 5000);

   [MLEst,Likelihood,exitflag_mle]=fminunc(objFunc,Initial,options);
    if model ==2 || model == 3
        MLEst(1:2) = exp(MLEst(1:2))./(1+exp(MLEst(1:2)));  
    elseif model == 1
        MLEst(1) = exp(MLEst(1))./(1+exp(MLEst(1)));
    end
end