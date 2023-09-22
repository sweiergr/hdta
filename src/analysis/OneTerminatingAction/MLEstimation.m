function [MLEst,Likelihood,exitflag_mle] = MLEstimation(choices,iX,supportX,nA,capPi2,gamma,model,theta,thetaY,nlconstraint,nLogLik,flowpayoffs,CCP_generator,backInduction)

    % This file implements the MLE estimation.
    
    % argument "model" indicates which model is estimated:
    % model == 1: estimate the consistent agent model (beta = betatilde=1)
    % model == 2: estimate the sophisticated agent (beta = betatilde <1)
    % model == 3: estimate the naive agent model (beta <betatilde = 1)
    
    % thetaY == 1 indicates theta is estimated
    
    if thetaY == 1
        ntheta = length(theta); % ntheta is the number of theta to be estimated
    elseif thetaY == 0
        ntheta = 0; % this means theta is set to the true value
    elseif thetaY ==2 || thetaY ==3
        ntheta = 2; % this is the case when normalized theta is estimated
        % that is, the utility of one action is normalized to zero
    end
    
    if model ==1 % estimate the consistent agent model
        if thetaY ==0
            objFunc = @(parameters)nLogLik(choices,iX,supportX,nA,capPi2,theta,...
                        exp(parameters(ntheta+1))/(1+exp(parameters(ntheta+1))),...
                        1,1,gamma,flowpayoffs,CCP_generator,backInduction);  
            Initial = [log(0.79/(1-0.79))];        
        elseif thetaY ==1
            objFunc = @(parameters)nLogLik(choices,iX,supportX,nA,capPi2,parameters(2:ntheta+1),...
                        exp(parameters(1))/(1+exp(parameters(1))),...
                        1,1,gamma,flowpayoffs,CCP_generator,backInduction);  
            Initial = [log(0.79/(1-0.79));theta'];       
        end
    elseif model ==2 % estimate the sophisticated agent model
        if thetaY ==0
            objFunc = @(parameters)nLogLik(choices,iX,supportX,nA,capPi2,theta,...
                        exp(parameters(ntheta+2))/(1+exp(parameters(ntheta+2))),exp(parameters(ntheta+1))/(1+exp(parameters(ntheta+1))),exp(parameters(ntheta+1))/(1+exp(parameters(ntheta+1))),gamma,flowpayoffs,CCP_generator,backInduction);        
            Initial = [log(0.39/(1-0.39)); log(0.79/(1-0.79))];
        elseif thetaY ==1
            objFunc = @(parameters)nLogLik(choices,iX,supportX,nA,capPi2,parameters(3:2+ntheta),...
                        exp(parameters(2))/(1+exp(parameters(2))),exp(parameters(1))/(1+exp(parameters(1))),exp(parameters(1))/(1+exp(parameters(1))),gamma,flowpayoffs,CCP_generator,backInduction);        
            Initial = [log(0.39/(1-0.39)); log(0.79/(1-0.79));theta'];
        elseif thetaY ==2 % estimation with normalization in terminating action case
            objFunc = @(parameters)nLogLik(choices,iX,supportX,nA,capPi2,[0;0;parameters(3:2+ntheta)],...
                        exp(parameters(2))/(1+exp(parameters(2))),exp(parameters(1))/(1+exp(parameters(1))),exp(parameters(1))/(1+exp(parameters(1))),gamma,flowpayoffs,CCP_generator,backInduction);        
            Initial = [log(0.39/(1-0.39)); log(0.79/(1-0.79));theta(3:4)'];
        elseif thetaY == 3 % estimation with another normalization case
            objFunc = @(parameters)nLogLik(choices,iX,supportX,nA,capPi2,[parameters(3:2+ntheta);0;0],...
                        exp(parameters(2))/(1+exp(parameters(2))),exp(parameters(1))/(1+exp(parameters(1))),exp(parameters(1))/(1+exp(parameters(1))),gamma,flowpayoffs,CCP_generator,backInduction);        
            Initial = [log(0.39/(1-0.39)); log(0.79/(1-0.79));theta(1:2)'];
        end
    elseif model == 3 % estimate the complete naive agent model
        if thetaY ==0
            objFunc = @(parameters)nLogLik(choices,iX,supportX,nA,capPi2,theta,...
                        exp(parameters(ntheta+2))/(1+exp(parameters(ntheta+2))),...
                        exp(parameters(ntheta+1))/(1+exp(parameters(ntheta+1))),...
                        1,gamma,flowpayoffs,CCP_generator,backInduction);  
            Initial = [log(0.39/(1-0.39)); log(0.79/(1-0.79))];
        elseif thetaY ==1
            objFunc = @(parameters)nLogLik(choices,iX,supportX,nA,capPi2,parameters(3:2+ntheta),...
                        exp(parameters(2))/(1+exp(parameters(2))),...
                        exp(parameters(1))/(1+exp(parameters(1))),...
                        1,gamma,flowpayoffs,CCP_generator,backInduction);  
            Initial = [log(0.39/(1-0.39)); log(0.79/(1-0.79));theta'];
        end
    end
    OptimizerOptions = optimset('Display','off','Algorithm','interior-point','AlwaysHonorConstraints','bounds',...
                                   'GradObj','off','TolFun',1E-10,'TolX',1E-10,'DerivativeCheck','off');
    [MLEst,Likelihood,exitflag_mle]=fminsearch(objFunc,Initial,OptimizerOptions);
    if model ==2 || model == 3
        MLEst(1:2) = exp(MLEst(1:2))./(1+exp(MLEst(1:2)));  
    elseif model == 1
        MLEst(1) = exp(MLEst(1))./(1+exp(MLEst(1)));
    end
end