%% Read Data
clc
clear
% Table for Bellman computation
Bellman = readtable("data_for_bellman_computations.csv");
% Only take naics_recode value of 1
Bellman = Bellman((Bellman.naics_recode==1),:);
% Use Omega1 to represent the fixed states
Omega1 = Bellman.omega1; % It turns out there are 45 possible Omega1
% Table for Omega
T1 = readtable('value_part2_problem1_thetaBGL.csv');
Omega = T1{:,[1,5:8]};
Omega(:,1)=Omega1;
% Table for Omega tilde
T2 = readtable("valuetilde_part2_problem1_thetaBGL.csv");
Omegat = T2{:,[1,5:8]};
Omegat(:,1)=Omega(:,1);

%% 1) Computation of the plant’s dynamic optimization decision for
% parameters

Coeff = [2,-0.5,-0.5,-5,-0.1]; %X,I,V,F,H
%Coeff = [2.872,-0.049,-0.077,-5.980,-0.065]; %BGL check
[NewV,Vtilde,Investprob]=Bellmanfun(Coeff,Omega,Omegat,Bellman);

disp(Investprob(:,4,:,:,:));

%% 2) Nested fixed point quasi-maximum likelihood estimation
% Read the sample
sample = readtable("analysis_data.csv");
% Only take naics_recode value of 1 and the plants are not comliant
sample = sample((sample.naics_recode==1 & sample.compliance==0),:);
somega = sample(:,{'omega1','DAV','ordered_violator','lag_investment','violation','investment'});

%% The quasi log-likelihood of parameter theta is
[loglike,ll] = LogLike(Investprob,somega);
disp('The log likelihood of theta hat is:');
disp(loglike);

%% 2b Maximize the likelihood
% Starting value is thetaBGL
X0 = [2.872,-0.049,-0.077,-5.980,-0.065];
% Create a function such that Omega,Omegat,Bellman,somega are parameters.
% Matlab just needs to maximize likelihood by choosing theta.
fun = @(x)-loglikefun(x,Omega,Omegat,Bellman,somega);
option = optimset('Display','Iter','MaxIter',1,'PlotFcns',@optimplotfval); 
% Do at most 5 iterations in fminsearch function
[x,fval,exitflag,output] = fminsearch(fun,X0,option);
disp('Theta_ML')
disp(x);

%% 2c Find the standard errors
% First estimate the variance covariance matrix
VarML = varml(x,ll,Omega,Omegat,Bellman,somega);
% Standard errors are just the square root of diagnal of this matrix
disp(VarML);
