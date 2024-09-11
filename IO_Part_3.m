%% Read Data
% This part is the same and Part 2
clc
clear
%mkdir('C:\Users\xuzik\Downloads\IO Assignment\Assignment\Part 2\');
addpath('C:\Users\xuzik\Downloads\IO Assignment\Assignment\Part 2\');
addpath('investprob_grids\');
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
% Table for transition probability
transitionprob = readtable("transitionprob_part3_problem1_DAVgrid2_thetaBGL.csv");
% Table for steady state probability
steadyprob = readtable("steadystate_part3_problem2_DAVgrid2_thetaBGL.csv");
%% 1) Calculation of Transition Matrix
% 1.1 Find the probability of investment probability
% Coeff = [2,-0.5,-0.5,-5,-0.1]; %X,I,V,F,H
Coeff = [2.872,-0.049,-0.077,-5.980,-0.065]; %BGL
[NewV,Vtilde,Investprob]=Bellmanfun(Coeff,Omega,Omegat,Bellman);

%% 1.3
% Given each potential state in omega, find the transition probability to
% each possible omegatilde
% Transition probability given regulation action
N = size(Omega,1); % The number of states
Tran = zeros(N,80,3);
% Rewrite the state transition probability in terms of omega to omegatilde,
% the 3 dimensions represent initial state omega,
% whether there is a violation,
% status
Tran_o2ot = zeros(N,2,3); % omega, vio, status
for j = 1:80
            % Given a state omega, the probability of regulation j and
            % transit to status 1,2,3 is Tran
            Tran(:,j,1) = Bellman.(15+(j-1)*7).*Bellman.(14+(j-1)*7);
            Tran(:,j,2) = Bellman.(16+(j-1)*7).*Bellman.(14+(j-1)*7);
            Tran(:,j,3) = Bellman.(17+(j-1)*7).*Bellman.(14+(j-1)*7);
            % violations
            if j<=20 || (j>=41 && j<=60) % vio = 0
            Tran_o2ot(:,1,:)=Tran_o2ot(:,1,:)+Tran(:,j,:);
            else % vio = 1
            Tran_o2ot(:,2,:)=Tran_o2ot(:,2,:)+Tran(:,j,:);
            end
end
%% Finally we can construct the transition probability
% Define an auxiliary transition matrix from omega to omegatilde
To2ot = zeros(45,20,3,2,2,...
    2,3); % first 5 dimensions are omega
for i = 1:N % for each state in omega
    To2ot(Omega(i,1),Omega(i,2)*2+1,Omega(i,3)+1,Omega(i,4)+1,Omega(i,5)+1,:,:)...
        =Tran_o2ot(i,:,:);
end

%% 1.3 Find transition probability
% Given omega1 = 1 (gravity = region = 1), DAVgrid=4 (DAV=2),
% violation = 1, ordered violator = 2 (HPV), 
% lag_investment = 1
% davgrid = 2; % test
omega1 = 1;
davgrid = 2;
status = 2; % HPV
laginv = 1;
vio = 1;
statefrom = [omega1,davgrid,status,laginv,vio]; % Initial state
Tables = cell(1,4);
Tables{1} = transitionprob; Tables{2}=Omega;
Tables{3} = Omegat; Tables{4} = To2ot; % Tables that required
tranprobnew = transitionprob; % replicate the reference table
% Then figure out the transition probability for each omegatilde 
% next period.
tran_ot2ot = tranmat(Investprob,Tables); % Transition Matrix
% Extract a transition vector from the matrix
transprob2 = tranvec(statefrom,tran_ot2ot,Tables); % Extract a transition vector
tranprobnew = addvars(tranprobnew,transprob2);
disp('The state transition probability for Omega1=1, DAVgrid=4, Violation=1,Odered Violator=2')
disp('and Lag Investment = 1 is:');
head(tranprobnew{:,["transprob","transprob2"]}); % test


%% 2) Computation of the steady state distribution
pi = steadypi(Investprob,Tables);
%head(pi);
%% 2.2 Report the state probabilities for omega1 = 1 and DAVgrid=4
ssprob2 = zeros(8,1);
M = size(transitionprob,1); % 161 number of states in each fixed state
Davgrid = transitionprob.DAVgridto;
Status = transitionprob.orderedvio_to;
Lag_inv = transitionprob.laginv_to;
Vio = transitionprob.violationto;
for i=1:8
    davgrid = 2;
    status = steadyprob.OrderedVio(i);
    laginv = steadyprob.Lag_Inv(i);
    vio = steadyprob.Violation(i);
    % search the index for this state
    for j=1:M
        if (Davgrid(j) == davgrid) &&...
           (Status(j) == status) &&...
           (Lag_inv(j) == laginv) &&...
           (Vio(j) == vio)

        ssprob2(i) = pi(1,j);
        end
    end
end
steadyprobnew = addvars(steadyprob,ssprob2);
disp('Steady state probability')
disp(steadyprobnew{:,["SSProb","ssprob2"]}); % test

%% 3）Computation of the model moment inputs across parameter grid values
% Read investment probabilities of each parameter grid values
%d = dir('investprob_grids\investprob*.csv'); % get a list of all csv files with investprob prefix
data = cell(1,500); % preallocate a cell array to store the data
M1 = zeros(45,161,500); % First set of moments
M2 = zeros(45,161,500); % Second set of moments
Investprobmat = zeros(45,161);
% Extract attributes from the investment prob
data{1} = readtable('investprob_grids\investprob1.csv');
omega1 = data{1}.omega1;
dav = data{1}.DAV;
status = data{1}.ordered_violator;
lag_inv = data{1}.lag_investment;
vio = data{1}.violation;
% Get the moments for each grid parameters
for i = 1:500
    tic
    filename = ['investprob_grids\investprob',num2str(i),'.csv'];
    data{i} = readtable(filename); % read each file and store it in the cell array
    % Calculate the first set of moments
    % Read the investment probability as a vector
    Investprobvec = data{i}.Investprob;
    % Then translate the vector into a 5-D matrix
    for j=1:N       
        Investprob(omega1(j),dav(j)*2+1,status(j)+1,lag_inv(j)+1,vio(j)+1)...
            =Investprobvec(j);
    end
    % M1 stores the steady state distribution in a 500 by 45 by 161 matrix,
    % 500 parameter grid points, 45 Omega1, and 161 Omega2 in each Omega1
    M1(:,:,i) = steadypi(Investprob,Tables);
    % Calculate the second set of moments
    for k=1:45
        for j=1:161
        % Reshape investment probability into a 45 by 161 matrix
        Investprobmat(k,j) =...
            Investprob(k,Davgrid(j)+1,Status(j)+1,Lag_inv(j)+1,Vio(j)+1);
        end
    end
    M2(:,:,i) = M1(:,:,i).*Investprobmat;
    toc
    disp(i);
end
% Keep only non-compliance state for M2
M2 = M2(:,2:end,:);
%% Test the moments for 2nd parameter vector
Moments = readtable("modelmoments_part3_problem3_param2.csv");
modelmoments2 = zeros(720,1);
% Read the state
    davgrid = Moments.davgrid;
    status = Moments.orderedvio;
    laginv = Moments.laginv;
    vio = Moments.violation;
    omega2 = zeros(720,1);
for i= 1:720
    for j = 1:161
        if (Davgrid(j) == davgrid(i)) &&...
           (Status(j) == status(i)) &&...
           (Lag_inv(j) == laginv(i)) &&...
           (Vio(j) == vio(i))
            %omega2 = 4+(mod(i-1,8))*20;
            omega2(i) = j;
            break
        end
    end
    if Moments.momset(i) == 1
        omega1 = idivide(i,int16(8),"ceil");
        modelmoments2(i) = M1(omega1,omega2(i),2);
    else
        omega1 = idivide(i-360,int16(8),"ceil");
        modelmoments2(i) = M2(omega1,omega2(i),2);
    end
end

Momemntsnew = addvars(Moments,modelmoments2);
head(Momemntsnew{:,["modelmoment1","modelmoments2"]});% test
%% Report the mean and standard deviation
% calculate mean
mean1 = mean(M1(:,:,1:10),1:2); % M1's mean, 45 by 161
mean2 = mean(M2(:,:,1:10),1:2); % M2's mean
% calculate standard deviation
std1 = std(M1(:,:,1:10),1,1:2); % M1's std
std2 = std(M2(:,:,1:10),1,1:2);
disp('mean of M1 is')
disp(mean1);
disp('mean of M2 is')
disp(mean2);
disp('Std of M1 is')
disp(std1);
disp('Std of M2 is')
disp(std2);
%% 4) Computation of the data moments 
% Read the panel data
sample = readtable("analysis_data.csv");
% Only take naics_recode value of 1 and the plants are not comliant
sample = sample((sample.naics_recode==1),:);
somega = sample(:,{'omega1','DAV','ordered_violator','lag_investment','violation','investment'});
%% Find data moments
MD1 = zeros(45,161);
MD2 = zeros(45,161);
DAVh = ceil(somega.DAV*2);
DAVl = floor(somega.DAV*2);
% Given a fixed state omega1
disp('time for calculating data moments')
tic
for i = 1:45
    % find the frequency of each initial variable state
    for j = 1:161
        % First set of data moments
        Nh = (somega.omega1==i & DAVh == Davgrid(j) &...
            somega.ordered_violator == Status(j) &...
            somega.lag_investment == Lag_inv(j) &...
            somega.violation == Vio(j)); % number of states when DAV= DAVh
        Nl = (somega.omega1==i & DAVl == Davgrid(j) &...
            somega.ordered_violator == Status(j) &...
            somega.lag_investment == Lag_inv(j) &...
            somega.violation == Vio(j)); % number of states when DAV= DAVl
        MD1(i,j) = sum(Nl+(somega.DAV*2-DAVl).*(Nh-Nl))...
            /nnz(somega.omega1==i);
        Mh = (somega.omega1==i & DAVh == Davgrid(j) &...
            somega.ordered_violator == Status(j) &...
            somega.lag_investment == Lag_inv(j) &...
            somega.violation == Vio(j) &...
            somega.investment == 1); % number of states when DAV=DAVh
        Ml = (somega.omega1==i & DAVh == Davgrid(j) &...
            somega.ordered_violator == Status(j) &...
            somega.lag_investment == Lag_inv(j) &...
            somega.violation == Vio(j) &...
            somega.investment == 1); % number of states when DAV=DAVl
        MD2(i,j) = sum(Ml+(somega.DAV*2-DAVl).*(Mh-Ml))/nnz(somega.omega1==i);
    end
end
toc
% Keep only 160 rows for MD2
MD2 = MD2(:,2:end);
%% Report the mean and standard deviation of moments
% mean of MD
disp("Mean of data moments are:");
disp(mean(MD1,"all"));
disp(mean(MD2,"all"));
disp("Standard deviation of data moments are:");
disp(mean(MD1,"all"));
disp(mean(MD2,"all"));

%% Verify moments using datamoments_part3_problem4_param2.csv
DataMoments = readtable("datamoments_part3_problem4_param2.csv");
datamoments2 = zeros(720,1);
% Read the state
    davgrid = DataMoments.davgrid;
    status = DataMoments.orderedvio;
    laginv = DataMoments.laginv;
    vio = DataMoments.violation;
    omega2 = zeros(720,1);
for i= 1:720
    for j = 1:161
        if (Davgrid(j) == davgrid(i)) &&...
           (Status(j) == status(i)) &&...
           (Lag_inv(j) == laginv(i)) &&...
           (Vio(j) == vio(i))
            %omega2 = 4+(mod(i-1,8))*20;
            omega2(i) = j;
            break
        end
    end
    if DataMoments.momset(i) == 1
        omega1 = idivide(i,int16(8),"ceil");
        datamoments2(i) = MD1(omega1,omega2(i));
    else
        omega1 = idivide(i-360,int16(8),"ceil");
        datamoments2(i) = MD2(omega1,omega2(i)-1);
    end
end

DataMomemntsnew = addvars(DataMoments,datamoments2);
head(DataMomemntsnew{:,["datamoment","datamoments2"]});% test

%% 5) Calculate parameter estimates
% Create moment_data_class_assignment_use.csv
% First reshape those moments
%Mv1 = reshape(M1,45*161,500);
%Mv2 = reshape(M2,45*160,500);
%MDv1 =reshape(MD1,[],1);
%MDv2 =reshape(MD2,[],1);
%moment_data_class_assignment_use = [MDv1, Mv1; MDv2, Mv2];
%writematrix(moment_data_class_assignment_use,'moment_data_class_assignment_use.csv');

% read matrix
weights = readmatrix("weightparams_class_assignment.csv");
disp(weights);