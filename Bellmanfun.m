function [NewV,Vtilde,Investprob] = Bellmanfun(Coeff,Omega,Omegat,Bellman)
tic
beta = 0.95^0.25;
gamma = 0;
N = size(Omega,1); % The number of states

% initialization
% Value function at the beginning of the period
NewV = zeros([45,20,3,2,2]); % Omega1, DAVgrid, violator status, lag_inv, lag2_inv
OldV = zeros([45,20,3,2,2]);
% Value function right after the regulator has moved.
Vtilde = zeros([45,20,3,2,2]); % Omega1, DAVgrid, violator status, lag_inv, violation 
% Investment probability
Investprob = zeros([45,20,3,2,2]);% Omega1, DAVgrid, violator status, lag_inv, violation
% Regulation probability
P_regu = zeros(N,80);
% Transition probability
Tran = zeros(N,80,3);
% Static utility given a state Omega and a regulation action
U = zeros(N,80,3);
EU = zeros(N,80);
% Vtilde under each regulation action
Vr = zeros(N,2,3);
% Expected Vtilde conditional on inspection, violation and fines
EVr = zeros(N,80);
% The difference between new and old value function
norm = 1; 

for j = 1:80
            % The probability of regulation action j is implemented is
            P_regu(:,j) = Bellman.(14+(j-1)*7);
            % Under regulation action j, the plant's violator status
            % becomes compliance, regular violator and HPV with probability
            % Tran0, Tran1 and Tran2 respectively.
            Tran(:,j,1) = Bellman.(15+(j-1)*7);
            Tran(:,j,2) = Bellman.(16+(j-1)*7);
            Tran(:,j,3) = Bellman.(17+(j-1)*7);
            % Static utility under Omega i and regulation action j and
            % transit to compliance.
            U(:,j,1) = Coeff(2)*Bellman.(11+(j-1)*7)+Coeff(3)*Bellman.(12+(j-1)*7)...
                + Coeff(4)*Bellman.(13+(j-1)*7);
            % transit to regular violator
            U(:,j,2) = Coeff(2)*Bellman.(11+(j-1)*7)+Coeff(3)*Bellman.(12+(j-1)*7)...
                + Coeff(4)*Bellman.(13+(j-1)*7);
            % transit to HPV
            U(:,j,3) = Coeff(2)*Bellman.(11+(j-1)*7)+Coeff(3)*Bellman.(12+(j-1)*7)...
                + Coeff(4)*Bellman.(13+(j-1)*7)+Coeff(5);
            % Expected utility conditional on Omega and inspection,
            % violation and fines
            EU(:,j) = dot(U(:,j,:),Tran(:,j,:),3); 
 end
 % EEU is the expected static utility
 EEU = zeros(N,1);
 for i = 1:N
 EEU(i) = dot(EU(i,:),P_regu(i,:),2);
 end




%% Iteration
iter = 0;
while norm > 1e-6
    iter = iter + 1;
% Loop through states Omegatilde. For each Omegatilde
%% Upper bound and lower bound of OldV(Omega) next period. If
% X=1 in this period, then Omega in next
% period has Omega.lag_inv=1 and Omega.lag2_inv=Omegat.lag_inv
Vinv_u = zeros(N,1);
Vinv_l = zeros(N,1);
Vnotinv_u = zeros(N,1);
Vnotinv_l = zeros(N,1);
omega1 = Omegat(:,1); % Fixed state
DAV = Omegat(:,2); % From 0 to 9.5
status = Omegat(:,3); % Violation status
lag_inv = Omegat(:,4); % Lag 1 quarter investment
vio = Omegat(:,5); % Violation
DAVgrid_new = DAV*0.9+vio; % DAV of next period, may not be an integer
% Linear interpolation on DAV
A = [ceil(DAVgrid_new*2)+1,20*ones(N,1)]; % This matrix is prepared for DAVu
DAVu = min(A,[],2);
DAVl = floor(DAVgrid_new*2)+1;
%% Loop 1
for i=1:N
Vinv_u(i) = OldV(omega1(i),DAVu(i),status(i)+1,2,lag_inv(i)+1); % +1 because matlab counts from 1 instead of 0.
Vinv_l(i) = OldV(omega1(i),DAVl(i),status(i)+1,2,lag_inv(i)+1);
Vnotinv_u(i) = OldV(omega1(i),DAVu(i),status(i)+1,1,lag_inv(i)+1);
Vnotinv_l(i) = OldV(omega1(i),DAVl(i),status(i)+1,1,lag_inv(i)+1);
end

Vinv = Vinv_l + (DAVgrid_new*2+1-DAVl).*(Vinv_u-Vinv_l);
Vnotinv = Vnotinv_l + (DAVgrid_new*2+1-DAVl).*(Vnotinv_u-Vnotinv_l);
NewVtil = log(exp(beta*Vinv-Coeff(1))+exp(beta*Vnotinv))+gamma;
Invp = exp(beta*Vinv-Coeff(1))./(exp(beta*Vinv-Coeff(1))+exp(beta*Vnotinv));

for i = 1:N
% Name the variables/vectors
omega1 = Omegat(i,1); % Fixed state
DAV = Omegat(i,2); % From 0 to 19
status = Omegat(i,3); % Violation status
lag_inv = Omegat(i,4); % Lag 1 quarter investment
vio = Omegat(i,5); % Violation
    if status ==0
            % If the plant is in compliance, then at the beginning of next
            % period, DAV = 0, status is compliance, two lags of investment
            % are also 0.
            NewVtil(i) = beta*OldV(omega1,1,1,1,1) + gamma;
            Vtilde(omega1,DAV*2+1,status+1,lag_inv+1,vio+1) = NewVtil(i);
    else
        % Update Vtilde
        Vtilde(omega1,DAV*2+1,status+1,lag_inv+1,vio+1) = NewVtil(i);
        % Investment probability given the plant in on Omega tilde.
        Investprob(omega1,DAV*2+1,status+1,lag_inv+1,vio+1) = Invp(i);
    end      
end

%% Vtilde in Loop 2
 for i = 1:N
 % Vtilde for each potential violator status and violation
 omega1 = Omega(i,1);
        DAV = Omega(i,2);
        lag_inv = Omega(i,4);
            Vr(i,1,1) = Vtilde(omega1,1,1,1,1);
            Vr(i,2,1) = Vtilde(omega1,1,1,1,1);
            Vr(i,1,2) = Vtilde(omega1,DAV*2+1,2,lag_inv+1,1);
            Vr(i,2,2) = Vtilde(omega1,DAV*2+1,2,lag_inv+1,2);
            Vr(i,1,3) = Vtilde(omega1,DAV*2+1,3,lag_inv+1,1);
            Vr(i,2,3) = Vtilde(omega1,DAV*2+1,3,lag_inv+1,2);   
 end

 % 80 regulation actions in terms of inspection, violation and fines
 for j = 1:80

            % Let vio be an indicator whether there is a violation in
            % regulation action j
         
               % Expected Vtilde conditional on Omega and inspection,
            % violation and fines
            if j<=20 || (j>=41 && j<=60) % vio = 0
            EVr(:,j) = dot(Vr(:,1,:),Tran(:,j,:),3);
            else % vio = 1
            EVr(:,j) = dot(Vr(:,2,:),Tran(:,j,:),3);
            end
 end



%% Loop through states Omega. For each Omega, find V(Omega) using equation
    for i = 1:N        
        % A1
        omega1 = Omega(i,1);
        DAV = Omega(i,2);
        status = Omega(i,3);
        lag_inv = Omega(i,4);
        lag2_inv = Omega(i,5);

        % Then we can calculate V(Omega) as the sum of expected static
        % utility and expected Vtilde conditional on Omega being the
        % initial state
        NewV(omega1,DAV*2+1,status+1,lag_inv+1,lag2_inv+1) =...
            dot(EVr(i,:),P_regu(i,:))+EEU(i);
    end
    
norm = max(abs(NewV - OldV),[],'all');
OldV = NewV; % Update value function

%disp([norm,iter])
end
toc
end

