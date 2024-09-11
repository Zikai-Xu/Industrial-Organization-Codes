function tran_ot2ot = tranmat(Investprob,Tables)
%TRANMAT Summary of this function goes here
%   Detailed explanation goes here
% Read tables
%transitionprob = Tables{1};
Omega = Tables{2};
Omegat = Tables{3};
To2ot = Tables{4};
%% 1.2
% Given investment, find the transition probability from omegatilde to
% omega with DAVh and omega with DAVl
N = size(Omega,1); % The number of states
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
% Transition probability from omegatilde to omega in the next period
P_dav = zeros(N,2); % First column represent DAVl, second represents DAVu
% probability that DAV becomes DAVl and DAVu respectively.
P_dav(:,2) = (DAVgrid_new*2+1-DAVl); 
P_dav(:,1) = 1-(2*DAVgrid_new+1-DAVl);


%% Find the transition probability from omegat to omegat
% the 5 dimensions represent: omegatilde, investment, DAV, violation,
% status
tran_ot2ot = zeros(45,20,3,2,2,...
    2,2,2,3);
for i=1:N % for each state in omega tilde
    % Investment probability at omegatilde i
    invp = Investprob(omega1(i),DAV(i)*2+1,status(i)+1 ...
        ,lag_inv(i)+1,vio(i)+1);
    for j = 1:2 % j=1,2 mean DAVl and DAVh respectively
        if j ==1
        % no investment
        tran_ot2ot(omega1(i),DAV(i)*2+1,status(i)+1 ...
        ,lag_inv(i)+1,vio(i)+1,1,j,:,:)=(1-invp).*P_dav(i,j).*...
            To2ot(omega1(i),DAVl(i),status(i)+1,1,lag_inv(i)+1,:,:);
        % investment
        tran_ot2ot(omega1(i),DAV(i)*2+1,status(i)+1 ...
        ,lag_inv(i)+1,vio(i)+1,2,j,:,:)=invp*P_dav(i,j)...
            *To2ot(omega1(i),DAVl(i),status(i)+1,2,lag_inv(i)+1,:,:);
        else
            % no investment
        tran_ot2ot(omega1(i),DAV(i)*2+1,status(i)+1 ...
        ,lag_inv(i)+1,vio(i)+1,1,j,:,:)=(1-invp)*P_dav(i,j)*...
            To2ot(omega1(i),DAVu(i),status(i)+1,1,lag_inv(i)+1,:,:);
        % investment
        tran_ot2ot(omega1(i),DAV(i)*2+1,status(i)+1 ...
        ,lag_inv(i)+1,vio(i)+1,2,j,:,:)=invp*P_dav(i,j)...
            *To2ot(omega1(i),DAVu(i),status(i)+1,2,lag_inv(i)+1,:,:);
        end
    end
end
end

