function pi = steadypi(Investprob,Tables)
%STEADYPI Summary of this function goes here
%   Detailed explanation goes here
% Read tables
transitionprob = Tables{1};
N1= 45; % size of fixed states
M = size(transitionprob,1); % 161 number of states in each fixed state
% Find P_ij given each fixed state
P = zeros(N1,M,M);
pi = zeros(N1,M);
b = zeros(M+1,1); b(end)=1;
Davgrid = transitionprob.DAVgridto;
Status = transitionprob.orderedvio_to;
Lag_inv = transitionprob.laginv_to;
Vio = transitionprob.violationto;
tran_ot2ot = tranmat(Investprob,Tables);
for k=1:N1
    % Omega1 = k
    for j=1:M
        statefrom = [k,Davgrid(j),Status(j),Lag_inv(j),Vio(j)];        
        P(k,j,:) = tranvec(statefrom,tran_ot2ot,Tables);
    end
    % Each P(i,:,:) is a 161 by 161 matrix
    % need to find pi(i,:)
    PP = reshape(P(k,:,:),M,M);
    A = [PP.'- eye(M);ones(1,M)]; % A is of size M+1 by M
    pi(k,:) = (A.'*A)\A.'*b;
end
end

