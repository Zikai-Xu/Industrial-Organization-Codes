% Parameters
beta = 0.9; % discount factor
A = [0.8,1.2];

% Steady State
% Euler equation gives 1 = beta(0.3Kss^{-0.7}+0.3)
Kss = ((1/beta/0.3)-1)^(-1/0.7);

% 500 grids around steady state
Kgrid = linspace(0,2*Kss,500);
N = 500*2; % number of states
% Initializing value function
V = zeros(500,2);
% First input is K and second input is z
% Policy function
Copt = zeros(500,2);
%% Iteration
iter = 0;
norm = 1;
oldV = V;
newV = V; % They are used to compute norm
while norm>1e-6
    iter = iter +1;
for i=1:N
    % Find the current state 
if i <=500
    kindex = i;
    z = 1;
    kt = Kgrid(kindex);
    At = A(z);
else
    kindex = i-500;
    z = 2;
    kt = Kgrid(kindex);
    At = A(z);
end

% Find the consumption given kt,zt and each kt+1
C = (At*kt^(0.3)+0.3*kt)*ones(500,1) - Kgrid.';
C(C==0) = 1e-6;
C(C<0) = NaN;
u = log(C); % flow utility at current state


% Find the new value function

aux = u+beta*(0.5*oldV(:,1)+0.5*oldV(:,2));
[newV(kindex,z),I] = max(aux);
% Optimal consumption rule
Copt(kindex,z) = At*kt^(0.3)+0.3*kt - Kgrid(I);
end
norm = max(abs(newV - oldV),[],'all');
oldV = newV; % Update value function
disp([norm,iter]);
end


