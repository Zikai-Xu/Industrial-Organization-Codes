function VarML = varml(x,ll,Omega,Omegat,Bellman,somega)
% Initialize the derivative of loglikelihood with respect to theta
N = size(somega,1);
dlogl = zeros(N,5);
I = eye(5)*1e-6;

for j = 1:5
    % plus 1e-6 to the jth element in theta
    newtheta = x + I(:,j); 
    [~,~,Investprob]=Bellmanfun(newtheta,Omega,Omegat,Bellman);
    [~,llnew] = LogLike(Investprob,somega); 
    dlogl(:,j) = (llnew-ll)./1e-6;% 7245 by 5 matrix

end
VarML = (dlogl.'*dlogl)^(-1);
end

