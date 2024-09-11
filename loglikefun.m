function loglike = loglikefun(Theta,Omega,Omegat,Bellman,somega)
%LOGLIKEFUN Summary of this function goes here
%   Detailed explanation goes here
[~,~,Investprob]=Bellmanfun(Theta,Omega,Omegat,Bellman);
[loglike,~] = LogLike(Investprob,somega);
end

