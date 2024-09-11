function [loglike,ll] = LogLike(Investprob,somega)
% Number of observations
N = size(somega,1);
loglike = 0;
ll = zeros(N,1);
for i = 1:N
    omega1 = somega.omega1(i);
    DAV = somega.DAV(i);
    DAVh = min(ceil(DAV*2)+1,20);
    DAVl = floor(DAV*2)+1;
    status = somega.ordered_violator(i)+1;
    laginv = somega.lag_investment(i)+1;
    vio = somega.violation(i)+1;
    X = somega.investment(i);
    Likelihoodh = X*Investprob(omega1,DAVh,status,laginv,vio) +...
        (1-X)*(1-Investprob(omega1,DAVh,status,laginv,vio));
    Likelihoodl = X*Investprob(omega1,DAVl,status,laginv,vio) +...
        (1-X)*(1-Investprob(omega1,DAVl,status,laginv,vio));
    Likelihood = Likelihoodl + (DAV*2+1-DAVl)*(Likelihoodh-Likelihoodl);
    ll(i) = log(Likelihood);
    loglike = loglike + log(Likelihood);
end
end

