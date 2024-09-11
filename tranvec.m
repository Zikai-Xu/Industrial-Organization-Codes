function P = tranvec(statefrom,tran_ot2ot,Tables)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Read tables
transitionprob = Tables{1};
% Given the initial state, find the distribution over state-to
DAVgridto = transitionprob.DAVgridto;
ordervioto = transitionprob.orderedvio_to;
lag_invto = transitionprob.laginv_to;
violationto = transitionprob.violationto;
M = size(transitionprob,1); %161
P = zeros(M,1);
% read the initial state
omega1 = statefrom(1);
davgrid = statefrom(2);
status = statefrom(3);
laginv = statefrom(4);
vio = statefrom(5);
% dav interpolate
davu = min(ceil((davgrid/2*0.9 + vio)*2),20);
davl = floor((davgrid/2*0.9 + vio)*2);
for i = 1:M
    if ordervioto(i)==0 % if compliance, then vioto = 0, DAVto = 0
        % lag_invto = 1
        P(i) = sum(tran_ot2ot(omega1,davgrid+1,status+1,laginv+1,vio+1,:,:,:,1),"all");
    elseif DAVgridto(i) == davl
        P(i) = tran_ot2ot(omega1,...
        davgrid+1,status+1,laginv+1,vio+1,lag_invto(i)+1,1,violationto(i)+1,ordervioto(i)+1);
    elseif DAVgridto(i) == davu
        P(i) = tran_ot2ot(omega1,...
        davgrid+1,status+1,laginv+1,vio+1,lag_invto(i)+1,2,violationto(i)+1,ordervioto(i)+1);
    end
end
end

