function [ dx, dy, ds, dtau, dkappa ] =  combinedir(blk, x, s, tau, kappa, ... 
   A, b, c, rp, rd, rg, mu, sigma, schur, Hx_inv, gx, Hx_inv_At)
% This function computes the Newton search direction in combined phase.
%   


eta = 1 - sigma;
k = length(blk);
n = sum(blk);
m = length(b);
rhsSchur = zeros(m,2);

Rd = eta*rd + s + sigma*mu*gx;

for i = 1:k
    ind_i = sum(blk(1:i-1))+1:sum(blk(1:i));
    rhsSchur = rhsSchur + Hx_inv_At{i}'*[Rd(ind_i) , c(ind_i)];
end


rhsSchur = rhsSchur + [eta*rp , b];
vq = schur\rhsSchur;
up = zeros(n,2);

clear ind_i
for i = 1:length(blk)
    ind_i = sum(blk(1:i-1))+1:sum(blk(1:i));
    up(ind_i,:) = Hx_inv_At{i}*vq - Hx_inv{i}*[Rd(ind_i) , c(ind_i)];
end

dtau = (c'*up(:,1) - b'*vq(:,1) + eta*rg - kappa + sigma*mu/tau)/...  % - dtau_p*dkappa_p/tau
         (-c'*up(:,2) + b'*vq(:,2) + kappa/tau);

dy = vq(:,1) + vq(:,2)*dtau;
dx = up(:,1) + up(:,2)*dtau;
ds = eta*rd - A'*dy + c*dtau;
dkappa = sigma*mu/tau - kappa - (kappa/tau)*dtau;  % - (kappa/tau)*dtau_c


end




