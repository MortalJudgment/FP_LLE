function r = sG(p,u) % LLE
f = nodalf(p,u); par = u(p.nu+1:end); % f,zeta,d
K = kron([[par(3),0];[0,par(3)]],p.mat.K); 
r = K*u(1:p.nu) - p.mat.M*f; 