function r=sG(p,u)  % PDE rhs for CH 
par=u(p.nu+1:end); eps=par(2); a=par(3); lam=par(4); s=par(5); u=u(1:p.nu); 
if p.nc.ilam(1)==3;[p.mat.K,p.mat.M]=LBc(p,a); end % update LBO
f=u+lam-u.^3; % ga, pause 
r=eps^2*p.mat.K*u-p.mat.M*f+s*p.mat.Krot*u;  
