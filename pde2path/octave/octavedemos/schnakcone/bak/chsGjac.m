function Gu=sGjac(p,u)  % AC Jacobian 
par=u(p.nu+1:end); eps=par(2); a=par(3);  s=par(5); 
u=u(1:p.nu); n=p.np; % split u into parameters and PDE variables 
fu=1-3*u.^2; Fu=spdiags(fu,0,n,n);  % convert derivatives into (sparse) matrix 
if p.nc.ilam(1)==3; [p.mat.K,p.mat.M]=LBc(p,a); end % update LBO
Gu=eps^2*p.mat.K-p.mat.M*Fu+s*p.mat.Krot;        % the Jacobian matrix 