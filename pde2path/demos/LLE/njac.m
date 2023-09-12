function [f1u,f1v,f2u,f2v]=njac(p,u) % Jacobian for LLE
% u1=u(1:p.np); u2=u(p.np+1:2*p.np); 
u1=u(1:p.nu/2); u2=u(p.nu/2+1:p.nu);
par = u(p.nu+1:end); % f,zeta,d

f1u = -par(2) + 3*u1.^2 + u2.^2;
f1v = -1 + 2*u1.*u2;
f2u = 1 + 2*u1.*u2;
f2v = -par(2) + u1.^2 + 3*u2.^2;