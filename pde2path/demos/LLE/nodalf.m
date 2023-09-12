function f=nodalf(p,u) % for LLE equation
u1 = u(1:p.np); u2 = u(p.np+1:2*p.np); 
par = u(p.nu+1:end); % f,zeta, d

f1 = -par(2)*u1 - u2 + (u1.^2 + u2.^2).*u1; 
f2 = u1 - par(2)*u2 + (u1.^2 + u2.^2).*u2 - par(1);
f = [f1; f2]; 