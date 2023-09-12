function Guuph=spjac(p,u)
u1 = u(1:p.np); u2 = u(p.np+1:2*p.np); par = u(2*p.nu+1:end);
n = p.np; 
f1uu = 6*u1; f1uv = 2*u2; f1vv = 2*u2;
f2uu = 2*u1; f2uv = 2*u1; f2vv = 6*u2;
ph1 = u(p.nu+1:p.nu+p.np); ph2 = u(p.nu+p.np+1:2*p.nu);
M1 = spdiags(f1uu.*ph1+f1uv.*ph2,0,n,n); 
M2 = spdiags(f1uv.*ph1+f1vv.*ph2,0,n,n); 
M3 = spdiags(f2uu.*ph1+f2uv.*ph2,0,n,n); 
M4 = spdiags(f2uv.*ph1+f2vv.*ph2,0,n,n);
Guuph = -p.mat.M*[[M1 M2];[M3 M4]]; 