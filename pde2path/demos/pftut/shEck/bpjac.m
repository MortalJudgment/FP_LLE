function duGuph=bpjac(p,u) % for SH, gen.form, but only top left block nonzero 
n=p.np; u1=u(1:p.np); par=u(2*p.nu+1:end); nup=par(2); ov=ones(n,1); 
f1uu=2*nup*ov-6*u1; % only nonzero entry 
f1uv=0*ov; f1vv=f1uv; f2uu=f1uv; f2uv=f1uv; f2vv=f1uv;
ph1=u(p.nu+1:p.nu+p.np); ph2=u(p.nu+p.np+1:2*p.nu);
M1=spdiags(f1uu.*ph1+f2uu.*ph2,0,n,n); M2=spdiags(f1uv.*ph1+f2uv.*ph2,0,n,n); 
M3=spdiags(f1uv.*ph1+f2uv.*ph2,0,n,n); M4=spdiags(f1vv.*ph1+f2vv.*ph2,0,n,n);
duGuph=-[[M1 M2]; [M3 M4]]*p.mat.M; 