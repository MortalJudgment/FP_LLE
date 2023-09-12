po=getpte(p); x=po(1,:)'; y=po(2,:)'; n=p.nu; M=p.mat.M; 
p.u(1:n)=x.^2.*y; plotsol(p,1,1,2); 
p.u(1:n)=1\(p.mat.Kxy1*p.u(1:n)); plotsol(p,10,1,2); axis tight; 