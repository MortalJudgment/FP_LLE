function E=chE(p,u)  % energy for CH  on cone 
n=p.nu; par=u(n+1:end); eps=par(2); a=par(3); ga=1; dS=surfelem(p,u); % dS=1; 
u=u(1:n); W=ga*0.25*(u.^2-1).^2; sig=sqrt(2*ga)/3; 
[K,~]=LBc(p,a); s2=(K*u).*u.*dS; s2=sum(s2,1); 
E1=p.mat.vM*(W.*dS)/eps+0.5*eps*s2; E=E1/(2*sig); 