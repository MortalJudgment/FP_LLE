function q=qf2(p,u) 
m=u(p.nu+1); a=u(p.nu+3); vol=pi; %a*sqrt(1+a^2); vol=a^2*(1+a^2)*pi; 
u=u(1:p.nu); q1=p.mat.vM*u/vol-m;
if p.nc.nq==1; q=q1; 
else uold=p.u(1:p.nu); uox=p.mat.Krot*uold;
    q2=uox'*(u-uold); q=[q1; q2]; 
end  