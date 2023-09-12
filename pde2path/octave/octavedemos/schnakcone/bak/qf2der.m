function qu=qf2der(p,u)  
a=u(p.nu+3); vol=pi; %vol=a*sqrt(1+a^2)*pi; 
q1u=p.mat.vM/vol; 
if p.nc.nq==1; qu=q1u; 
else uold=p.u(1:p.nu); uox=p.mat.Krot*uold; q2u=uox'; qu=[q1u;q2u]; 
end
