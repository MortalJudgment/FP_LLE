function psol3Dm_NT(p,sol0,sol,wnr,cmp,tit) % plot sol0 and sol into one fig, e.g., for Skiba
sl=length(sol.t); figure(wnr); clf; hold on;  [po,~,~]=getpte(p); 
nx=p.np; x=po(1,1:nx)'; v1=ones(nx,1); t=v1*sol.par(1); 
for i=1:sl
  if cmp==0;  par=p.u(p.nu+1:end); u=[sol.u(:,i);par]; z=p.fuha.con(p,u); 
  else cs=(cmp-1)*p.np; z=sol.u(cs+1:cs+nx,i); end 
  plot3(x,sol.t(i)*t,z(1:nx),'r'); 
end
sl0=length(sol0.t); 
nx=p.np; x=po(1,1:nx)'; v1=ones(nx,1); t=v1*sol.par(1); 
for i=1:sl0
  if cmp==0; par=p.u(p.nu+1:end); u=[sol0.u(:,i);par]; z=p.fuha.con(p,u); 
  else cs=(cmp-1)*p.np; z=sol0.u(cs+1:cs+nx,i); end 
  plot3(x,sol0.t(i)*t,z(1:nx)); 
end
view(15,40); axis tight; grid on; 
set(gca,'FontSize',p.plot.fs); title(tit);