function dia=rocdiagn(p,wnr,varargin) % standard diagnostics, plot u_1, lam_1
% use as template to plot adapted diagnostics 
f1=1; f2=1; try; f1=varargin{1}; f2=varargin{2}; catch; end; 
s0=p.oc.s0; s1=p.oc.s1; u0=p.oc.u0; u1=p.oc.u1;
sol=p.cp; T=sol.par(1); 
t=T*sol.t; y=sol.u; sl=length(t); n=s0.nu/4; dia=zeros(8,sl); 
for i=1:sl
    dia(1,i)=norm(y(1:n,i)-u0(1:n),'inf'); % diff R from init 
    dia(2,i)=norm(y(1:n,i)-u1(1:n),'inf'); % diff R from target
    dia(3,i)=norm(y(n+1:2*n,i)-u0(n+1:2*n),'inf'); % diff B from IC
    dia(4,i)=norm(y(n+1:2*n,i)-u1(n+1:2*n),'inf'); % diff B from target  
    dia(5,i)=norm(y(2*n+1:3*n,i)-u1(2*n+1:3*n),'inf'); % diff x from target  
    dia(6,i)=norm(y(3*n+1:4*n,i)-u1(3*n+1:4*n),'inf'); % diff y from target   
end
figure(wnr); set(gca,'FontSize',s0.plot.fs); 
if 0
plot(t,dia(1,:),'k',t,dia(2,:),'b',t,dia(3,:),'m',t,dia(4,:),'g','LineWidth',2); 
legend('||R-R_0||_\infty','||R-R_1||_\infty','||B-B_0||_\infty','||B-B_1||_\infty'); 
else
 plot(t,dia(2,:),'k',t,dia(4,:),'b',t,f1*dia(5,:),'m',t,f2*dia(6,:),'g','LineWidth',2); 
legend('||R-R_1||_\infty','||B-B_1||_\infty',...
    [mat2str(f1) '*||x-x_1||_\infty]'],[mat2str(f2) '*||y-y_1||_\infty]']);    
end
axis tight; xlabel('t'); set(gca,'fontsize',14); 
rho=s1.u(s1.nu+p.oc.rhoi); 
[Jtran,jcv,jcvd]=jcaiT(s1,sol,rho); Jtran,rho
Jh=s1.branch(23,end); 
jp=Jtran+exp(-rho*T)*Jh; 
figure(12); clf; 
 plot(t,jcv,'k',t,jcvd,'b'); axis tight; legend('J_{ca}','e^{-\rho t}J_{ca}'); 
 title(['J=' mat2str(jp,4)]); 
axis tight; xlabel('t'); set(gca,'fontsize',14); 
end 