%% Demo for imperfect pitchforks 
close all; keep pphome; 
%% c2: init (generic), then specific settings (could also be set in init)
p=[]; par=[1 -0.2 1 0.1 0]; p=acinit(p,5,40,par); p.nc.dsmax=0.05; 
ulam=[0 0.5]; p.usrlam=ulam; p=setfn(p,'i1b'); cont(p,10);
%% use deflation
p=loadp('sq1st','pt10'); x=getpte(p); x=x'; x=[x; x]; x=p.mat.drop*x; 
x2=max(x(:,1)); y2=max(x(:,2)); 
p.osol=p.u; u0=p.u; nu=p.nu; p.nc.nq=0; 
u0(1:p.nu)=p.u(1:p.nu)-0.5*cos(pi/(x2)*x(:,1)); % set initial guess for Newton 
u0(1:p.nu)=0*p.u(1:p.nu)+1; 
%p.u=u0; plotsol(p); pause 
p.sw.jac=2; p.deflq=2; p.deflal1=0e-2; p.deflal2=1; %p.nc.del=1e-6; 
p.nd=1; [u, p]=deflsol(p,u0); p.u=u; plotsol(p);
%% run deflated Newton, use result for next deflation 
wi=1:2:3; % long
wi=4:6; % shorter
wi=9:15; 
fc=rand(length(ww),1)-0.5; ug=zeros(nu,1); 
ww=[[0.25;0] [0; 0.25] [0.25; 0.25] ...
    [0.5;0] [0; 0.5] [0.5; 0.5] ... % 4-6
    [0.5;0.25] [0.25; 0.5] [1;0] [0;1] [1;1] ... % 9-11
    [1; 0.25] [0.25;1] [1; 0.5] [0.5;1] ... % 12-15
    [2;0] [0;2] [2;1] [1;2] [2;2]]; 
wnr=sum(ww.^2,1); 
%ww=[[1;0] [0;1] [1;1] [2;0] [0;2] [2;1] [1;2] [2;2]]; 
p.nc.imax=20; 
for i=wi;  % get ICs from random Fourier
     sx=x2*(rand(1)-0.5); sy=y2*(rand(1)-0.5); 
    for j=1:nu/2       
        ug(j)=ug(j)+fc(i)*cos((x(j,:)-[sx, sy])*ww(:,i));
        ug(j+nu/2)=ug(j+nu/2)-wnr(i)*fc(i)*cos((x(j,:)-[sx, sy])*ww(:,i));
    end
end 
p.nc.tol=1e-10; p.deflq=2; 
u0(1:p.nu)=0*p.u(1:p.nu)+ug; p.u=u0; plotsol(p); pause
[u, p]=deflsol(p,u0); p.u=u; plotsol(p); 
p.nc.dsmax=0.1; p1=p; % save p in p1 to start cont. runs
%%
p.nd, for i=1:p.nd; p.u=p.osol(:,i); plotsol(p); title(mat2str(i)); pause; end
%% run cont on obtained solutions
p=postdefl(p1,8,'j2a',0.1); p=cont(p,20); 
%%
p=postdefl(p1,2,'j2b',-0.1); p=cont(p,20);
p=postdefl(p1,4,'j3a',0.1); p.sw.bifcheck=0; p=cont(p,20);
p=postdefl(p1,4,'j3b',-0.1); p.sw.bifcheck=0; p=cont(p,20);
%% 
f=3; c=0; figure(f); clf;
plotbra('j2a',f,c,'cl','r'); 
plotbra('j2b',f,c,'cl','r');
plotbra('j3a',f,c,'cl','b'); 
plotbra('j3b',f,c,'cl','b');
axis([-0.2 0.51 0.1 3.5]); 
%% jaccheck 
p=loadp('i1','pt6'); x=getpte(p); x=x'; x2=max(x); % plotsol(p); pause;
p.osol=p.u; u0(1:p.nu)=p.u(1:p.nu)+0.2*cos(pi/(2*x2)*x)+0.05; p.u=u0; 
p.nc.del=1e-6; p.deflq=2; p.deflal1=1e-2; p.deflal2=1; p.nd=1; 
p.fuha.sGb=p.fuha.sG; p.fuha.sGjacb=p.fuha.sGjac;  % mod rhs
p.fuha.sG=@deflsG; p.fuha.sGjac=@deflsGjac; 
[Gu, Gun]=jaccheck(p); Gud=abs(Gu-Gun); e1=max(max(Gud)); 
spy(Gud>e1/1000); 
p.fuha.sG=p.fuha.sGb; p.fuha.sGjac=p.fuha.sGjacb; 


