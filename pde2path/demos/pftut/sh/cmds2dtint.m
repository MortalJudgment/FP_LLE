%% patterns from initial guesses; long 2D hex-domain; init and zero-branch 
p=[]; lx=4*pi;nx=round(3*lx);ly=2*pi/sqrt(3);lam=0.2;nu=1.3;par=[lam; nu]; sw.sym=2;
sw.ref=1; ndim=2; p=shinit(p,nx,lx,ly,ndim,par,sw); dir='tint'; p=setfn(p,dir); 
po=getpte(p); x=po(1,:)'; y=po(2,:)'; % extract coord from p 
%% a 'good' initial guess for hex2str front, and hex2str front from Newton loop 
p.u(1:p.np)=(cos(x)+cos(x/2).*cos(sqrt(3)*y/2)).*(x<=0)+2*cos(x).*(x>0); 
spl(p,''); title('initial guess 1'); r=norm(resi(p,p.u),'inf'); 
[u,res,iter]=nloop(p,p.u); p.u(1:p.nu)=u(1:p.nu); 
fprintf('initial res=%g, res=%g after %i iteration\n',r,res,iter); 
spl(p,''); title('solution 1'); pause, clf(2); p=pmcont(p,20); % plot, then cont 
%% a 'bad' initial guess for hex2str front, Newton loop goes to stripes 
u0=(cos(x)+cos(x/2).*cos(sqrt(3)*y/2)).*(x<=0)+4*cos(x).*(x>0); p.u(1:p.np)=u0; 
spl(p,''); title('initial guess 2'); r=norm(resi(p,p.u),'inf'); 
[u,res,iter,Gu,Glam,p]=nloop(p,p.u); p.u(1:p.nu)=u(1:p.nu); 
fprintf('initial res=%g, res=%g after %i iteration\n',r,res,iter); 
spl(p,''); title('solution 2'); 
%% to obtain hex2str front, do a few steps with tintxs: preparations 
p.u(1:p.np)=u0; t1=0;ts=[];nc=0;dt=0.01;nt=500;pmod=20;smod=100; p.mat.Kadv=0; 
%% the tint loop, repeat this cell until residual is small (here just once) 
[p,t1,ts,nc]=tintxs(p,t1,ts,dt,nt,nc,pmod,smod,@nodalf);  
%% plot time series of res 
tss=5; plot(ts(1,tss:end),ts(2,tss:end));axis tight;legend('res');xlabel('t'); 
%% Newton loop after tint, then cont 
[u,res,iter]=nloop(p,p.u); p.u(1:p.nu)=u(1:p.nu);
fprintf('initial res=%g, res=%g after %i iteration\n',r,res,iter); 
plotsol(p,1,1,2); pause; clf(2); p=pmcont(p,20);