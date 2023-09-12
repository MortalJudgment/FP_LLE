close all; keep pphome; 
%% CH on cone 
dir='h'; nx=40; %nx=20;  
dsmin=0.005; dsmax=0.05; eps=0.1; m=-0.7; a=0.75; lam=0; s=0; ell=1.05; 
p=[]; par=[m eps a lam s]; p=chcinit(p,1,nx,par,ell); huclean(p); 
p.sw.verb=2; p.nc.dsmax=dsmax; p.nc.dsmin=dsmin; p=setfn(p,dir); 
p.file.smod=5; p.nc.ilam=[1 4]; p.sw.bifcheck=2;  p.tau=0*p.u;plotsol(p,1,1,1);
p.sw.foldcheck=0;  p0=p; 
%% dummy refine near 0 
p.nc.ngen=2; p.nc.maxt=20000; p.nref=300; p.fuha.e2rs=@e2rs_rad; p=oomeshada(p); p0=p; 
%% cont hom branch, with mass constraint 
p=p0; p.nc.mu1=1; p.nc.nq=1; p=findbif(p,4); %p=cont(p,8); %p=cont(p,20); 
%% 1st branches 
p=swibra(dir,'bpt1','b1',-0.01); p.nc.dsmax=0.25; p=cont(p,40); cplot(p); 
%%
p=swibra(dir,'bpt3','b3',0.01); p.nc.dsmax=0.25; p=cont(p,40); cplot(p); 
%%
figure(3); c=6; clf; plotbra('h','pt15',3,c,'cl','k','lsw',0);
plotbra('b1','pt40',3,c,'cl','b');  
%plotbra('b2','pt25',3,c,'cl','r','lab',[16,20]); 
plotbra('b3','pt25',3,c,'cl','m','lab',[15,21]); 
ylabel('E_\epsilon'); %axis([-0.9 0.1 0.8 6]); 
%%
p=loadp('b1','pt14'); plotsol(p); nolti; cplot(p); pause 
p=loadp('b1','pt21'); plotsol(p); nolti; cplot(p); pause 
p=loadp('b1','pt29'); plotsol(p); nolti; cplot(p); pause 
p=loadp('b1','pt37'); plotsol(p); nolti; cplot(p);
%%
p=loadp('b3','pt18'); plotsol(p); nolti; cplot(p); 
%% meshada and cont in eps
v=[30 60]; p=swiparf('b1','pt21','b1e1',[2 4]); p.sol.ds=-0.01; plotsol(p,1,1,1); view(v); 
p.fuha.e2rs=@e2rs; E1=chE(p,p.u); p.nc.ngen=2; p.nc.sig=0.5; p=oomeshada(p); view(v);  E2=chE(p,p.u); E1,E2
clf(2); p.nc.lammax=1; p.nc.lammin=0.05; p.usrlam=0.1; p=cont(p,10); 
%% meshada and cont in eps
v=[30 60]; p=swiparf('b1','pt29','b1e2',[2 4]); p.sol.ds=-0.01; plotsol(p,1,1,1); view(v); 
p.fuha.e2rs=@e2rs; E1=chE(p,p.u); p.nc.ngen=2; p.nc.sig=0.5; p=oomeshada(p); view(v);  E2=chE(p,p.u); E1,E2
clf(2); p.nc.lammax=1; p.nc.lammin=0.05; p.usrlam=0.1; p=cont(p,10); 
%% meshada and cont in eps
p=swiparf('b3','pt18','b3e',[2 4]); p.sol.ds=-0.01; getaux(p)',pause 
%p.fuha.e2rs=@e2rs; p.nc.ngen=2; p.nc.sig=0.5; p=oomeshada(p); view(v); p.nc.tol=1e-6; 
p.nc.lammax=1; p.nc.lammin=0.05; p.usrlam=0.1; p=cont(p,10); 
%% bd and soln plots in eps 
figure(3); clf; c=6; plotbra('b1e1','pt8',3,c,'cl','b','lab',7); 
plotbra('b1e2','pt7',3,c,'cl','r','lab',7);  
plotbra('b3e','pt7',3,c,'cl','m','lab',7);  
xlabel('\epsilon'); ylabel('E_\epsilon'); 
%% soln plots, eps-cont 
p=loadp('b1e1','pt7'); plotsol(p); cplot(p); pause 
p=loadp('b1e','pt7'); plotsol(p); cplot(p);  pause 
p=loadp('b3e','pt7'); plotsol(p); cplot(p); 
%% analytical expr. for interface lengths (straight lines, circles) 
a=0.75; 2*sqrt(1+a^2), 2*pi/sqrt(2)*a;  
a=0.5:0.01:1; figure(10); clf; plot(a,2*sqrt(1+a.^2),'r'); 
%% cont in a=base-radius, larger a, somewhat boring (straight lines win) 
p=swiparf('b1','pt21','b1a1',[3 4]); p.usrlam=[0.3 0.5 1 2]; clf(2); p.sol.ds=0.1; getaux(p)', 
p.nc.lammax=2; p.nc.lammin=0.25;  p=cont(p,11); 
p=swiparf('b1','pt27','b1a2',[3 4]); p.usrlam=[0.3 0.5 1 2]; clf(2); p.sol.ds=0.1; getaux(p)', 
p.nc.lammax=2; p.nc.lammin=0.25;  p=cont(p,11); 
%% cont in a=base-radius, smaller a 
p=swiparf('b1','pt21','b1a1b',[3 4]); clf(2); p.usrlam=[0.3 0.5 1 2]; p.sol.ds=-0.01; getaux(p)', 
p.nc.lammax=2; p.nc.lammin=0.25; p.nc.dsmax=0.05; p=cont(p,20); 
p.nc.dsmax=0.1; p.nc.tol=1e-6; p=cont(p,10); 
p=swiparf('b1','pt27','b1a2b',[3 4]); p.usrlam=[0.3 0.5 1 2]; p.sol.ds=-0.01; getaux(p)', 
p.nc.lammax=2; p.nc.lammin=0.25; p.nc.dsmax=0.05; p=cont(p,11); 
%% BD plot, a-cont 
figure(3); clf; c=6; plotbra('b1a1b','pt34',3,c,'cl',p2pc('b2')); 
plotbra('b1a2b',3,c,'cl',p2pc('r2')); ylabel('E_\epsilon'); 
%% soln plots 
p=loadp('b1a1b','pt11'); plotsol(p); nolti; cplot(p); pause 
p=loadp('b1a1b','pt22'); plotsol(p); nolti; cplot(p); pause 
p=loadp('b1a1b','pt27'); plotsol(p); nolti; cplot(p); pause 
p=loadp('b1a1b','pt32'); plotsol(p); nolti; cplot(p); pause 
p=loadp('b1a2b','pt9'); plotsol(p); nolti; cplot(p); pause 
p=loadp('b1a2b','pt17'); plotsol(p); nolti; cplot(p);
