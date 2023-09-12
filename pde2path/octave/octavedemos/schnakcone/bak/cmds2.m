close all; keep pphome; 
%% Schnakenberg cone, a=0.5 
dir='h2'; nx=50; dsmin=0.01; dsmax=0.1; ell=1.2; eps=0.05; 
lam=3.225; sig=0; D=60; a=0.5; s=0; 
p=[]; par=[lam sig D a s eps]; 
p=schnakcinit(p,nx,par,ell); p.np, p.pm.resfac=1e-4; p.sol.ds=-dsmin; p.sw.verb=2; 
p.nc.dsmax=dsmax; p.nc.dsmin=dsmin; p=setfn(p,dir); p.nc.mu2=0.01; p.nc.tol=1e-6; p0=p; 
% dummy refine near 0 
p.nc.ngen=2; p.nc.maxt=20000; p.nref=300; p.fuha.e2rs=@e2rs_rad; p=oomeshada(p); p0=p; 
%%
p=p0; p.nc.neig=20; p.sw.jac=1; p=cont(p,10);
%% BP1, 
bp='bpt1'; aux=[]; aux.m=4; aux.isotol=0.25; aux.besw=0; 
p0=qswibra(dir,bp,aux); p0.nc.dsmin=0.01; p0.sw.bifcheck=0; p0.nc.dsmax=0.05; 
p0.pm.resfac=1e-2; p0.sol.ds=0.1; 
%% generating 3 branches via gentau 
p=gentau(p0,[1],'b1'); p=cont(p,5); cplot(p); p=pmcont(p,20); 
p=gentau(p0,[0 1],'b2'); p=cont(p,5); cplot(p); pause; p=pmcont(p,20); 
p=gentau(p0,[0 0 0 1],'b3'); p=cont(p,5); cplot(p); pause; p=pmcont(p,20); 
%% plot BD 
dir='h2'; 
fnr=3; figure(fnr); clf; pcmp=7; plotbra(dir,'pt20',fnr,pcmp,'cl','k','lsw',0);
plotbra('b1',fnr,pcmp,'cl','r','lab',10,'lp',15); 
plotbra('b2',fnr,pcmp,'cl','m','lab',10); 
plotbra('b3',fnr,pcmp,'cl',p2pc('b3'),'lab',10); 
ylabel('max u_1'); set(gca,'YTick',[3.2 3.6 4 4.4]);
axis([3 3.225 3 4.4]); 
%% solution plots 
p=loadp('b1','pt10'); plotsol(p); nola; cplot(p,1); pause 
p=loadp('b2','pt10'); plotsol(p); nola; cplot(p,1); pause 
p=loadp('b3','pt10'); plotsol(p); nola; cplot(p,1); 
%% 
p=swiparf('b1','pt10','b1a',4); p.usrlam=[0.25 0.5]; p.sol.ds=0.05; getaux(p)', 
p.nc.ntot=100; clf(2); p.nc.tol=1e-6; p.sw.foldcheck=0; 
p.nc.lammax=2; p.nc.lammin=0.25; p.nc.dsmax=0.1; p=pmcont(p,40); 
%% 
p=swiparf('b2','pt10','b2a',4); p.usrlam=[0.25 0.5]; p.sol.ds=0.01; getaux(p)', 
p.nc.ntot=100; clf(2); p.nc.tol=1e-8; p.sw.foldcheck=0; 
p.nc.lammax=2; p.nc.lammin=0.25; p.nc.dsmax=0.5; p=pmcont(p,40); 
%% 
p=swiparf('b3','pt10','b3a',4); p.usrlam=[0.25 0.5]; p.sol.ds=-0.01; getaux(p)', 
p.nc.ntot=100; clf(2); p.nc.tol=1e-8; p.sw.foldcheck=0; 
p.nc.lammax=2; p.nc.lammin=0.25; p.nc.dsmax=0.1; p=pmcont(p,40); 
%%
fnr=3; figure(fnr); clf; pcmp=0; 
plotbra('b1a','pt40',fnr,pcmp,'cl','r','lab',[0 30 40]);
plotbra('b2a','pt20',fnr,pcmp,'cl','m','labi',20); 
plotbra('b3a','pt40',fnr,pcmp,'cl','b','labi',40); 
ylabel('||u_1||_2'); 
%% solution plots 
dir='b1a'; iset=[0 30 40]; %dir='b2a'; iset=[0 20]; %dir='b3a'; iset=[0 40];
for i=1:length(iset)
    pt=['pt' mat2str(iset(i))]; 
    p=loadp(dir,pt); plotsol(p); cplot(p,2); pause 
end 
