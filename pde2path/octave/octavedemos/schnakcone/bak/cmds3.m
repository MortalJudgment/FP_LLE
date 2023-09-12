close all; keep pphome; 
%% Schnakenberg cone, a=0.4 
dir='h3'; nx=40; dsmin=0.01; dsmax=0.1; ell=1.2; eps=0.05; 
lam=3.225; sig=0; D=60; a=0.4; s=0; 
p=[]; par=[lam sig D a s eps]; 
p=schnakcinit(p,nx,par,ell); p.np, p.pm.resfac=1e-4; p.sol.ds=-dsmin; p.sw.verb=2; 
p.nc.dsmax=dsmax; p.nc.dsmin=dsmin; p=setfn(p,dir); p.nc.mu2=0.01; p.nc.tol=1e-6; p0=p; 
%% dummy refine near 0 
p.nc.ngen=2; p.nc.maxt=20000; p.nref=300; p.fuha.e2rs=@e2rs_rad; p=oomeshada(p); p0=p; 
%%
p=p0; p.nc.neig=20; p.sw.jac=1; 
p=cont(p,10);
%% 1st branch 
p=swibra(dir,'bpt1','b1',0.01); p.nc.dsmax=0.25; p=cont(p,10); cplot(p); 
%% BP1, 
bp='bpt1'; aux=[]; aux.m=4; aux.isotol=0.25; 
aux.besw=0; %aux.ali=[1 3 5 7 9 11 13]; aux.besw=1; 
p0=qswibra(dir,bp,aux); p0.nc.dsmin=0.01; p0.sw.bifcheck=0; p0.nc.dsmax=0.05; 
p0.pm.resfac=1e-2; p0.sol.ds=0.05; 
%% centered, 
p=gentau(p0,[1],'c1a');  p.sol.ds=0.1; p=cont(p,3); cplot(p); pause; p=pmcont(p,15); cplot(p);
%% centered with stripe 
p=gentau(p0,[0 1],'c1b');  p.sol.ds=-0.1; p=cont(p,3); cplot(p); pause; p=pmcont(p,15); cplot(p);
%% somewhat stripe like 
p=gentau(p0,[0 0 1],'c1c');  p.sol.ds=-0.1; p=cont(p,3); cplot(p); pause; p=pmcont(p,15); cplot(p);
%% cntered again 
p=gentau(p0,[0 0 0 1],'c1d');  p.sol.ds=-0.1; p=cont(p,3); cplot(p); pause; p=pmcont(p,15); cplot(p);
%% plot BD 
fnr=3; figure(fnr); clf; pcmp=7; plotbra(dir,'pt10',fnr,pcmp,'cl','k','lsw',0);
plotbra('c1a','pt10',fnr,pcmp,'cl','r','lab',10); 
plotbra('c1b',fnr,pcmp,'cl','m','lab',10); 
plotbra('c1c',fnr,pcmp,'cl',p2pc('b1'),'lab',10); 
plotbra('c1d',fnr,pcmp,'cl',p2pc('b3'),'lab',10); 
ylabel('max u_1'); set(gca,'YTick',[3.2 3.6 4 4.4]);
%% solution plots 
p=loadp('c1a','pt10'); plotsol(p); cplot(p); pause 
p=loadp('c1b','pt10'); plotsol(p); cplot(p); pause 
p=loadp('c1c','pt10'); plotsol(p); cplot(p); pause 
p=loadp('c1d','pt10'); plotsol(p); cplot(p); 
