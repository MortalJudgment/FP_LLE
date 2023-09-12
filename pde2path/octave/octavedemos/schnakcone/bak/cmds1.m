close all; keep pphome; 
%% Schnakenberg cone, a=0.75 
dir='h1'; nx=40; dsmin=0.01; dsmax=0.1; ell=1.2; eps=0.05; 
lam=3.225; sig=0; D=60; a=0.75; s=0; 
p=[]; par=[lam sig D a s eps]; 
p=schnakcinit(p,nx,par,ell); p.np, p.pm.resfac=1e-4; p.sol.ds=-dsmin; p.sw.verb=2; 
p.nc.dsmax=dsmax; p.nc.dsmin=dsmin; p=setfn(p,dir); p.nc.mu2=0.01; p.nc.tol=1e-6; p0=p; 
% dummy refine near 0 
p.nc.ngen=2; p.nc.maxt=20000; p.nref=300; p.fuha.e2rs=@e2rs_rad; p=oomeshada(p); p0=p; 
%%
p=p0; p.nc.neig=20; p.sw.jac=1; p=cont(p,10);
%% BP1, 
bp='bpt1'; aux=[]; aux.m=6; aux.besw=0;
p0=qswibra(dir,bp,aux); p0.nc.dsmin=0.01; p0.sw.bifcheck=0; p0.nc.dsmax=0.05; 
p0.pm.resfac=1e-2; p0.sol.ds=0.05; 
%% centered 
p=gentau(p0,[1],'a1a');  p.sol.ds=0.1; p=cont(p,20); cplot(p); p=pmcont(p,20); cplot(p);
%% centered 
p=gentau(p0,[0 1],'a1b');  p.sol.ds=0.1; p=cont(p,5); cplot(p);  p=pmcont(p,20); cplot(p);
%% somewhat stripe like 
p=gentau(p0,[0 0 1],'a1c');  p.sol.ds=0.1; p=cont(p,5); cplot(p);p=pmcont(p,20); cplot(p);
%% cntered again 
p=gentau(p0,[0 0 0 1],'a1d');  p.sol.ds=0.1; p=cont(p,20); cplot(p); p=pmcont(p,20); cplot(p);
%% plot BD 
dir='h1'; 
fnr=3; figure(fnr); clf; pcmp=7; plotbra(dir,'pt10',fnr,pcmp,'cl','k','lsw',0);
plotbra('a1a',fnr,pcmp,'cl','r','lab',10); 
plotbra('a1b','pt20',fnr,pcmp,'cl','m','lab',10); 
plotbra('a1c',fnr,pcmp,'cl',p2pc('b1'),'lab',20); 
plotbra('a1d',fnr,pcmp,'cl',p2pc('b3'),'lab',20); 
ylabel('max u_1'); set(gca,'YTick',[3.2 3.6 4 4.4]);
%% solution plots 
p=loadp('a1a','pt20'); plotsol(p); cplot(p); pause 
p=loadp('a1b','pt20'); plotsol(p); cplot(p); pause 
p=loadp('a1c','pt20'); plotsol(p); cplot(p); pause 
p=loadp('a1d','pt20'); plotsol(p); cplot(p); 
%% cont in a=base-radius, smaller a 
p=swiparf('a1a','pt10','a1a-10a',4); clf(2); p.usrlam=[0.25 0.5]; p.sol.ds=-0.01; getaux(p)', 
p.nc.lammax=2; p.nc.lammin=0.25; p.nc.dsmax=0.05; p=cont(p,20); 
%% cont in a=base-radius, smaller a 
p=swiparf('a1b','pt10','a1b-10a',4); p.usrlam=[0.25 0.5]; p.sol.ds=-0.01; getaux(p)', 
p.nc.lammax=2; p.nc.lammin=0.25; p.nc.dsmax=0.05; p=cont(p,20); 
%% 
figure(3); clf; plotbra('a1a-10a','pt20'); 
%%
p=loadp('a1a-10a','pt10'); plotsol(p); cplot(p); pause 
p=loadp('a1a-10a','pt20'); plotsol(p); cplot(p); 