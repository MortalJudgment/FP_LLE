close all; keep pphome; 
%% Schnakenberg cone, a=0.4 
dir='h4'; nx=40; dsmin=0.01; dsmax=0.1; ell=1.2; eps=0.05; 
lam=3.225; sig=0; D=60; a=0.4; s=0; 
p=[]; par=[lam sig D a s eps]; 
p=schnakcinit(p,nx,par,ell); p.np, p.pm.resfac=1e-4; p.sol.ds=-dsmin; p.sw.verb=2; 
p.nc.dsmax=dsmax; p.nc.dsmin=dsmin; p=setfn(p,dir); p.nc.mu2=0.01; p.nc.tol=1e-6; p0=p; 
%% dummy refine near 0 
p.nc.ngen=2; p.nc.maxt=20000; p.nref=300; p.fuha.e2rs=@e2rs_rad; p=oomeshada(p); p0=p; 
%%
p=p0; p.nc.neig=20; p.sw.jac=1; p=cont(p,10);
%% BP1, 
bp='bpt1'; aux=[]; aux.m=6; aux.isotol=0.25; 
aux.besw=0; %aux.ali=[1 3 5 7 9 11 13]; aux.besw=1; 
p0=qswibra(dir,bp,aux); p0.nc.dsmin=0.01; p0.sw.bifcheck=0; p0.nc.dsmax=0.05; 
p0.pm.resfac=1e-2; p0.sol.ds=0.05; 
%% centered, 
p=gentau(p0,[1],'d1a');  p.sol.ds=0.1; p=cont(p,5); cplot(p); p=pmcont(p,20); cplot(p);
%% centered with stripe 
p=gentau(p0,[0 1],'d1b'); p.nc.intol=-1e-5;  p.sol.ds=0.1; p=cont(p,5); cplot(p); p=pmcont(p,20); cplot(p);
%% somewhat stripe like 
p=gentau(p0,[0 0 1],'d1c');  p.sol.ds=0.1; p=cont(p,5); cplot(p); p=pmcont(p,20); cplot(p);
%% cntered again 
p=gentau(p0,[0 0 0 1],'d1d');  p.sol.ds=0.1; p=cont(p,5); cplot(p); p=pmcont(p,20); cplot(p);
%% plot BD 
fnr=3; figure(fnr); clf; pcmp=0; plotbra(dir,'pt10',fnr,pcmp,'cl','k','lsw',0);
plotbra('d1a',fnr,pcmp,'cl','r','lab',20); 
plotbra('d1b',fnr,pcmp,'cl','m','lab',20); 
plotbra('d1c',fnr,pcmp,'cl',p2pc('b1'),'lab',10); 
plotbra('d1d',fnr,pcmp,'cl',p2pc('b3'),'lab',20); 
ylabel('max u_1'); %set(gca,'YTick',[3.2 3.6 4 4.4]);
%% solution plots 
p=loadp('d1a','pt20'); plotsol(p); cplot(p); pause 
p=loadp('d1b','pt20'); plotsol(p); cplot(p); pause 
p=loadp('d1c','pt20'); plotsol(p); cplot(p); pause 
p=loadp('d1d','pt20'); plotsol(p); cplot(p); 
%% later BP, 
bp='bpt3'; aux=[]; aux.m=6; aux.isotol=0.25; 
aux.besw=0; %aux.ali=[1 3 5 7 9 11 13]; aux.besw=1; 
p0=qswibra(dir,bp,aux); p0.nc.dsmin=0.01; p0.sw.bifcheck=0; p0.nc.dsmax=0.05; 
p0.pm.resfac=1e-2; p0.sol.ds=0.05; 
%% centered, 
p=gentau(p0,[1],'d1e');  p.sol.ds=0.1; p=cont(p,5); cplot(p); p=pmcont(p,20); cplot(p);
%% cont in a=base-radius, seems to give closed loop
p=swiparf('d1b','pt20','d1b-20a',4); p.usrlam=[0.25 0.5]; p.sol.ds=-0.01; getaux(p)', 
p.nc.ntot=100; clf(2); p.nc.tol=1e-6; p.sw.bifcheck=0; p.sw.foldcheck=0; 
p.nc.lammax=2; p.nc.lammin=0.05; p.nc.dsmax=0.4; p=cont(p,40); 
%%
fnr=3; figure(fnr); clf; pcmp=0; plotbra('d1b-20a','pt39',fnr,pcmp,'cl','m','lab',[5 10 20 40 35]);
%% solution plots 
p=loadp('d1b-20a','pt5'); plotsol(p); cplot(p); pause 
p=loadp('d1b-20a','pt10'); plotsol(p); cplot(p); pause 
p=loadp('d1b-20a','pt20'); plotsol(p); cplot(p); pause 
p=loadp('d1b-20a','pt30'); plotsol(p); cplot(p); pause 
p=loadp('d1b-20a','pt35'); plotsol(p); cplot(p); pause; 
%%
p=loadp('d1b-20a','pt35'); plotsol(p); cplot(p); pause; 
p=loadp('d1b-20a','pt39'); plotsol(p); cplot(p); 
%% 
p=swiparf('d1d','pt20','d1d-20a',4); p.usrlam=[0.25 0.5]; p.sol.ds=0.01; getaux(p)', 
p.nc.ntot=100; clf(2); p.nc.tol=1e-6; 
p.nc.lammax=2; p.nc.lammin=0.25; p.nc.dsmax=0.5; p=pmcont(p,40); 
figure(3); clf; plotbra('a1a-10a','pt20'); 