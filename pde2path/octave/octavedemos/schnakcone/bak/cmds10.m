close all; keep pphome; 
%% Schnakenberg sphere, R=10
R=10; dir='tr10'; nx=41; ny=31; dsmin=0.01; dsmax=0.021; 
p=[]; par=[3.225 0 60 R 0]; lx=pi; del=1e-3; ly=pi/2-del; ref=1; ye=1.2; 
p=schnakSinit(p,[lx,ly],nx,par,ref,ye); p.np, p.pm.resfac=1e-4; p.sol.ds=-dsmin; p.sw.verb=2; 
p.nc.dsmax=dsmax; p.nc.dsmin=dsmin; p=setfn(p,dir); p.nc.mu2=0.01; p.nc.tol=1e-6; p0=p; 
%%
p=cont(p,30);
%% BP1, l=6, using qswibra, first with aux.besw=0 to inspect kernel, then with
% appropriate aux.ali to reduce 'active' kernel to l+1 vectors and then solve 
bp='bpt1'; aux=[]; aux.m=13; aux.besw=0; aux.isotol=1e-8; % small isotol due to high-dim kernel! 
aux.ali=[1 3 5 7 9 11 13]; aux.besw=1; % comment out this line for kernel inspection
p0=qswibra(dir,bp,aux); p0.nc.dsmin=0.01; p0.sw.bifcheck=0; p0.nc.dsmax=0.051; 
p0.pm.resfac=1e-2; p0.sol.ds=0.05; 
%%
p=seltau(p0,1,'a1b',2); p.sol.ds=-0.05; p=cont(p,3); spplot(p); pause; 
p.nc.nq=1;p.nc.ilam=[1 5]; p=pmcont(p,20); spplot(p);
%% ico, short stable part after fold 
p=seltau(p0,2,'a2a',2); p=cont(p,3); spplot(p); pause; 
p.nc.nq=1;p.nc.ilam=[1 5]; p=pmcont(p,20); spplot(p);
%%
p=seltau(p0,2,'a2b',2); p.sol.ds=-0.05; p=cont(p,5); spplot(p); pause; 
p.nc.nq=1; p.nc.ilam=[1 5]; p=pmcont(p,20); spplot(p);
%% stripes
%p=gentau(p0,[0 0 0 0 0 0 0 0 0 0 0 0 1],'10-13a');
p=gentau(p0,[0 0 0 0 0 0 1],'a13a'); % if qswibra with ali
p1=p; pause; p=pmcont(p,15); spplot(p);
%%
p=p1; p.sol.ds=-0.05; pause; p=pmcont(p,15); spplot(p);
%% octa, pitch
%p=gentau(p0,[0 0 0 0 0 0 0 0 0 0 0 0 1],'10-13a');
p=gentau(p0,[0 0 0 1],'a3a');  p=cont(p,4); spplot(p); pause; p.nc.tol=1e-5;
p.nc.nq=1; p.nc.ilam=[1 5]; p.sol.ds=-0.05; p=pmcont(p,20); spplot(p);
%%
p=gentau(p0,[0 0 0 1],'a3b'); p.sol.ds=-0.05; p1=p; p=cont(p,3); spplot(p); pause; 
p.nc.nq=1; p.nc.ilam=[1 5]; p=pmcont(p,10); spplot(p);
%% BP2, l=5; using cswibra with besw=0 to just compute many 'almost' kernel vectors ...
bp='bpt2'; aux=[]; aux.m=13; aux.isotol=1e-8; aux.besw=0; 
aux.ali=[1 2 4 6 8 10]; aux.besw=1; % l=5
p0=cswibra(dir,bp,aux); p0.nc.dsmin=0.01; p0.sw.bifcheck=0; p0.nc.dsmax=0.051; 
p0.pm.resfac=1e-2; p0.sol.ds=0.025; 
%%
p=seltau(p0,3,'b1',3); pause; p.sol.ds=-0.05; p=cont(p,3); spplot(p); pause; 
p.nc.nq=1;p.nc.ilam=[1 5]; p=pmcont(p,10); spplot(p);
%%
p=gentau(p0,1,'b2'); pause; p=pmcont(p,20); spplot(p); 
%% plot BD 
fnr=3; figure(fnr); clf; pcmp=6; plotbra(dir,'pt30',fnr,pcmp,'cl','k','lsw',0);
plotbra('a1a','pt19',fnr,pcmp,'cl',p2pc('r1'),'lab',10); 
plotbra('a1b',fnr,pcmp,'cl',p2pc('r1'),'lab',10); 
plotbra('a2a',fnr,pcmp,'cl',p2pc('b2'),'lab',19); 
plotbra('a2b',fnr,pcmp,'cl',p2pc('b2'),'lsw',0); 
plotbra('a3a','pt20',fnr,pcmp,'cl','m','lab',11); 
plotbra('a13a',fnr,pcmp,'cl',p2pc('o2'),'lab',8); 
plotbra('a13b',fnr,pcmp,'cl',p2pc('o2'),'lab',13); 
plotbra('b1',fnr,pcmp,'cl',p2pc('g1'),'lab',6); 
plotbra('b2','pt10',fnr,pcmp,'cl',p2pc('g2'),'lab',7); 
axis([2.98 3.23 3.01 4.8]); 
ylabel('max u_1'); set(gca,'YTick',[3.2 3.6 4 4.4]);
%% solution plots 
spplotf('a1a','pt10'); pause; spplotf('a1b','pt10'); pause; 
spplotf('a2a','pt19'); pause; spplotf('a3a','pt11'); pause; 
spplotf('a13a','pt8'); pause; spplotf('a13b','pt13'); pause
spplotf('b1','pt6'); pause; spplotf('b2','pt7'); 
%% cont in R, penta-branch, increasing R
p=swiparf('a2a','pt19','a2aRb',[4 5]); huclean(p); p.sw.bifcheck=0; 
p.sol.ds=0.1; p.nc.mu1=1e-2; p.pm.resfac=1e-4;   
p=pmcont(p,10); p.nc.tol=1e-4; p=pmcont(p,10);
%% cont in R, penta-branch, decreasing R 
p=swiparf('a2a','pt19','a2aRa',[4 5]); huclean(p); p.sw.bifcheck=0; 
p.sol.ds=-0.1; p.nc.tol=1e-5; p.pm.resfac=1e-3; p=cont(p,20);
%% cont in R, stripes 
p=swiparf('a13a','pt8','a13aRb',4); huclean(p); p.sw.bifcheck=2; 
p.sol.ds=0.1; p.nc.tol=1e-4; p.pm.resfac=1e-4;  p=pmcont(p,30);
%%
p=swiparf('a13a','pt8','a13aRa',4); huclean(p); p.sw.bifcheck=2; 
p.sol.ds=-0.1; p.nc.tol=1e-6; p.pm.resfac=1e-4;   p=pmcont(p,20);
%% secondary
p=swibra('a13aRa','bpt1','scnd1',-0.01); p.nc.tol=1e-5; p=pmcont(p,20);
%%
fnr=3; figure(fnr); clf; pcmp=6; plotbra('a2aRa','pt20',fnr,pcmp,'cl','b','lab',20);
plotbra('a2aRb','pt50',fnr,pcmp,'cl','b','lab',50);
plotbra('a13aRa','pt20',fnr,pcmp,'cl',p2pc('o2'),'lab',20);
plotbra('a13aRb','pt30',fnr,pcmp,'cl',p2pc('o2'),'lab',30);
plotbra('scnd1','pt20',fnr,pcmp,'cl',p2pc('o3'),'lp',22,'lab',20);
ylabel('max u_1'); 
%%
spplotf('a2aRa','pt20'); pause; spplotf('a2aRb','pt50'); pause 
spplotf('a13aRa','pt20'); pause; spplotf('a13aRb','pt30'); pause;
spplotf('scnd1','pt20');
