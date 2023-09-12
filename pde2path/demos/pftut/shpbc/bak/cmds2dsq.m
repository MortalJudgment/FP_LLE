%% commands for SH 2nd order sys on rect. for hexagons
keep pphome; close all; p=[]; 
%% init and zero-branch 
p=[]; lx=2*pi; nx=round(8*lx); ly=lx; dir='sq'; ndim=2; lam=-0.001; nu=0; 
par=[lam; nu; 0; 0];  % s1,s2 for PC 
sw.ref=0; sw.sym=2; p=shinit(p,nx,lx,ly,ndim,par,sw); huclean(p); p=setfn(p,dir); 
p.sol.ds=0.005; p.nc.dsmin=0.005; p.sol.dsmax=0.02; 
p.file.smod=10; p.sw.bifcheck=2; p=cont(p,25); 
%%  BP1, stripes and spots 
aux=[];  aux.soltol=1e-8; aux.isotol=1e-3; 
aux.m=4;  aux.besw=0; 
aux.ali=[1 4]; aux.besw=1; % 'active' list
p0=cswibra(dir,'bpt1',aux); p0.sw.bifcheck=2; p0.nc.tol=1e-6; p0.nc.dsmax=0.11; 
%% stripes 
p=seltau(p0,1,'sq1st');p.sol.ds=0.05; p=cont(p,3); p=qxon(p); p=cont(p,20); 
%% spots
p=seltau(p0,3,'sq1sp'); p.sol.ds=0.1; p=cont(p,5); p=qxyon(p); p=pmcont(p,20); 
%%  BP2, kernel is 8 dim, => stripes, rhombs 1,2, simpl.sq, super-sq, anti-sq 
aux.m=4;  aux.besw=0; aux.ali=[]; aux.isotol=1e-3; 
aux.ali=[1 4]; aux.besw=1; aux.ral=1; 
p0=cswibra(dir,'bpt2',aux); p0.sw.bifcheck=0; p0.nc.tol=1e-5;  p0.nc.dsmax=0.11; 
%% (obl.) stripes 
p=seltau(p0,1,'sq2st',3); p.sol.ds=0.02; p=cont(p,3); p=qxon(p); p=pmcont(p,25); 
%% spots
p=gentau(p0,1,'sq2sp'); p.sol.ds=0.1; p=cont(p,5); p=qxyon(p);   p=pmcont(p,25); 
%%  BP3, diagonal spot and stripes 
aux=[]; aux.besw=0; 
aux.isotol=1e-3; aux.ali=[1 4]; aux.besw=1; aux.ral=1; 
p0=cswibra(dir,'bpt3',aux);  p0.nc.dsmax=0.11; 
%%
p=seltau(p0,3,'sq3sp'); p.sol.ds=0.05; p=cont(p,2); p=qxyon(p); p.nc.ilam=[1 3 4]; p=cont(p,10); 
p=seltau(p0,1,'sq3st',3); p.sol.ds=0.05; p=cont(p,2); p=qxon(p); p.nc.ilam=[1 3]; p=cont(p,10); 
%% BD-plots 
fnr=3; figure(fnr); cmp=5; clf; 
plotbra('sq',fnr,cmp,'cl','k','lsw',0); plotbra('sq1st',fnr,cmp,'cl','b'); 
plotbra('sq1sp',fnr,cmp,'cl','r'); plotbra('sq2s',fnr,cmp,'cl','b');
plotbra('sq2sp',fnr,cmp,'cl',p2pc('r1')); plotbra('sq3st',fnr,cmp,'cl','b'); 
plotbra('sq3sp',fnr,cmp,'cl',p2pc('r1')); 
title('\nu=0'); xlabel('\lambda'), ylabel('max(|u|)'); set(gca,'XTick',[0 0.5 1]); 
%% soln plots 
plotsol('sq1sp','pt18'); pause; plotsol('sq1s','pt15'); pause; 
plotsol('sq2s','pt21'); pause; plotsol('sq2sp','pt15'); 
pause; plotsol('sq3s','pt10'); 
%% with fourier 
d='sq2sp'; pt='pt10'; p=loadp(d,pt); plotsol(p); fouplot(p,10,1,[3 5],[d '/' pt],1); 