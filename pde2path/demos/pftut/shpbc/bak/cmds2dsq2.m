%% commands for SH 2nd order sys on rect. for hexagons
keep pphome; close all; p=[]; 
%% init and zero-branch 
p=[]; lx=2*pi; nx=round(6*lx); ly=lx; dir='sq2'; ndim=2; lam=-0.001; nu=0; 
par=[lam; nu; 0; 0];  % s1,s2 for PC 
sw.ref=0; sw.sym=2; p=shinit(p,nx,lx,ly,ndim,par,sw); huclean(p); p=setfn(p,dir); 
p.sol.ds=0.005; p.nc.dsmin=0.005; p.sol.dsmax=0.02; 
p.file.smod=10; p.sw.bifcheck=2; p=cont(p,25); 
%%  BP1, stripes and spots 
aux=[]; aux.soltol=1e-16; aux.isotol=1e-12; aux.m=4;  aux.besw=0; 
aux.ali=[1 2]; aux.besw=1; %aux.ral=1;  % 'active' list
p0=cswibra(dir,'bpt1',aux); p0.sw.bifcheck=2; p0.nc.tol=1e-6; p0.nc.dsmax=0.11; 
%% stripes 
p=seltau(p0,1,'sqb1s',3); p.sol.ds=0.05; p=cont(p,3); p=qyon(p); p=cont(p,20); 
%% spots
p=seltau(p0,2,'sqb1sp'); p.sol.ds=0.1; p=cont(p,5); p=qxyon(p); p=pmcont(p,20); 
%%  BP2, oblique stripes and spots 
aux.m=8;  aux.besw=0; aux.ali=[]; aux.isotol=1e-18; 
aux.ali=[1 2 4 6]; aux.besw=1; 
p0=cswibra(dir,'bpt2',aux); %p0.sw.bifcheck=0; 
p0.nc.tol=1e-6;  p0.nc.dsmax=0.11; 
%% (obl.) stripes 
p=seltau(p0,1,'sqb2s',3); p.sol.ds=0.02; p=cont(p,5); p=qyon(p); p=cont(p,25); 
%% spots
p=seltau(p0,8,'sqb2sp',3); p.sol.ds=0.02; p=cont(p,5); p=qxyon(p); p=cont(p,25); 
%%  BP3, vert. stripes 
aux=[]; aux.m=4;  aux.besw=0; aux.ali=[]; aux.isotol=1e-18; 
aux.ali=[1 4]; aux.besw=1; 
p0=cswibra(dir,'bpt3',aux);  p0.nc.dsmax=0.11; 
%% (obl.) stripes 
p=seltau(p0,1,'sqb3s',3); p.sol.ds=0.02; p=cont(p,5); p=qyon(p); p=cont(p,25); 
%% spots
p=seltau(p0,2,'sqb3sp',3); p.sol.ds=0.02; p=cont(p,5); p=qxyon(p); p=cont(p,25); 
%%
p=gentau(p0,1,'sqb3s'); p.sol.ds=0.05; p=cont(p,2); p=qxon(p); p.nc.ilam=[1 3]; p=cont(p,10); 
%% BD-plots 
fnr=3; figure(fnr); cmp=5; clf; 
plotbra('sqb',fnr,cmp,'cl','k','lsw',0); plotbra('sqb1s',fnr,cmp,'cl','b'); 
plotbra('sqb1sp',fnr,cmp,'cl','r'); plotbra('sqb2s',fnr,cmp,'cl',p2pc('b2'));
plotbra('sqb2sp',fnr,cmp,'cl',p2pc('r1')); plotbra('sqb3s',fnr,cmp,'cl',p2pc('b3')); 
 plotbra('sqb3sp',fnr,cmp,'cl',p2pc('o1')); 
axis([0 1.05 0 0.9]); box on;
title('\nu=0'); xlabel('\lambda'), ylabel('max(|u|)'); set(gca,'XTick',[0 0.5 1]); 
%% soln plots 
plotsol('sq1sp','pt18'); pause; plotsol('sq1s','pt15'); pause; 
plotsol('sq2s','pt22'); pause; plotsol('sq2sp','pt15'); 
pause; plotsol('sq3s','pt10'); 
%% with fourier 
d='sq2sp'; pt='pt10'; p=loadp(d,pt); plotsol(p); fouplot(p,10,1,[3 5],[d '/' pt],1); 