function p=chcinit(p,r1,nr,par,ell) 
p=stanparam(p); screenlayout(p); p.nc.neq=1; p.sw.sfem=-1; 
p.fuha.sG=@sG; p.fuha.sGjac=@sGjac; p.sw.jac=1; p.cf=2; 
%pde=diskpdeo(r1,nr); 
%pde=diskpdeo2(r1,nr,4); 
pde=diskpdeo2b(r1,nr,4,ell); 
p.pdeo=pde; p.np=pde.grid.nPoints; p.nu=p.np; 
p.plot.pstyle=2; p.plot.cm='cool'; p.plot.axis='image'; p.plot.nola=1; 
p.nc.neig=40; p.sol.xi=1/(p.nu);  hom=ones(p.np,1); p.u=[hom; par']; 
p=oosetfemops(p); p.vol=p.mat.vM*ones(p.nu,1); p.nc.ilam=[1 4]; 
p.u=[par(1)*hom; par']; 
p.nc.nsteps=20; p.sw.foldcheck=1; 
p.plot.auxdict={'m','eps','a','s'}; % parameter names 
p.nc.lammax=4; p.sol.ds=0.05; p.nc.dsmax=0.05; p.nc.nsteps=100; 
p.usrlam=[-0.5 0]; plotsol(p); p.plot.bpcmp=6; 
p.nc.nq=1; p.fuha.qf=@qf2; p.fuha.qfder=@qf2der; p.fuha.outfu=@chbra; 
r=resi(p,p.u); res=norm(r,'inf'); fprintf('ini-res=%g\n',res); %p.sw.jac=0; 
[u,res,iter,Gu,Glam,p]=nloop(p,p.u); fprintf('1st-res=%g\n',res); 