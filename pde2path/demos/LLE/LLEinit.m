function p=LLEinit(p,dom,nx,par,varargin)
% initialization
p = stanparam(p); screenlayout(p); 
p.nc.neq = 2; % number of equations in the model
p.sw.sfem = -1; % type of numerical calculation, here OOPDE
p.sw.spjac=1; % use analytical Jacobian for spectral point cont (fold cont)
p.plot.auxdict={'f','\zeta','d','||u_1||_{\infty}','min(|u_1|)'};
%description of the model 
p.fuha.sG=@sG; % the model itself
p.fuha.sGjac=@sGjac; % the Jacobian of the model
p.fuha.spjac=@spjac; % Jacobian for spectral point

switch length(dom)      % control mesh on the right dimension
  case 1; lx=dom; p.pdeo=stanpdeo1D(lx,2*lx/nx); p.vol=2*lx; 
  case 2; ny=round(dom(2)/dom(1)*nx); lx=dom(1); ly=dom(2); p.vol=4*lx*ly; 
      pde=stanpdeo2D(lx,ly,nx,ny,varargin{1}); p.pdeo=pde; p.plot.pstyle=2; 
  case 3; lx=dom(1); ly=dom(2); lz=dom(3); h=2*lx/(nx-1);  p.vol=8*lx*ly*lz; 
      pde=stanpdeo3D(lx,ly,lz,h,varargin{1}); p.pdeo=pde; p.plot.pstyle=2; 
      p.plot.EdgeColor='none'; 
end

p.np = p.pdeo.grid.nPoints; %number grid points of each component
p.nu = p.np*p.nc.neq;         %total grid points
p.nc.neig = 20; p.nc.nsteps = 50; 
p = setfemops(p); p.sol.xi = 1/p.nu; p.file.smod = 10; p.sw.para = 2; p.sw.foldcheck = 1; 
p.nc.ilam=1; p.sol.ds = -0.1; p.nc.dsmax = 0.1; p.nc.lammin = 2; p.sw.bifcheck = 2; 
lam = par(1); u = lam*ones(p.np,1); v = (1/lam)*ones(p.np,1); 
p.u=[u;v;par']; % initial solution guess with parameters
[po,tr,ed]=getpte(p); p.mesh.bp=po; p.mesh.be=ed; p.mesh.bt=tr; p.nc.ngen=1; 
p.plot.bpcmp=4; p.plot.axis='image'; p.plot.cm='hot'; 
p.nc.resfac=1e-3; p.pm.mst=4; 