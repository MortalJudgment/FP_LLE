function p=lleinit(p,nno,rhof,lx,ds,par)
% This function is used to set the standard parameters of the pde2path
% programm
%
% Parameters in :
% p : pde2path structure
% nno : nb of nodes for the FEM discretization
% rhof : sqaure modulus of the flat solution of the LLE
% lx : the LLE is considered over the interval [-lx,lx]
% ds : starting stepsize in the continuation method
%
% Parameters out :
% All passed through the p struct 

p = stanparam(p) ; % infuses p with standard parameter settings
screenlayout(p) ; % open, clear and arrange the common figures

% special parameters related to this model
% basics
p.nc.neq=2; % number of equations in the model
p.sw.sfem=-1; % type of numerical calculation, here OOPDE
p.sw.spjac=1; % use analytical Jacobian for spectral point cont (fold cont)
p.plotudict={'Re(\psi)','Im(\psi)'};
p.plot.auxdict={'F','\alpha','\beta','||\psi||^2_2'};
% description of the model
p.fuha.sG=@sG; % the model itself
p.fuha.sGjac=@sGjac; % the Jacobian of the model
p.fuha.spjac=@spjac; % Jacobian for spectral point cont ( fold cont )

% domain and mesh
p.vol=2*lx;
p.pdeo=stanpdeo1D(lx,2*lx/nno); % mesh [ - lx , lx ] , max mesh pt 2* lx / nno

p.np=p.pdeo.grid.nPoints; % set the number of mesh points
p.nu=p.np*p.nc.neq; % set the number of unknowns (=2*(mesh points ), as 2 components ) 
p=setfemops(p); % compute FEM - operators

% bifurcation parameter , continuation basics and first guess for solution
p.sol.xi=1/p.nu; % weight in arclength - continuation
p.nc.ilam=1;  % primary bifurcation parameter located at p . u ( p . nu + p . nc . ilam )
p.sol.ds=ds;% starting (current) stepsize
p.nc.dsmax=10*ds;% maximal stepsize
p.nc.dsmin=0;% minimal stepsize

% construction the trivial solution
F=par(1);
alpha=par(2);
psif=F/(1+1i*(alpha-rhof));
nno=p.np; % nb of nodes in the FEM mesh
ov=ones(nno,1); % dummy for the 1 function
u1=real(psif)*ov; % initial guess for u_1
u2=imag(psif)*ov;  % initial guess for u_2
p.u=[u1;u2;par']; % initial solution guess with parameters

