function p=lleinit(p,nno,rhof,lx,ds,bparameter,par)
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
%
screenlayout(p) ; % open, clear and arrange the common figures
p.plot.ifig=0; % do not plot tangent vector
p.plot.pfig=1; % figure for the real part of the solution (p.plot.pfig+1 for the imaginary part)
p.plot.brfig=3; % figure for plotting the bifurcation paths during computation progress 
p.plot.brafig=4; % figure for the final plot of the bifurcation diagram
p.plotudict={'Re(\psi)','Im(\psi)'};
p.plot.auxdict={'F','\alpha','\beta','h','L^2 norm'};
p.plot.pstyle = -1; % flag to call userplot in plotsol
% special parameters related to this model
% basics
p.nc.neq=2; % number of equations in the model
p.sw.sfem = -1; % use the sG / sGjac and oosetfemops setting ( without OOPDE !)
p.sw.spjac=1; % (1/0) use analytical/numerical Jacobian for spectral point cont (fold cont) when 1 
p.sw.jac=1; % (1/0) use analytical/numericalJacobian for G
p.sw.qjac=0; % (1/0) use analytical/numerical Jacobian for g

% description of the model
p.fuha.sG=@sG; % the model itself
p.fuha.sGjac=@sGjac; % the Jacobian of the model

% domain and mesh
p.vol=lx;
p.x=linspace(-lx/2,lx/2,nno);
p.np=nno; % set the number of mesh points
p.nu=p.np*p.nc.neq; % set the number of unknowns (=2*(mesh points ), as 2 components ) 
 

% bifurcation parameter , continuation basics and first guess for solution
p.sol.xi=1/p.nu; % weight in arclength - continuation
p.nc.ilam=bparameter;  % primary bifurcation parameter located at p.u(p.nu+p.nc.ilam )
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

