function p=LLEinit(p, nnodes, rho, lx, ds, par)
%% setting standard parameters
    p=stanparam(p); % infuses p with standard parameter settings
    screenlayout(p); % open, clear and arrange the common figures

%% special parameters related to this model
    % basics
    p.nc.neq=2; % number of equations in the model
    p.sw.sfem=-1; % type of numerical calculation, here OOPDE
    p.sw.spjac=1; % use analytical Jacobian for spectral point cont (fold cont)
    % names for cmp of stanbra
    p.plot.auxdict={'F','\alpha','\beta','L2-norm'};
    % description of the model
    p.fuha.sG=@sG; % the model itself
    p.fuha.sGjac=@sGjac; % the Jacobian of the model
    p.fuha.spjac=@spjac; % Jacobian for spectral point cont (fold cont)

%% domain and mesh
    p.vol=2*lx;
    p.pdeo = stanpdeo1D(lx,2*lx/nnodes);
    p.np=p.pdeo.grid.nPoints; % number of meshpoints
    p.nu=p.np*p.nc.neq; % number of unknowns (=2*(mesh points), as 2 components)
    p=setfemops(p);     % compute FEM-operators
%% bifurcation parameter, continuation basics and first guess for solution
    p.nc.ilam=1; % primary bifurcation parameter located at p.u(p.np+p.nc.ilam)
    p.sol.xi=1/p.nu; % weight in arclength-continuation
    p.sol.ds=ds; % starting stepsize
    p.nc.dsmax=10*ds; % maximal stepsize
    p.nc.dsmin=0; % minimal stepsize
    % construction the trivial solution
    F=par(2); alpha=par(1);
    psif=F/(1+1i*(alpha-rho));
    ov=ones(p.np,1); % dummy for the 1 function
    u1=real(psif)*ov; % initial guess for u_1
    u2=imag(psif)*ov;  % initial guess for u_2
    p.u=[u1;u2;par']; % initial solution guess with parameters
end