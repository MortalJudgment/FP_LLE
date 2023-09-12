%% 1 - creating a clear working space
close all ; keep pphome ;
%% 2 - initialising the problem
dirname='test'; % diectory where the results will be stored
% if the folder exists, it will be removed otherwise the new results are added to the previous ones  
% if isfolder(dirname)
%     disp(['Folder ',dirname,' exists and will be removed'])
%     yn=input('Confirm this action (y/n) : ','s');
%     if ~strcmp(yn,'y')
%         return;
%     else
%         rmdir(dirname, 's')
%     end
% end
% Parameters of the LLE
lletype=0;
Fcste=1.6;
alpha=-0.5;
beta=-0.2;
bparameter = 2; % bifurcation parameter 1 = F, 2 = alpha, 3 = beta
lx = pi ; % set domain length [ - lx/2 , lx/2 ]  
% Main discretization parameters
nbstepstb = 200; % nb of steps along the trivial branch
nbstepsbb = 200; % nb of steps along the bifurcation branches
ds = 5E-3; % starting stepsize in the continuation method
nnodes = 400; % nb of nodes for the FEM discretization
% Compute flat solutions of the LLE 
polyn3=[1+2*lletype,-2*alpha*(1+2*lletype),(1+alpha^2),-Fcste^2];
racines=roots(polyn3);
Trho=[];
for j=1:3
    if isreal(racines(j)) && racines(j)>=0
        Trho=[Trho,racines(j)];
    end
end
Trho=sort(Trho); % square modulus of the flat solutions
Tpsif=Fcste./(1+1i*(alpha-(1+2*lletype)*Trho)); % Flat solutions
disp(['There are ',num2str(length(Trho)),' flat solutions '])
for j=1:length(Trho)
    disp([num2str(j),' - square modulus of the flat sol. : rho = ',num2str(Trho(j))]);
    disp(['   - L2 norm of the flat sol.  = ',num2str(sqrt(2*lx*Trho(j)))]);
end
j=input('Which one to use ? index value = ');
rhof=Trho(j);
% Initialisation of pde2path
p =[];
%par =[Fcste, alpha, beta]; 
par =[Fcste, alpha, beta, lx,]; 
p = lleinit(p, nnodes, rhof, lx, ds, bparameter, par) ;
p.fuha.outfu=@llebra; % add L^2 norm of the solution to the output in p.branch
p.plot.pmod =10; % shows each 10 th solution in fig 2 only
p.file.smod =10; % stores each 10 th solution only
p = setfn(p,dirname) ;  % set output directory to the name indicated in variable dirname
p.sw.foldcheck =1; % enables/disable detection of folds
p.sw.bifcheck =1; % enables/disable detection of bifs
                  % 0/1/2 for bif.detection off/ via LU decomposition (default) /via counting eigenvalues,

p.plot.bpcmp=5; % plot the 4th user defined component of p.branch 
                % (the L^2 norm of the complex solution here as defined in
                % llebra.m)
                % Note that the default value is 0 and the L^2 norm of the
                % 1st compoment in then plot

% check the jacobian
% p=oosetfemops(p);
% [Gua, Gun]=jaccheck(p);


%% 3 - contiunation of the trivial branch
close(6);
p = cont(p ,nbstepstb) ; % continuation for a maximum of 300 steps 
                  % with stepsize p.sol.ds, starting from p.u, direction tau
                                  
nbpt=p.file.bcount-1; % nb of bifurcation points detected
nfpt=p.file.fcount-1; % nb of fold (saddle) points detected
%
disp(['nb of b.p. found = ',num2str(nbpt)]);     
disp('  | j | Lamnda | L^2 norm');
ind=1;
for j=1:size(p.branch,2)
   if p.branch(2,j)==1
       disp([num2str(ind), ' ', num2str(j),' ', num2str(p.branch(4,j)) ' ', num2str(p.branch(11,j))]);
       ind=ind+1;
   end
end
%% 4 - contiunation of the branches that bifurcate from the trivial branch


nbpt=input('How many branches to be investigated ? ans = ');
TJ=input('List of indexes of bifurcation points to investigate : ','s');
TJ=str2num(TJ);
disp(['nb of fold points found = ',num2str(nfpt)]);
nbpt=length(TJ);

nbptbranch=zeros(1,nbpt);
nfptbranch=zeros(1,nbpt);
Tnsteps=nbstepsbb*ones(1,nbpt);
for j=1:nbpt
    bradirname=[dirname,'/branche',num2str(TJ(j))];
    bptname=['bpt',num2str(TJ(j))];
    p = swibra(dirname, bptname, bradirname) ; % switch to new branch
    p = cont(p ,Tnsteps(j)) ; % continuation for a maximum nsteps(j) steps
    nbptbranch(j)=p.file.bcount-1; % nb of bifurcation points detected in the curent branch
    nfptbranch(j)=p.file.fcount-1; % nb of fold points detected in the curent branch
end

for j=1:nbpt
    disp(['On branch #',num2str(TJ(j)),' : ', 'B.P. detected = ',num2str(nbptbranch(j)), ...
        '   F.P. detected = ',num2str( nfptbranch(j))]);
end

%
%     bradirname=[dirname,'/branchefp'];
%     bptname='fpt1';
%     p = swibra(dirname, bptname, bradirname) ; % switch to new branch
%     p = cont(p ,200) ; % continuation for a maximum nsteps(j) steps


%% 5 - plot bifurcation diagram
figure(p.plot.brafig);clf;
varargin=['tyun','--'];
plotbra(dirname) ;

% Array of Matlab colors to plot the various branches in different colors
col=['b','g','r','c','m','y'];
j=floor(nbpt/6)+1;
Tcol=repmat(col,1,j);

for j=1:nbpt
    bradirname=[dirname,'/branche',num2str(TJ(j))];
    plotbra(bradirname,'cl',Tcol(j)) ; %
end
plotbra(dirname) ; % plot trivial branch
grid
