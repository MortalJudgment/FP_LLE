%% 1 - creating a clear working space
close all ; keep pphome ;
%% 2 - initialising the problem
nbstepstb=300; % nb of steps along the trivial branch
nbstepsbb=200; % nb of steps along the bifurcation branches
ds=0.01; % starting stepsize in the continuation method
nnodes=250; % nb of nodes for the FEM discretization
dirname='test';
% Parameters of the LLE
Fcste=8;
alpha=10;
beta=0.4;
lx = pi/2 ; % set domain length [ - lx , lx ]  
% Compute flat solutions of the LLE 
lletype = 0; % 0 - original LLE, 1 - LLE for Fabry-Perot resonators 
polyn3=[1+2*lletype,-2*alpha*(1 + 2*lletype),(1+alpha^2),-Fcste^2];
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
par =[Fcste, alpha, beta]; 
p = schnakinit(p, nnodes, rhof, lx, ds, par) ;
p.fuha.outfu=@llebra; % add L^2 norm of the solution to the output in p.branch
p.plot.pmod =10; % shows each 10 th solution in fig 2 only
p.file.smod =10; % stores each 10 th solution only
p = setfn(p,dirname) ;  % set output directory to the name indicated in variable dirname
p.sw.foldcheck =1; % enables/disable detection of folds
p.sw.bifcheck =1; % enables/disable detection of bifs

p.plot.bpcmp=4; % plot the 4th user defined component of p.branch 
                % (the L^2 norm of the complex solution here as defined in
                % llebra.m)
                % Note that the default value is 0 and the L^2 norm of the
                % 1st compoment in then plot

%% 3 - contiunation of the trivial branch
p = cont(p ,nbstepstb) ; % continuation for a maximum of 300 steps 
                  % with stepsize p.sol.ds, starting from p.u, direction tau
                
                
                
%% 4 - contiunation of the branches that bifurcate from the trivial branch
p=swibra(dirname,'bpt1','branch1'); % switch to branch
p=cont(p,120); % cont for max 200 steps
%% 5 - continue on another branch from trival branch
% p=swibra(dirname,'bpt2','branch2'); % switch to branch
% p=cont(p,200); % cont for max 200 steps
p=swibra(dirname,'bpt3','branch3'); % switch to branch
p=cont(p,150); % cont for max 200 steps
% p=swibra(dirname,'bpt5','branch5'); % switch to branch
% p=cont(p,200); % cont for max 200 steps
p=swibra(dirname,'bpt6','branch6'); % switch to branch
p=cont(p,100); % cont for max 200 steps
% p=swibra(dirname,'bpt7','branch7'); % switch to branch
% p=cont(p,200); % cont for max 200 steps
p=swibra(dirname,'bpt12','branch12'); % switch to branch
p=cont(p,500); % cont for max 200 steps

figure(3) ;
clf ;
plotbra(dirname) ; 
plotbra('branch1','cl','r');
% plotbra('branch2','cl',[0.6350 0.0780 0.1840]); 
plotbra('branch3','cl','b');
% plotbra('branch5','cl',[0.8500 0.3250 0.0980]);
plotbra('branch6','cl','g');
% plotbra('branch7','cl','k');
plotbra('branch12','cl','m');
grid
