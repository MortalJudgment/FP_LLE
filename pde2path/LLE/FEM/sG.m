function r=sG(p,u) 
% Computation of the residual of the PDE
% matrices p.mat.K and p.mat.M here which are the FEM and mass matrices
% will be defined in oosetfemops.m

nno=p.np; % nb of nodes in the FEM mesh
nnov=p.nu; % nb of nodal value for the 2 components 

u1 = u(1:nno) ; % solution component 1
u2 = u(nno +1:2*nno) ; % solution component 2
par = u(nnov+1:end) ; % parameters par=[F,alpha,beta]

F=par(1);
alpha=par(2);
beta=par(3);

f1 = -u2 - alpha*u1 + (u1.^2 + u2.^2).*u1;
f2 = u1 - alpha*u2 + (u1.^2 + u2.^2).*u2 - F;
f = [f1;f2];
D=[[-beta/2,0];[0,-beta/2]]; 
K = kron(D, p.mat.K) ; % assemble full FEM matrix
r = K *[ u1 ; u2 ] - p.mat.M * f ; % the residual

end
