function Gu=sGjac(p,u)
% Defines the analytical expression of the Jacobian,
J=p.np-1;
h=p.vol/J;
par = u(p.nu +1: end) ; % parameters par=[F,alpha,beta]
lletype = p.type;
nno = p.np ; % nb of nodes in the FEM mesh

[f1u,f1v,f2u,f2v]=njac(p,u,par); % the Jacobian 

Fu=[[spdiags(f1u,0,nno,nno),spdiags(f1v,0,nno,nno)];
     [spdiags(f2u,0,nno,nno),spdiags(f2v,0,nno,nno)]];
 
beta=par(3);
D=[[beta/2,0];[0,beta/2]]; 
if lletype == 0
    Gu = kron(D,p.mat.K) + p.mat.M * Fu ; % assemble the Jacobian
else
    Qu = nlterm(p,u);
    Gu = kron(D,p.mat.K) + p.mat.M * Fu + lletype*h/pi*Qu; % assemble the Jacobian
end
end

function [f1u,f1v,f2u,f2v]= njac(p,u,par)
% Jacobian for LLE

nno = p.np ; % nb of nodes in the FEM mesh
u1 = u (1:nno) ; % solution component 1
u2 = u(nno +1:2*nno) ; % solution component 2

ov = ones(nno,1); % dummy for the 1 function
alpha = par(2);
% entries of the jacobian
f1u = -alpha*ov + 3*u1.^2+u2.^2  ;
f1v = -ov + 2*u1.*u2;

f2u = ov + 2*u1.*u2;
f2v = -alpha*ov + 3*u2.^2+u1.^2;

end

function Qu = nlterm(p,u)
% nonlinear term

nno = p.np ; % nb of nodes in the FEM mesh
u1 = u(1:nno) ; % solution component 1
u2 = u(nno +1:2*nno) ; % solution component 2

In = eye(2*nno); 

NN1 = 2*(u(1:2*nno)'*u(1:2*nno))*In + 4*u(1:2*nno)*u(1:2*nno)';
T = zeros(2*nno,2*nno);
T(:,1) = u1(1)*u(1:2*nno);
T(:,nno) = u1(end)*u(1:2*nno);
T(:,nno+1) = u2(1)*u(1:2*nno);
T(:,2*nno) = u2(end)*u(1:2*nno);
NN2 = (u1(1)^2 + u1(end)^2 + u2(1)^2 + u2(end)^2)*In + 2*T;
Qu = NN1 + NN2;
end

