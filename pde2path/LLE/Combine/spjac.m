function Gvvph=spjac(p,u)
u1 = u(1:p.np);         % solution component 1
u2 = u(p.np+1:2*p.np); 	% solution component 2
n = p.np; % number of function points per component
% second order derivations of the model
f1uu=6*u1; f1uv=2*u2; f1vv=2*u1;
f2uu=2*u2; f2uv=2*u1; f2vv=6*u2;
% implementation of the derivations as sparse matrices
ph1=u(p.nu+1:p.nu+p.np);
ph2=u(p.nu+p.np+1:2*p.nu);
M1=spdiags(f1uu.*ph1+f1uv.*ph2,0,n,n); 
M2=spdiags(f1uv.*ph1+f1vv.*ph2,0,n,n); 
M3=spdiags(f2uu.*ph1+f2uv.*ph2,0,n,n); 
M4=spdiags(f2uv.*ph1+f2vv.*ph2,0,n,n);
Gvvph=-p.mat.M*[[M1 M2]; [M3 M4]]; 
end

% function Guuph=spjac(p,u)
% % Defines the analytical Jacobian for spectral continuation used for
% % fold and branch point continuation. 
% % It is not necessary to define it, but, as it has to be calculated 
% % numerically otherwise, speeds up fold continuation a lot.   
% 
% nno=p.np; % number of function points per component = nb of nodes in the FEM mesh
% nnov=p.nu; % nb of nodal value for the 2 components 
% 
% u1=u(1:nno); % first component
% u2=u(nno+1:2*nno); % second component
% 
% % second order derivations of the model
% f1uu=-2*u2; 
% f1uv=-2*u1; 
% f1vv=-6*u2;
% f2uu=6*u1;
% f2uv=2*u2; 
% f2vv=2*u1;
% % implementation of the derivations as sparse matrices
% ph1=u(nnov+1:nnov+nno); ph2=u(nnov+nno+1:2*nnov);
% M1=spdiags(f1uu.*ph1+f1uv.*ph2,0,nno,nno); 
% M2=spdiags(f1uv.*ph1+f1vv.*ph2,0,nno,nno); 
% M3=spdiags(f2uu.*ph1+f2uv.*ph2,0,nno,nno); 
% M4=spdiags(f2uv.*ph1+f2vv.*ph2,0,nno,nno);
% Guuph=-p.mat.M*[[M1 M2]; [M3 M4]]; 
% end