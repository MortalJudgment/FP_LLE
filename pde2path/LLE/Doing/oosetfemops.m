function p=oosetfemops(p)
% [p.mat.K,M,~]=p.pdeo.fem.assema(p.pdeo.grid,1,1,1); % FEM/mass matrices
% % mass matrix adaption for the problem, as it is a system
% p.mat.M=kron([[1,0];[0,1]],M);
n = p.np; % number of function points per component
h = p.vol/(n-1);
e = ones(n,1) ;
A = spdiags( [-e, 2*e, -e], [-1 0 1] , n,n ) ;
% Change point near boundary
A(1,1) = 1; A(2,1) = -1;
A(end,end) = 1; A(end-1,end) = -1;
p.mat.K = 1/h^2*A;

M = eye(n);
p.mat.M = kron([[1,0];[0,1]],M);
end 
