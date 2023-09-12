function Gu=sGjac(p,u)
    par=u(p.nu+1:end); % parameters
    n=p.np;
    [f1u,f1v,f2u,f2v]=njac(p,u,par); % the Jacobian, see below
    Fu=[[spdiags(f1u,0,n,n),spdiags(f1v,0,n,n)];
       [spdiags(f2u,0,n,n),spdiags(f2v,0,n,n)]];
%     Gu=kron([[par(3)/2,0];[0,par(3)/2]],p.mat.K) - p.mat.M*Fu; % assemble the Jacobian
    Gu=kron([[-par(3)/2,0];[0,-par(3)/2]],p.mat.K) - Fu; % assemble the Jacobian
end 
function [f1u,f1v,f2u,f2v]=njac(p,u,par) % Jacobian for Schnakenberg
    u1=u(1:p.np);           % solution component 1
    u2=u(p.np+1:2*p.np);    % solution component 2
    n = p.np; % number of function points per component
    ov = ones(n,1); % dummy for the 1 function
% % entries of the jacobian
    f1u = -par(2)*ov + 3*u1.^2 + u2.^2;
    f1v = -ov + 2*u1.*u2;
    f2u = ov + 2*u1.*u2;
    f2v = -par(2)*ov + u1.^2 + 3*u2.^2;
end
