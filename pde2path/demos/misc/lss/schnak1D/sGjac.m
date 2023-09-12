function Gu=sGjac(p,u)
    par=u(p.nu+1:end); % identify parameters
    n=p.np;
    [f1u,f1v,f2u,f2v]=njac(p,u,par); % loading the non-linear jacobian in
                                     % nodal form, see below
    Fu=[[spdiags(f1u,0,n,n),spdiags(f1v,0,n,n)];
       [spdiags(f2u,0,n,n),spdiags(f2v,0,n,n)]]; % nonlinear part of the jacobian
    Gu=kron([[1,0];[0,par(3)]],p.mat.K)-p.mat.M*Fu; % assemble the jacobian
end 
function [f1u,f1v,f2u,f2v]=njac(p,u,par) % Jacobian for Schnakenberg, nodal
                                         % version
    u1=u(1:p.np); % identify solution component 1
    u2=u(p.np+1:2*p.np); % identify solution component 2
    % entries of the jacobian
    f1u=-1+2*u1.*u2+2*par(2)*(u1-u2.^(-1));
    f1v=u1.^2+2*par(2)*u2.^(-2).*(u1-u2.^(-1)); 
    f2u=-2*u1.*u2-2*par(2)*(u1-u2.^(-1)); 
    f2v=-u1.^2-2*par(2)*u2.^(-2).*(u1-u2.^(-1)); 
end