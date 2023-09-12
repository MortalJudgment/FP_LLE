function out=llebra(p,u)
% Compute additionnal output to be used for the bifurcation diagram
% We add to p.branch the equ. parameters and L2 norm of the solution u1 + i
% u2 evaluated by the rectangle quadrature rule
out=[
    u(p.nu+1:end) ; % parameters 
    sqrt(norm(u(1:p.nu))^2*p.vol/p.np); % L2 norm
    ];