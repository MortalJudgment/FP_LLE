function out=llebra(p,u)
% output to bifurcation diagram function :
% L2 norm of the solution u1 + i u2 and equ. parameters
out=[
    u(p.nu+1:end) ; % parameters 
    sqrt(norm(u(1:p.nu))^2*p.vol/p.np); % L2 norm
    ];