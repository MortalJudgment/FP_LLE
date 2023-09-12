function q=qf(p,u) 
m=u(p.nu+1); u=u(1:p.nu); q=p.mat.vM*u/p.Om-m; end  % mass constraint 