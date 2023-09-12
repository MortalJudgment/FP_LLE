function out=chbra(p,u) %  
E=chE(p,u); 
out=[u(p.nu+1:end); E; max(abs(u(1:p.np))); min(abs(u(1:p.np)))]; 