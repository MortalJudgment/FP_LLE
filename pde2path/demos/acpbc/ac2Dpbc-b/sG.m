function r=sG(p,u)  % AC with periodic BC 
par=u(p.nu+1:end); sx=par(4); sy=par(5); 
up=u(1:p.nu); % params, and u on periodic domain 
u=p.mat.fill*up; % extend ('fill') u to full domain 
f=par(2)*u+u.^3-par(3)*u.^5; % nonlinearity on extended domain 
F=p.mat.M0*f; % map back to active nodes of periodic domain 
r=par(1)*p.mat.K*up-F-sx*p.mat.Kx*up-sy*p.mat.Ky*up; 