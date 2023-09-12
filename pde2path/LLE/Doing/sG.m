function r=sG(p,u)
    u1 = u(1:p.np);         % solution component 1
    u2 = u(p.np+1:2*p.np);  % solution component 2
    par= u(p.nu+1:end);      % parameters
    % Orginal function G of equations
    f1 = -u2 - par(2)*u1 + (u1.^2 + u2.^2).*u1;
    f2 = u1 - par(2)*u2 + (u1.^2 + u2.^2).*u2 - par(1);
    f = [f1;f2];
    K = kron([[-par(3)/2,0];[0,-par(3)/2]],p.mat.K) ; % assemble full FEM matrix
%     r = K*[u1;u2] - p.mat.M*f; 
    r = K*[u1;u2] - f; 
end