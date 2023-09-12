function r=sG(p,u)
    u1 = u(1:p.np);         % solution component 1
    u2 = u(p.np+1:2*p.np);  % solution component 2
    par= u(p.nu+1:end);      % parameters
    beta = par(3);
    if p.case == 1  %F is fixed
        alpha = par(1); F = par(2); 
    else            %alpha is fixed
        alpha = par(2); F = par(1); 
    end
    % Orginal function G of equations
    f1 = -u2 - alpha*u1 + (u1.^2 + u2.^2).*u1;
    f2 = u1 - alpha*u2 + (u1.^2 + u2.^2).*u2 - F;
    f = [f1;f2];
    K=kron([[-beta/2,0];[0,-beta/2]],p.mat.K); % assemble full FEM matrix
    r=K*[u1;u2]-p.mat.M*f; % the residual
end
