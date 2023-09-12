function p=oosetfemops(p)
J=p.np-1;
h=p.vol/J;
B=[-ones(1,J+1);2*ones(1,J+1);-ones(1,J+1)];
B=transpose(B);
p.mat.K=spdiags(B,[-1,0,1],J+1,J+1);
p.mat.K(1,1)=1 ; p.mat.K(J+1,J+1)=1 ; % Neumann B.C.
p.mat.K=p.mat.K/h^2;
p.mat.M=kron([[1,0];[0,1]],speye(J+1));
end 
