function slsolplotT(p,v) % just to make the script shorter
sol=p.cp; dia=sldiagnT(p,15); s1=p.ocopt.s1;
psol3DT(s1,sol,1,1,v,[]); zlabel('P'); xlabel('x'); % plot P 
psol3DT(s1,sol,2,0,v,[]); zlabel('k'); xlabel('x');% plot k 
