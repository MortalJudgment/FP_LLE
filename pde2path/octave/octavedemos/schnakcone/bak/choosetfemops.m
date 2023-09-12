function p=oosetfemops(p) % set FEM operators 
fem=p.pdeo.fem; gr=p.pdeo.grid; a=p.u(p.nu+3); 
po=getpte(p); x=po(1,:)'; y=po(2,:)'; [K,M]=LBc(p,a); 
Krot=convection(fem,gr,[-y;x]);p.mat.M=M; p.mat.K=K; p.mat.Krot=Krot;
p.mat.vM=sum(M,1); 