function Gu=getGupde(p,u,r)
% GETGUPDE: get jacobian for pde part used in full jacobian of getGu.
%
%  Gu=getGupde(p,u,r)
% Here "r" is the residual at u used in numjac
% Imporant switches: p.sw.jac, p.sw.sfem
%
% See also getGu, filltrafo, stanparam, numjac
global pj; % numjacfac numjacG; % for call of numjac
if (p.sw.jac==1); % pde-part anal. 
    if p.sw.sfem==0 % use full jac
    bc=p.fuha.bcjac(p,u); upde=p.mat.fill*u(1:p.nu); 
    [cj,aj,bj]=p.fuha.Gjac(p,u); zerov=zeros(p.nc.neq,1); 
    [Gu,dum]=assempde(bc,p.mesh.p,p.mesh.e,p.mesh.t,cj,aj,zerov, upde); 
    if(any(bj)); Kadv=assemadv(p.mesh.p,p.mesh.t,bj); 
        Gu=Gu-Kadv; end 
    Gu=filltrafo(p,Gu);
    return;
    else % use "simple jac". Note that p.mat.fill is built in here already
    Gu=p.fuha.sGjac(p,u); 
    end
else
    if p.S==0; % first call, generate sparsity S 
    thresh=p.nc.del*ones(p.nu,1); pj=p; pj.u=u;
    M=getM(p); Ms=p.mat.M1; 
    np=p.nu/p.nc.neq; 
    M(p.nus+1:p.nu,1:p.nus)=ones(p.nu-p.nus,p.nus); % infl. of u(surf) on w (bulk), brute 
    %M(p.nus+1:2*p.nus+2,1:p.nus)=ones(p.nus+2,p.nus); % finer 
    M(1:p.nus,p.nus+1:2*p.nus+2)=ones(p.nus,p.nus+2); 
   % try spmat=p.mat.spmat; catch spmat=ones(p.nc.neq); end;   M=kron(spmat,M(1:np,1:np));  
    S=M>0; u=u(1:p.nu); S(end,:)=ones(1,p.nu); 
    [Gu,njfac,njG]=numjac('resinj',0,u(:),r(1:p.nu),thresh,[],0,S,[]); 
    p.S=abs(Gu)>0; 
    else
       [Gu,njfac,njG]=numjac('resinj',0,u(:),r(1:p.nu),thresh,[],0,p.S,[]); 
    end
  %  figure(11); spy(p.S); figure(12); spy(Gu); pause
end