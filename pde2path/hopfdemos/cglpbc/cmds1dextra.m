%% bif to SWs/TWs at HBP3; 
para=4; ds=0.1; figure(2); clf; aux=[]; aux.tl=30; aux.dlam=0; dir='01D'; hp='hpt3'; nsteps=20;
aux.z=[-1i 1]; ndir='1dtw3'; % gives rTW, speed=om/2
p=hoswibra(dir,hp,ds,para,ndir,aux); 
p.nc.dsmax=0.2; p.hopf.bisec=5; 
p.sw.bifcheck=0; pause; p=cont(p,nsteps); % switch off bif-detection and cont 
%% TWswibra, speed s=om/k, HBP3  
aux.z=[-1i 1]; kwnr=2; p=twswibra(dir,'hpt3',spar,kwnr,'1dtw3b',aux); 
p.u0(1:p.nu)=p.tau(1:p.nu); p.u0=p.u0'; 
p.u0x=p.mat.Kx*p.u0; p.u(1:p.nu)=p.u(1:p.nu)+0.05*p.tau(1:p.nu); 
p.nc.nq=1; p.nc.ilam=[1;6];  % 1 phase-cond, speed as second parameter
p.fuha.qf=@qf; p.sw.qjac=1; p.fuha.qfder=@qjac; p.sw.verb=3; p.nc.tol=1e-6; 
p.sw.bprint=2; clf(2); p.nc.dsmax=0.1; p.sol.ds=0.1; p=cont(p,20);
