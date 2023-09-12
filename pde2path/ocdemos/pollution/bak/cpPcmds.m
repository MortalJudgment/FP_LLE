
hdir='h1'; hpt='pt8'; fdir='FSS'; fpt='pt9'; 
[muv1, muv2,ind,h]=floqpsap(hdir,hpt); fprintf('d(u_H)=%i\n',ind); 
q=loadp(hdir,hpt); % load hopf orbit
q.hopf.tom.Nmax=1e10; % switch of mesh control of mtom 
huclean(q); hoplot(q,1,1); 
if 0 % refine hopf orbit in time steps for better projection
q=uhopftref(q,1.5); % refine mesh
[y,T,lam,res,iter,A,q]=honloop(q,q.hopf.y,q.hopf.T,q.u(q.nu+q.nc.ilam)); % recompute sol on finer grid
q.hopf.y=y; q.hopf.T=T;
res, hoplot(q,2,1); 
end
%% Cell2: computation the of canonical path from FSS/pt9 to h2/pt17, al=0.95
p=[]; p=ocinit(p,fdir,fpt,hdir,hpt); % load standard options
p.oc.s1=q; % overwrite starting states with finer CPS from Cell 1
p.oc.nti=1e3; p.oc.rhoi=1; p.oc.wT=8; p.oc.mtom=0;
p.tomopt.err=1e-6; p.oc.derr=1; %p.oc.derr=1e-2;
p2=isc(p,0:0.1:0.1); %0.5);
%%
p3=isc(p2,0.25:0.25:0.5);
%%
p4=isc(p3,0.7:0.2:0.9);
%%
polldiagn(p3,5,20,10,1,1);
%% CPS to hCPS
hdir='h2'; hpt='pt17'; fdir='FSS'; fpt='pt13'; 
[muv1, muv2,ind,h]=floqpsap(hdir,hpt); fprintf('d(u_H)=%i\n',ind); 
q=loadp(hdir,hpt); % load hopf orbit
q.hopf.tom.Nmax=1e10; % switch of mesh control of mtom 
huclean(q); hoplot(q,1,1); 
if 1 % refine hopf orbit in time steps for better projection
q=uhopftref(q,1.5); % refine mesh
[y,T,lam,res,iter,A,q]=honloop(q,q.hopf.y,q.hopf.T,q.u(q.nu+q.nc.ilam)); % recompute sol on finer grid
q.hopf.y=y; q.hopf.T=T;
res, hoplot(q,2,1); 
end
%% Cell2: computation the of canonical path from FSS/pt9 to h2/pt17, al=0.95
p=[];
p=ocinit(p,fdir,fpt,hdir,hpt); % load standard options
p.oc.s1=q; % overwrite starting states with finer CPS from Cell 1
%p.oc.s1=leastdis(p.oc.s1,p.oc.s0.u);
p.oc.nti=1e3; p.oc.rhoi=1; p.oc.wT=10; p.oc.mtom=0;
p.tomopt.err=1e-6; p.oc.derr=1; %p.oc.derr=1e-2;
p5=isc(p,0:0.1:0.1); %0.5);
%%
p6=isc(p5,0.2:0.2:0.6); %
%%
cpplot(p6,4); 
%%
polldiagn(p6,11,20,10,1,1);
