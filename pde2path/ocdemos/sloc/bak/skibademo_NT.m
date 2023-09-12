% Skiba example, continues cpdemo1D.m
% Either run Cells 3 and 4 from cpdemo1D or, if storing was enabled there,
% load the continuation from p1/pt68 to f1/pt13
% load('skibap.mat')
vv=[]; % vector for objective values of paths to f2/pt12 for alphas in alv
alv=[];
sol={};
doplot=1;
for k=1:35 % loop all computed al values from old path
    p2=[];
    p2=ocinit(p2,'p1','pt68','f2','pt12');
    p2.opt.s0.u(1:42)=p.hist.u{k}(:,1); % reset starting states
    p2.opt.rhoi=1; p2.opt.nti=10; p2.opt.T=100; % set mandantory options
    alvin=[0.1 0.25 0.5 0.75 1];
    p2=NTisc(p2,alvin); % continuation call
    alv=[alv p.hist.alpha(k)]; vv=[vv p2.hist.vv(end)]; % store new al, objectiv value
    Jd=p.hist.vv(k)-vv(end);
    sol{k}=p2.sol.cp; % store solutions
    if abs(Jd)<0.05; doplot=asknu('plot path?',doplot); % Skiba point(s) found
        if doplot==1 % plot the paths to FSC and FSM
            sol1.u=p.hist.u{k}; sol1.t=p.hist.tt{k}; sol1.par=p.hist.par{k};
            psol3Dm_NT(p.ocopt.s1,sol{k},sol1,1,1,[]); view(v); zlabel('P'); % plot both paths
            psol3Dm_NT(p.ocopt.s1,sol{k},sol1,2,0,[]); view(v); zlabel('k'); pause
        end
    end
end
%% plot of objective values over alpha
figure(6); clf; plot(p.hist.alpha(1,:),p.hist.vv(1,:),'-*'); set(gca,'FontSize',14); 
xlabel('\alpha','FontSize',14); ylabel('J_{a}','FontSize',14); hold on;
plot(alv(1,:),vv(1,:),'-*'); set(gca,'FontSize',14); 
xlabel('\alpha','FontSize',14); ylabel('J_{a}','FontSize',14);