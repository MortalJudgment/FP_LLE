function userplot(p,wnr) 
%
% Specific plot function : in the case of a synchronized solution, the
% solution is symetrized with respect to the y axis to depicted a solution
% over [-pi,pi]
% 

par = p.u(p.nu+1:end) ; % parameters par=[F,alpha,beta,epsilon,lx]
lx=par(4);
u1=p.u(1:p.np);
u2=p.u(p.np+1:p.nu);
if round(lx)==3 % the solution is symetrized
    nbpt=2*p.np-1;
    U1=[fliplr(transpose(u1(2:end))),transpose(u1)];
    U2=[fliplr(transpose(u2(2:end))),transpose(u2)];
    x=linspace(-pi,pi,nbpt);
else
    U1=transpose(u1);
    U2=transpose(u2);
    x=linspace(-pi,pi,p.np);
end
%
% plot real part
%
figure(wnr);
clf;
plot(x,U1,'-g');
% title([p.file.dir + ': pt '+ mat2str(p.file.count-1), 'Re']);
title('Re')
bb=axis; 
bby(1)=bb(3);bby(2)=bb(4);
if bb(4)-bb(3) < 1
    mbb=0.5*(bb(4)+bb(3));
    bby(1)=mbb-0.5;
    bby(2)=mbb+0.5;
end
xlim([-pi,pi]);
ylim([bby(1),bby(2)]);
ytickformat('%.1f')
set(gca,'FontSize',11);
xticks([-pi 0 pi])
xticklabels({'-\pi','0','\pi'})
%
% plot imaginary part
%
figure(wnr+1);
clf;
plot(x,U2,'-r');
% title([p.file.dir + ': pt '+ mat2str(p.file.count-1), 'Im']);
title('Im')
bb=axis;
bby(1)=bb(3);bby(2)=bb(4);
if bb(4)-bb(3) < 1
    mbb=0.5*(bb(4)+bb(3));
    bby(1)=mbb-0.5;
    bby(2)=mbb+0.5;
end
xlim([-pi,pi]);
ylim([bby(1),bby(2)]);
ytickformat('%.1f')
set(gca,'FontSize',11);
xticks([-pi 0 pi])
xticklabels({'-\pi','0','\pi'})
