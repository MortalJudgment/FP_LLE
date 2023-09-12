cp=p6.cp; ts=4; nx=20; np=21; x=linspace(0,1,nx+1)'; T=real(cp.par); tv=T*cp.t(1:ts:end); 
u=cp.u; tl=length(tv); ov=ones(1,nx+1); ga=q.u(q.nu+5); 
for k=1:5 % x-t-plots of solns   
    figure(k); clf; hold on;    
    if k<5; z=u((k-1)*np+1:k*np,1:ts:end)'; % states and co-states 
    else % control 
        l1=u(2*np+1:3*np,1:ts:end)';
        z=-(1+l1)./ga; % control 
    end
    surf(x,tv,z);
    ylabel('t'); view(70,30); 
    switch k; case 1; title('v');  case 2; title('w'); 
    case 3; title('\lambda');  case 4; title('\mu'); case 5; title('q'); 
    end 
    set(gca,'fontsize',14); 
end
%% plot of J_c over R,B 
figure(5);clf;
jc=[];
for k=1:size(p2.cp.u,2)
    p2.oc.s1.np=2;
    jc=[rocjcf(p2.oc.s1,[p2.cp.u(:,k);p2.oc.s1.u(p2.oc.s1.nu+1:end)]),jc];
end
for k=1:10:size(p2.cp.u,2)
    plot3(p2.cp.u(1,end-k:end),p2.cp.u(3,end-k:end),sum(jc(:,end-k:end),1)*1/2); 
    %pause(0.1);
end
xlabel('R_1'); ylabel('B_1'); zlabel('J_c');