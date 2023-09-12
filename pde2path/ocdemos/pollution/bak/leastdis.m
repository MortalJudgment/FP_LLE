function p=leastdis(p,u0)
y=p.hopf.y;
dis=zeros(size(y,2),1);
for k=1:size(dis,1)
    dis(k)=norm(y(1:size(y,1)/2,k)-u0(1:size(y,1)/2)); % measure distant according to state only
end
[sm,idx]=min(dis);
fprintf('minimal distance: %g \n',sm(1));
yshift=circshift(y(:,1:end-1),size(y,2)-idx(1),2);
fprintf('Numer of mesh-point-shifts: %i \n',size(y,2)-idx(1));
p.hopf.y=[yshift yshift(:,1)];
