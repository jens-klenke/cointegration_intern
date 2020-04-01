
 
function B = Bc(u,d);
T=rows(u);
rho=(1+d/T);
v=zeros(size(u));
v(:,1)=u(:,1);
for t=2:T
    v(t,:)=rho*v(t-1,:)+u(t,:);
end

B=v/sqrt(T)
end