function u=gen_ar1(T,rho,eps)
% phi=toeplitz([1 cumprod(rho*ones(1,T-1))],[1 zeros(1,T-1)]);
% u=phi*eps;

% Alternative way

u=zeros(1,T);
u(1)=eps(1);
for t=2:T
    u(t)=rho*u(t-1)+eps(t);
end
u=u';