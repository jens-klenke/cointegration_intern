function Y=DGP(START,Coeff, nlag, err, ecvectors, r)
% Coefficients must be ordered as: differences variable1 lag1 ...
% variable1 lag_n, variable2 lag_1, ..., variable_J lag_n, levels variable1
% ...variable_J, intercept; THE DIFFERENT VARIABLES ARE IN COLUMNS
T=size(err,1);
nvar=size(err,2);
Y=zeros(T+nlag+1,nvar);
Y(1:nlag+1,:)=START;
dY=Y(2:nlag+1,:)-Y(1:nlag,:);
dx=zeros(T,nvar);
dY=dY(end:-1:1,:); % sort in reverse order, newest on top
Fix=Coeff(end,:);
Coeff(end,:)=[];
if r>0
    for t=1:T
        Z=[reshape(dY,numel(dY),1); reshape(Y(nlag+t,:)*ecvectors(1:r,:)', r,1); 1];
        dx = Z'*Coeff+err(t,:);
        dY=[dx;dY(1:end-1,:)];
        Y(nlag+1+t,:)=Y(nlag+t,:)+dx;
    end
else
    for t=1:T
        %Z=[dY(:); 1];
        dx(t,:)= dY(:)'*Coeff+Fix+err(t,:);
        dY(2:end,:)=dY(1:end-1,:);
        dY(1,:)=dx(t,:);
        %dY=[Temp;dY(1:end-1,:)];
        %dx(t,:) =Temp;
        %Y(nlag+1+t,:)=Y(nlag+t,:)+dx(t,:);
    end
end
dx(1,:)=dx(1,:)+Y(nlag+1,:);
X=cumsum(dx);
Y=[Y(1:nlag+1,:);X];