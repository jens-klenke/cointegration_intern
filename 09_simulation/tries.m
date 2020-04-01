T = 10;
k = 3;
lambda=(1/T:1/T:1)

u=randn(T,k+1); % Draw random Shocks

W1=cumsum(u(:,1:k))/sqrt(T); % Simulate Wiener Process
u12=sqrt(R2run./(1-R2run)).*u(:,1:k)*ones(k,1)/sqrt(k) + u(:,k+1);
J12=Bc(u12,c_run); %Ohrnsetin Uhlenbeck Process

% det 2
W1d = W1 - (ones(T,1)*mean(W1))
J12dc = J12 - (ones(T,1)*mean(J12))

% det 3
W1d=W1-(4-6*repmat(lambda',1,k)).*(ones(T,1)*mean(W1))-(12*repmat(lambda',1,k)-6).*(ones(T,1)*mean(repmat(lambda',1,k).*W1));
J12dc=J12-(4-6*lambda').*(ones(T,1)*mean(J12))-(12*lambda'-6).*(ones(T,1)*mean(lambda'.*J12));

Wdc=[W1d J12dc];
%Common terms
WdcDW2=mean(Wdc(1:T-1,:).*(repmat(u(2:T,k+1),1,k+1)));
WdcWdci=inv(1/T.^2*Wdc'*Wdc);
W1dW1di=inv(1/T*W1d(1:T-1,:)'*W1d(1:T-1,:));
W1dJ12dc=mean(W1d(1:T-1,:).*repmat(J12dc(1:T-1,:),1,k));
J12dc_sq=mean(J12dc(1:T-1).^2);
J12DW2=mean(J12dc(1:T-1).*u(2:T,k+1));

%Boswijk
BoswijkStat(j)=c_run.^2 .*J12dc_sq + 2*c_run.*sqrt(T)*J12DW2 + WdcDW2*WdcWdci*WdcDW2'

%Johansen
Gc=mean(Wdc.*repmat(J12dc,1,size(Wdc,2)))'*([zeros(1,k) c_run])/sqrt(T); % the square root is to make it fit with other extra power, see notes */
JohansenStat(j)=max(eig(Wdc_dW_pr*WdcWdci*dW_Wdc_pr+Gc'*WdcWdci*dW_Wdc_pr+dW_Wdc_pr'*WdcWdci*Gc+Gc'*WdcWdci*Gc));

%Engel

 EngleGrangerStat(j) = c_run.*sqrt(etadc'*Adc*etadc)./sqrt(etadc'*Dmat*etadc) +...
                        (etadc'*Wdc_dWtilde*etadc)/(sqrt(etadc'*Dmat*etadc)*sqrt(etadc'*Adc*etadc))
                   
                   










