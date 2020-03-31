T = 10;
k = 1;
u=randn(T,k+1); % Draw random Shocks

W1=cumsum(u(:,1:k))/sqrt(T); % Simulate Wiener Process
                    u12=sqrt(R2run./(1-R2run)).*u(:,1:k)*ones(k,1)/sqrt(k) + u(:,k+1);
                    J12=Bc(u12,c_run); %Ohrnsetin Uhlenbeck Process