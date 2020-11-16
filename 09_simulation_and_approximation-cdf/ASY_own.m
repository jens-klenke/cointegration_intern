
                
                % Loop over repetitions
                parfor j=1:rep;
                    u=randn(T,k+1); % Draw random Shocks
                    W1=cumsum(u(:,1:k))/sqrt(T); % Simulate Wiener Process
                    u12=sqrt(R2run./(1-R2run)).*u(:,1:k)*ones(k,1)/sqrt(k) + u(:,k+1);
                    J12=Bc(u12,c_run); %Ohrnsetin Uhlenbeck Process
                    % Corrections according to case
                    
                    if dets==1; % No Constant, no trend
                        W1d=W1;
                        J12dc=J12;
                    elseif dets==2; %Constant, no trend
                        W1d = W1 - (ones(T,1)*mean(W1));
                        J12dc = J12 - (ones(T,1)*mean(J12));
                    elseif dets==3; % Constant and Trend
                        W1d=W1-(4-6*repmat(lambda',1,k)).*(ones(T,1)*mean(W1))-(12*repmat(lambda',1,k)-6).*(ones(T,1)*mean(repmat(lambda',1,k).*W1));
                        J12dc=J12-(4-6*lambda').*(ones(T,1)*mean(J12))-(12*lambda'-6).*(ones(T,1)*mean(lambda'.*J12));
                    end;
                    Wdc=[W1d J12dc];
                    
                    
                    
                    % -----------------------------Common Terms------------
                    WdcDW2=mean(Wdc(1:T-1,:).*(repmat(u(2:T,k+1),1,k+1)));
                    WdcWdci=inv(1/T.^2*Wdc'*Wdc);
                    W1dW1di=inv(1/T*W1d(1:T-1,:)'*W1d(1:T-1,:));
                    W1dJ12dc=mean(W1d(1:T-1,:).*repmat(J12dc(1:T-1,:),1,k));
                    J12dc_sq=mean(J12dc(1:T-1).^2);
                    J12DW2=mean(J12dc(1:T-1).*u(2:T,k+1));
                    % -------------------------------- Boswijk------------------------ 
                    BoswijkStat(j)=c_run.^2 .*J12dc_sq + 2*c_run.*sqrt(T)*J12DW2 + WdcDW2*WdcWdci*WdcDW2';
                    if (cc == 1 && rr == 1);
                        NullStatBoswijk(j)=BoswijkStat(j);
                    else
                        BoswijkPValue(j)=1-minindc(abs(BoswijkStat(j)-NullDistrBoswijk))/rep+10.^(-1000);
                    end;
                    
                    % -------------------------------- Johansen -------------------------------- */
                    
                    Gc=mean(Wdc.*repmat(J12dc,1,size(Wdc,2)))'*([zeros(1,k) c_run])/sqrt(T); % the square root is to make it fit with other extra power, see notes */
                    Wdc_dW_pr=1/T*u(2:T,:)'*Wdc(1:T-1,:);
                    dW_Wdc_pr=1/T*Wdc(1:T-1,:)'*u(2:T,:);
                    
                    JohansenStat(j)=max(eig(Wdc_dW_pr*WdcWdci*dW_Wdc_pr+Gc'*WdcWdci*dW_Wdc_pr+dW_Wdc_pr'*WdcWdci*Gc+Gc'*WdcWdci*Gc));
                    if (cc == 1 && rr == 1);
                        NullStatJohansen(j)=JohansenStat(j);
                    else
                        JohansenPValue(j)=1-minindc(abs(JohansenStat(j)-NullDistrJohansen))/rep+10.^(-1000);
                    end;
                    
                    % -------------------------------- Engle-Granger -------------------------------- */
                    
                    etadc = [(-W1dW1di*mean(W1d(1:T-1,:).*repmat(J12dc(1:T-1,:),1,k))'); 1];
                    Adc = 1/T*Wdc(1:T-1,:)'*Wdc(1:T-1,:);
                    Dmat=[eye(k) (sqrt(R2run./(1-R2run)).*ones(k,1)/sqrt(k));...
                        (sqrt(R2run./(1-R2run)).*ones(1,k)/sqrt(k)) (1+R2run./(1-R2run))];
                    utilde = [u(:,1:k) u12];
                    Wdc_dWtilde=1/sqrt(T)*Wdc(1:T-1,:)'*utilde(2:T,:);

                    EngleGrangerStat(j) = c_run.*sqrt(etadc'*Adc*etadc)./sqrt(etadc'*Dmat*etadc) +...
                        (etadc'*Wdc_dWtilde*etadc)/(sqrt(etadc'*Dmat*etadc)*sqrt(etadc'*Adc*etadc));
                    
                    if (cc == 1 && rr == 1);
                        NullStatEngleGranger(j)=EngleGrangerStat(j);
                    else
                        EngleGrangerPValue(j)=minindc(abs(EngleGrangerStat(j)-NullDistrEngleGranger))/rep;
                    end;
                    
                    % -------------------------------- ECR (Banerjee) -------------------------------- */
                    zaehler= sqrt(T)*(J12DW2'-W1dJ12dc*W1dW1di*mean(W1d(1:T-1,:).*repmat(u(2:T,k+1),1,k))');
                    nenner=(sqrt(J12dc_sq'-W1dJ12dc*W1dW1di*W1dJ12dc'));
                    ErrCorrStat(j) = c_run.*sqrt(J12dc_sq'-W1dJ12dc*W1dW1di*W1dJ12dc') + zaehler/nenner;
                    
                    
                    if (cc == 1); % here we condition on each R2 because the r2=0 c.v. turned out not to be good for R2 very large >0.75*/
                        NullStatErrCorr(j)=real(ErrCorrStat(j)); % a 0 plus imaginary occurred once... */
                    else
                        ErrCorrPValue(j)=minindc(abs(ErrCorrStat(j)-NullDistrErrCorr))/rep;
                    end;
                    
                end; %rep loop end
                
                %% Write Null Distributions and Critical Values for underlying tests

                NullDistrBoswijk=sort(NullStatBoswijk);
                    BoswijkPValue=1-rankindx(NullStatBoswijk,1)/rep+10.^(-1000);
                

                NullDistrJohansen=sort(NullStatJohansen);
                    JohansenPValue=1-rankindx(NullStatJohansen,1)/rep+10.^(-1000);
                
                NullDistrEngleGranger=sort(NullStatEngleGranger);
                    EngleGrangerPValue=rankindx(NullStatEngleGranger,1)/rep+10.^(-1000);
                
                NullDistrErrCorr=sort(NullStatErrCorr);
                ErrCorrPValue=rankindx(NullStatErrCorr,1)/rep+10.^(-1000);
                

                 
                % Define Statistics
                FisherStatEJ=-2*(log(EngleGrangerPValue)+log(JohansenPValue));
                
                 
            end; % R2 
        end; % c 