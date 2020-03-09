%% Output
clc
clear
critic=[0.1 0.05 0.01];
TJ=[50 75 100 125 150];
rlabel={'EG-Process, T=50','EG-Process, T=75', 'EG-Process, T=100', 'EG-Process, T=125','EG-Process, T=150','J-Process, T=50', ...
    'J-Process, T=75','J-Process, T=100','J-Process, T=125', 'J-Process, T=150',...
    'Mixed-Process, T=50', 'Mixed-Process, T=75','Mixed-Process, T=100','Mixed-Process, T=125', 'Mixed-Process, T=150'}
clabel={'Bayer Hanck, bootstrap', 'Bayer Hanck, bootstrap 4', 'Bayer Hanck, Hartung',  'Bayer Hanck, Hartung 4','Bayer Hanck, naive', 'Johansen lmax, bootstrap', 'Engle-Granger, bootstrap', 'Phillips', 'Johansen lmax, asy', 'Engle-Granger, asy'}
names={' DGP(B)' ' DGP(A) ' ' DGP(C)'}
close all
for z=1:3
    criti=critic(z);
    A=zeros(3*length(TJ),10);
    C=A;
    for k=1:length(TJ)
        t=TJ(k);
        savename=['results_T_all_tests_varlag_2=' num2str(t)];
        load(savename)
        A(k,:)=mean([p_hut_size<criti; p1_fac_size<criti; p_hut_hart_size<criti; p_hart_hut_fix_size<criti; ...
            min(p1_size,p2_size)<criti; p1_size<criti; p2_size<criti; p14_size<criti; ps1size_plmax<criti; ps1size_pcadf<criti]');
        A(length(TJ)+k,:)=mean([p2_hut_size<criti; p2_fac_size<criti; p2_hut_hart_size<criti; p2_hart_hut_fix_size<criti; ...
            min(p21_size,p22_size)<criti; p21_size<criti; p22_size<criti; p24_size<criti; ps2size_plmax<criti; ps2size_pcadf<criti]');
        A(2*length(TJ)+k,:)=mean([p3_hut_size<criti; p3_fac_size<criti; p3_hut_hart_size<criti; p3_hart_hut_fix_size<criti; ...
            min(p31_size,p32_size)<criti; p31_size<criti; p32_size<criti; p32_size<criti; ps3size_plmax<criti; ps3size_pcadf<criti]');
        C(k,:)=mean([p_hut_power<criti; p1_fac_power<criti; p_hut_hart_power<criti; p_hart_hut_fix_power<criti; ...
            min(p1_power,p2_power)<criti; p1_power<criti; p2_power<criti; p14_power<criti; ps1power_plmax<criti; ps1power_pcadf<criti]');
        C(length(TJ)+k,:)=mean([p2_hut_power<criti; p2_fac_power<criti;  p2_hut_hart_power<criti; p2_hart_hut_fix_power<criti; ...
            min(p21_power,p22_power)<criti; p21_power<criti; p22_power<criti; p24_power<criti; ps2power_plmax<criti; ps2power_pcadf<criti]');
        C(2*length(TJ)+k,:)=mean([p3_hut_power<criti; p3_fac_power<criti;  p3_hut_hart_power<criti; p3_hart_hut_fix_power<criti; ...
            min(p31_power,p32_power)<criti; p31_power<criti; p32_power<criti; p32_power<criti; ps3power_plmax<criti; ps3power_pcadf<criti]');
    end
    tablename1=['Size_for_alpha=' num2str(criti*100) '%.tex']
    tablename2=['Power_for_alpha=' num2str(criti*100) '%.tex']
    figurename=['Power for \alpha=' num2str(criti*100) '%']
    figurename2=['Size for \alpha=' num2str(criti*100) '%']
    matrix2latex(A, tablename1, 'rowLabels', rlabel,'columnLabels', clabel);
    matrix2latex(C, tablename2, 'rowLabels', rlabel,'columnLabels', clabel);
    figure
    for run=2:-1:1
        subplot(1,2,3-run)
        plot(TJ,C(5*(run-1)+1:5*run,3),'r','LineWidth',2)
        hold on
        plot(TJ,C(5*(run-1)+1:5*run,1),'r--','LineWidth',2)
        hold on
        plot(TJ,C(5*(run-1)+1:5*run,6),':d','LineWidth',2)
        plot(TJ,C(5*(run-1)+1:5*run,7),':s','LineWidth',2)
        legend( '\tau^*','\chi^*','\lambda_{max}^*', 'AEG^*')
        xlabel('time series length')
        title([figurename names(run)])
    end
    figure
    for run=2:-1:1
        subplot(1,2,3-run)
        plot(TJ,A(5*(run-1)+1:5*run,3),'r','LineWidth',2)
        hold on
        plot(TJ,A(5*(run-1)+1:5*run,1),'r--','LineWidth',2)
        hold on
        plot(TJ,A(5*(run-1)+1:5*run,6),':d','LineWidth',2)
        plot(TJ,A(5*(run-1)+1:5*run,7),':s','LineWidth',2)
        ylim([0 0.2])
        legend( '\tau^*','\chi^*','\lambda_{max}^*', 'AEG^*')
        xlabel('time series length')
        title([figurename2 names(run)])
    end
end

