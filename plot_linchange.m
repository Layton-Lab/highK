% Plot simulation results from driver_lindepPT.m
% where the eta_ptKreab is changed linearly starting from day 0 to day n

clear all;

%% load data
f1 = './MultiDaySim/15-Jan-2024_driver_multiday_insulin-1_Kamt_meal-104_TGFeff-1_alphaTGF-0.11694_etaPTKreab-0.36_ndays-50_notes-PTTGF.mat';
f5 = './MultiDaySim/17-Jan-2024_driver_lindep_insulin-1_Kamt_meal-104_etaPT_final-0.36_etaPT_0-0.67t_eta_days-30_ndays-50_notes-30days.mat';
f4 = './MultiDaySim/17-Jan-2024_driver_lindep_insulin-1_Kamt_meal-104_etaPT_final-0.36_etaPT_0-0.67t_eta_days-20_ndays-50_notes-20days.mat';
f3 = './MultiDaySim/17-Jan-2024_driver_lindep_insulin-1_Kamt_meal-104_etaPT_final-0.36_etaPT_0-0.67t_eta_days-10_ndays-50_notes-10days.mat'; %'./MultiDaySim/22-Oct-2023_driver_multiday_insulin-1_Kamt_meal-104_TGFeff-3_alphaTGF-0.11694_etaPTKreab-0.36_ndays-50_notes-PTonly.mat';
f2 = './MultiDaySim/17-Jan-2024_driver_lindep_insulin-1_Kamt_meal-104_etaPT_final-0.36_etaPT_0-0.67t_eta_days-5_ndays-50_notes-5days.mat'; %'./MultiDaySim/22-Oct-2023_driver_multiday_insulin-1_Kamt_meal-26_TGFeff-3_alphaTGF-0.11694_etaPTKreab-0.67_ndays-50_notes-control.mat';

dat1 = load(f1);
dat2 = load(f2);
dat3 = load(f3);
dat4 = load(f4);
dat5 = load(f5);


lab1 = 'High K^+ - PT + TGF effects';
lab2 = 'High K^+ - linear change (5 days)';
lab3 = 'High K^+ - linear change (10 days)'; 
lab4 = 'High K^+ - linear change (20 days)';
lab5 = 'High K^+ - linear change (30 days)';

%% All the days
T_all1 = []; Y_all1 = [];
T_all2 = []; Y_all2 = [];
T_all3 = []; Y_all3 = [];
T_all4 = []; Y_all4 = [];
T_all5 = []; Y_all5 = [];
for ii = 1:dat1.n_days
    temp = dat1.Tvals{ii}./60 + (24)*(ii - 1);
    T_all1 = [T_all1; temp];
    Y_all1 = [Y_all1; dat1.Yvals{ii}];
    temp = dat2.Tvals{ii}./60 + 24*(ii - 1);
    T_all2 = [T_all2; temp];
    Y_all2 = [Y_all2; dat2.Yvals{ii}];
    temp = dat3.Tvals{ii}./60 + 24 * (ii - 1);
    T_all3 = [T_all3; temp];
    Y_all3 = [Y_all3; dat3.Yvals{ii}];
    temp = dat4.Tvals{ii}./60 + 24 * (ii - 1);
    T_all4 = [T_all4; temp];
    Y_all4 = [Y_all4; dat4.Yvals{ii}];
    temp = dat5.Tvals{ii}./60 + 24 * (ii - 1);
    T_all5 = [T_all5; temp];
    Y_all5 = [Y_all5; dat5.Yvals{ii}];
end
% convert to days
T_all1 = T_all1./24; T_all2 = T_all2./24; 
T_all3 = T_all3./24; T_all4 = T_all4./24;
T_all5 = T_all5./24;

%% Make figures
figure(2)
clf;
nr = 1; nc = 2;
f.labs = 18; f.xlab = 18; f.ylab = 18; f.gca = 18; f.leg = 16; f.title = 22;
lw = 3; lwgray = 4.5; lsgray = ':';
ls1 = '-'; ls2 = '-'; ls3 = '-'; ls4 = '-'; ls5 = '-';
cmap = parula(6);
cmap2 = spring(3);
c1 = cmap(1,:); c2 = cmap(2,:);
c3 = cmap(3,:);
c4 = cmap(4,:);c5 = cmap(5,:);
cgraymap = gray(6);
cgray = cgraymap(1,:);
subplot(nr,nc,1)
hold on
plot(T_all1,Y_all1(:,2)/dat1.pars.V_plasma, 'linewidth',lw,'linestyle', ls1, 'color',c1)
plot(T_all2,Y_all2(:,2)/dat2.pars.V_plasma, 'linewidth',lw,'linestyle', ls2, 'color',c2)
plot(T_all3,Y_all3(:,2)/dat3.pars.V_plasma, 'linewidth',lw,'linestyle', ls3, 'color',c3)
plot(T_all4,Y_all4(:,2)/dat4.pars.V_plasma, 'linewidth',lw,'linestyle', ls4, 'color',c4)
plot(T_all5,Y_all5(:,2)/dat5.pars.V_plasma, 'linewidth',lw,'linestyle', ls5, 'color',c5)
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
set(gca, 'fontsize', f.gca)
xlabel('Time (days)', 'fontsize', f.xlab)
ylabel('Plasma [K^+] (mmol/L)', 'fontsize', f.ylab)
ylim([3.4,6.5])
%title('Plasma [K^+]', 'fontsize', f.title)
grid on
%legend({lab1, lab2,lab3,lab4}, 'fontsize', f.leg, 'location', 'northwest')
legend({lab1, lab2,lab3,lab4,lab5}, 'fontsize', f.leg, 'location', 'northwest')

subplot(nr,nc,2)
hold on
plot(T_all1,Y_all1(:,4)/dat1.pars.V_muscle,'linewidth',lw,'linestyle', ls1, 'color',c1)
plot(T_all2,Y_all2(:,4)/dat2.pars.V_muscle,'linewidth',lw,'linestyle', ls2, 'color',c2)
plot(T_all3,Y_all3(:,4)/dat3.pars.V_muscle,'linewidth',lw,'linestyle', ls3, 'color',c3)
plot(T_all4,Y_all4(:,4)/dat4.pars.V_muscle,'linewidth',lw,'linestyle', ls4, 'color',c4)
plot(T_all5,Y_all5(:,4)/dat5.pars.V_muscle,'linewidth',lw,'linestyle', ls5, 'color',c5)
yline(120,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(140,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
set(gca, 'fontsize', f.gca)
ylim([115,180])
xlabel('Time (days)', 'fontsize', f.xlab)
ylabel('Intracellular [K^+] (mmol/L)', 'fontsize', f.ylab)
%title('Intracellular [K^+]', 'fontsize', f.title)
grid on

legend({lab1, lab2,lab3,lab4,lab5}, 'fontsize', f.leg, 'location', 'northwest')
%legend({lab1, lab2,lab3,lab4}, 'fontsize', f.leg, 'location', 'northwest')

AddLetters2Plots(figure(2), {'(A)', '(B)'},...
                'HShift', -0.05, 'VShift', -0.06, ...
                'fontsize', f.labs)
%% plasma K only (zoomed)
figure(10)
clf;
hold on
plot(T_all1,Y_all1(:,2)/dat1.pars.V_plasma, 'linewidth',lw,'linestyle', ls1, 'color',c1)
plot(T_all2,Y_all2(:,2)/dat2.pars.V_plasma, 'linewidth',lw,'linestyle', ls2, 'color',c2)
plot(T_all3,Y_all3(:,2)/dat3.pars.V_plasma, 'linewidth',lw,'linestyle', ls3, 'color',c3)
plot(T_all4,Y_all4(:,2)/dat4.pars.V_plasma, 'linewidth',lw,'linestyle', ls4, 'color',c4)
plot(T_all5,Y_all5(:,2)/dat5.pars.V_plasma, 'linewidth',lw,'linestyle', ls5, 'color',c5)
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
set(gca, 'fontsize', f.gca)
xlabel('Time (days)', 'fontsize', f.xlab)
ylabel('Plasma [K^+] (mmol/L)', 'fontsize', f.ylab)
title('Plasma [K^+]', 'fontsize', f.title)
xlim([40,50])
grid on

legend({lab1, lab2,lab3,lab4,lab5}, 'fontsize', f.leg, 'location', 'northwestoutside')
%legend({lab1, lab2,lab3,lab4}, 'fontsize', f.leg, 'location', 'northwestoutside')