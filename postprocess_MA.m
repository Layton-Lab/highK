clear all;

fname = './Sensitivity/2024-01-17_MA_ee_var-Kplas_r-100_notes-new50days.csv'; % 50 days TODO -- update with new ranges
T_Kplas = readtable(fname, 'ReadRowNames', true);
EEplas = table2array(T_Kplas);
fname = './Sensitivity/2024-01-17_MA_ee_var-Kmusc_r-100_notes-new50days.csv'; % 50 days TODO -- update with new ranges
T_Kmusc = readtable(fname, 'ReadRowNames', true); 
EEmusc = table2array(T_Kmusc);

%% set parameter names
parnames = T_Kplas.Properties.VariableNames;
pname_plt = {};
for ii = 1:length(parnames)
    pname_plt{ii} = change_parname(parnames{ii});
end

%% Process elementary effects data
% compute mu, mustar, sigma, MI for each parameter
% row 1: mu, row 2: mustar, row3: sigma, row4: MI
plasvals = zeros(4, length(parnames));
muscvals = zeros(4, length(parnames));

plasvals(1,:) = mean(EEplas);
plasvals(2,:) = mean(abs(EEplas));
plasvals(3,:) = std(EEplas);
plasvals(4,:) = sqrt(plasvals(1,:).^2 + plasvals(3,:).^2);

muscvals(1,:) = mean(EEmusc);
muscvals(2,:) = mean(abs(EEmusc));
muscvals(3,:) = std(EEmusc);
muscvals(4,:) = sqrt(muscvals(1,:).^2 + muscvals(3,:).^2);

%% Make figures
% Morris plots
cmap = turbo(length(pname_plt));
marksize=10; ms = 'diamond';
fx = 35; fy = 35; fleg = 16; ft = 22;
ftxt = 15; fgca = 18;
figure(1)
clf;
nr = 1; nc = 2;
% Kplas
subplot(nr,nc,1)
hold on
for jj = 1:length(pname_plt)
    id = jj;
    plot(plasvals(2,id), plasvals(3,id),...
        'markersize', marksize, 'marker', ms,...
        'color', cmap(jj, :), 'MarkerFaceColor', cmap(jj,:),...
        'linestyle', 'none', 'DisplayName', pname_plt{id})
    if or(plasvals(2,id) > 0.1, plasvals(3,id) > 0.1)
        if id == 1 % kgut
            dx = 0.1; dy = -0.1;
        elseif id == 2 % GFRbase
            dx = 0.025;dy = -0.1;
        elseif id == 3 % eta_ptKreab_base
            dx = -0.1; dy = 0.15;
        elseif id == 4 % eta_LOHKreab
            dx = 0.1; dy = 0.1;
        elseif id == 5 % eta_ptKreab
            dx = 0.1; dy = 0.1;
        elseif id == 9 % alpha_TGF
            dx = -0.15; dy = 0.14;
        elseif id == 15 % alpha_al
            dx = 0.05; dy = -0.05;
        elseif id == 16 % beta_al
            dx = 0.05; dy = -0.05;
        elseif plasvals(3,id) > 1.0
            dx = 0.1; dy = -0.1;
        else
            dx = 0.1; dy = -0.1;
        end
        text(plasvals(2,id) + dx, plasvals(3,id) + dy, pname_plt{id},...
            'fontsize', ftxt)
    end

end
set(gca, 'fontsize', fgca)
xlabel('\mu*', 'fontsize', fx)
ylabel('\sigma', 'fontsize', fy)
xlim([0.0, 4.0])
ylim([0.0, 4.0])
title('Plasma [K^+]', 'fontsize', ft)
grid on

% Kmusc
subplot(nr,nc,2)
hold on
for jj = 1:length(pname_plt)
    id = jj;
    plot(muscvals(2,id), muscvals(3,id),...
        'markersize', marksize, 'marker', ms,...
        'color', cmap(jj,:), 'MarkerFaceColor', cmap(jj,:),...
        'linestyle', 'none', 'DisplayName', pname_plt{id})
    if or(muscvals(2,id) >10, muscvals(3,id) > 10)
        if id == 1 % kgut
            dx = 3; dy = 6;
        elseif id == 2 % GFRbase
            dx = -1; dy = -6;
        elseif id == 9 % alphaTGF
            dx = -15; dy = 7;
        elseif id == 11 % C_al_base
            dx = -5; dy = 8;
        elseif id == 12 % mKALDO
            dx = 1; dy = -4;
        elseif id == 14 % Ainsulin
            dx = -4; dy = -6;
        elseif id == 15 % alpha_al
            dx = 0.0; dy = -5;
        elseif id == 16 % beta_al
            dx = 0.0; dy = -5;
        elseif muscvals(2,id) > 150
            dx = 0; dy = 5;
        else
            dx = 0.0; dy = 5.0;
        end
        text(muscvals(2,id) + dx, muscvals(3,id) + dy, pname_plt{id},...
            'fontsize', ftxt)
    end
end
set(gca, 'fontsize', fgca)
xlabel('\mu*', 'fontsize', fx)
ylabel('\sigma', 'fontsize', fy)
ylim([0.0, 200])
xlim([0.0, 200])
% xticks(0.0:25:150)
% yticks(0.0:25:150)
title('Intracellular [K^+]', 'fontsize', ft)
grid on

legend('fontsize',fleg)

AddLetters2Plots(figure(1), {'(A)', '(B)'},...
                'HShift', -0.05, 'VShift', -0.06, ...
                'fontsize', 20)

%% Morris Index Plots
figure(2)
clf;
nr = 2; nc = 1;
cmap = turbo(18);
marksize=15; ms = 'o';
fx = 22; fy = 22; fleg = 16; ft = 22;
ftxt = 16; fgca = 18;
cid = 16; cid2 = cid;
[~, Mids] = sort(muscvals(4,:), 'descend'); % sort Morris Indices by muscle

% Kplas
subplot(nr, nc, 1)
plot(plasvals(4,Mids), 'linestyle', 'none', 'markersize', marksize,...
                'marker', ms, 'color', cmap(cid,:),...
                'markerfacecolor', cmap(cid,:))
set(gca, 'fontsize', 16)
xticks(1:length(Mids))
xlim([1,length(Mids)])
xlabel('Parameter', 'fontsize', fx)
ylabel('Morris Index', 'fontsize', fy)
title('Plasma [K^+]', 'fontsize', ft)
xticklabels(pname_plt(Mids))
grid on

% Kmusc
subplot(nr,nc,2)
plot(muscvals(4,Mids), 'linestyle', 'none', 'markersize', marksize,...
        'marker', ms, 'color', cmap(cid,:),...
        'markerfacecolor', cmap(cid,:))
set(gca, 'fontsize', 16)
xticks(1:length(Mids))
xlim([1,length(Mids)])
xlabel('Parameter', 'fontsize', fx)
ylabel('Morris Index', 'fontsize', fy)
title('Intracellular [K^+]', 'fontsize', ft)
xticklabels(pname_plt(Mids))
grid on
AddLetters2Plots(figure(2), {'(A)', '(B)'},...
                'HShift', -0.05, 'VShift', -0.06, ...
                'fontsize', 20)

