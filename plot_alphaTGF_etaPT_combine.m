clear all;

date_string = '06-Feb-2024';
notes = 'allsims';

pars = set_params(); % get parameter values

eta_ptKreab_vals = [0.1800,0.2400,0.3000,0.3600,0.4300,0.4900,0.5500,0.6100,0.6700];
alpha_fold = [0.5000, 0.6250, 0.7500, 0.8750, 1.0000, 1.1250, 1.2500, 1.3750, 1.5000];
alpha_TGF_vals = alpha_fold * pars.alpha_TGF;

finalKplas = zeros(length(eta_ptKreab_vals), length(alpha_TGF_vals));
finalKmusc = zeros(length(eta_ptKreab_vals), length(alpha_TGF_vals));

ylabs = cell(size(eta_ptKreab_vals)); % ylabels
xlabs = cell(size(alpha_TGF_vals)); % xlabels

for ii = 1:length(eta_ptKreab_vals)
    for jj = 1:length(alpha_TGF_vals)
        eta_ptKreab = eta_ptKreab_vals(ii);
        alpha_TGF = alpha_TGF_vals(jj);

        [Mplas_final, MIC_final] = get_sim_data(eta_ptKreab, alpha_TGF, ...
                                                date_string, notes);

        finalKplas(ii,jj) = Mplas_final/pars.V_plasma;
        finalKmusc(ii,jj) = MIC_final/pars.V_muscle;

        % x and y labels
        ylabs{ii} = num2str(eta_ptKreab);
        xlabs{jj} = num2str(alpha_fold(jj));


    end % for ii
end % for jj

figure(1)
cmap = turbo;
fsize = 16;
min_fold = min(min(min(finalKplas))/4.2, min(min(finalKmusc))/130);
max_fold = max(max(max(finalKplas))/4.2, max(max(finalKmusc))/130);
clf;
subplot(1,2,1)
h1 = heatmap(finalKplas, ...
                'colormap', cmap);
h1.XData = xlabs;
h1.YData = ylabs;
h1.Title = 'Final plasma [K^+]';
h1.XLabel = '\alpha_{TGF} fold change';
h1.YLabel = '\eta_{pt-Kreab}';
h1.FontSize = fsize;
h1.ColorLimits = [3.75,7.25]; %[min_fold, max_fold] * 4.2;

subplot(1,2,2)
h2 = heatmap(finalKmusc, ...
                'colormap',cmap);
h2.XData = xlabs;
h2.YData = ylabs;
h2.Title = 'Final intracellular [K^+]';
h2.XLabel = '\alpha_{TGF} fold change';
h2.YLabel = '\eta_{pt-Kreab}';
h2.FontSize = fsize;
h2.ColorLimits = [120, 275];

%%-------------------
% Functions
%%-------------------
function [Mplas_final, MIC_final] = get_sim_data(eta_PTKreab, alphaTGF, ...
                                                date_string, notes)
    n_days = 50;
    MealInsulin = 1;
    Kamt_high = 4 * 78 / 3; % high K intake, per meal
    Kamt_meal = Kamt_high;
    TGF_eff = 1;
    fname = strcat('./MultiDaySim/', ...
                    date_string, ...
                    '_driver_multiday',...
                    '_insulin-', num2str(MealInsulin),...
                    '_Kamt_meal-', num2str(Kamt_meal),...
                    '_TGFeff-', num2str(TGF_eff),...
                    '_alphaTGF-', num2str(alphaTGF),...
                    '_etaPTKreab-', num2str(eta_PTKreab),...
                    '_ndays-', num2str(n_days),...
                    '_notes-', notes,...
                    '.mat');
    dat = load(fname);
    Y_lastday = dat.Yvals{50}; % last day simulation

    Mplas_final = Y_lastday(end, 2);
    MIC_final = Y_lastday(end,4);
end