% High K simulation with dynamic PT
clear all;

% initialize parameter values
fprintf('loading params \n')
pars = set_params();
[params, parnames] = pars2vector(pars, 0);

%----------------
% Simulation settings
%------------------
MealInsulin = 1; % insulin + KCl meals
len_meal = 30; % length of meal in mins
doFF=1; % do GI FF effect
Kamt_high = 4*78/3;  % high K intake
Kamt_meal = Kamt_high;

n_days = 50;

TGF_eff = 1; %1; % 1: PT + GFR, 2: GFR only, 3: PT only


alpha_TGF = pars.alpha_TGF;

% dynamic eta_PT
dyn_etaPT = 1;
etaPT_0 = 0.67; % starting eta_ptKreab
etaPT_final = 0.36; % final eta_ptKreab % Wang 2023 on high K diet

t_eta_days = 1; % how many days to high K adjustment


%---------------------
%---------------------
% simulation settings
opts.do_insulin = MealInsulin;
opts.do_FF = doFF;
opts.do_TGFeff = [TGF_eff, alpha_TGF, etaPT_0];

% dynamic eta_ptKreab settings
t_etaPT_mins = t_eta_days * 60 * 24; % convert from days
t0 = 0;
m = (etaPT_final - etaPT_0) / (t_etaPT_mins - t0); % slope
b = etaPT_0 - m * t0; % intercept
opts.do_dyn_etaPT = [dyn_etaPT, etaPT_final, etaPT_0, t_etaPT_mins];

%% set initial conditions
temp = load('./SS/SS_4vars.mat');
SS = temp.SS;
[IC, ~, ~] = getSS(SS, params, ...
    'do_insulin', opts.do_insulin,...
    'do_FF',opts.do_FF,...
    'do_figs', 0); % start at SS

%% Intake for multiple days
MealTimes = [6, 12, 18]*60; % meal times
Meal_Kamts = Kamt_meal * ones(size(MealTimes)); % amounts of K per meal

Tvals = {}; Yvals = {};
fprintf('day 1 \n')
[T, Y] = one_day_meals(IC, Meal_Kamts, len_meal, MealTimes, params, opts);
Tvals{1} = T; Yvals{1} = Y; % day 1 simulation

for ii = 2:n_days
    if mod(ii, 10) == 0
        fprintf('day %s \n', num2str(ii))
    end
    IC = Yvals{ii-1}(end,:);
    [T, Y] = one_day_meals(IC, Meal_Kamts, len_meal, MealTimes, params, opts);
    Tvals{ii} = T; Yvals{ii} = Y;
end

%-------------
% save results
%-------------
save_res = input('save results? (0/1) ');
if save_res
    notes = input('notes: ');
    fname = strcat('./MultiDaySim/', date, '_driver_lindep',...
        '_insulin-', num2str(MealInsulin),...
        '_Kamt_meal-', num2str(Kamt_meal),...
        '_etaPT_final-', num2str(etaPT_final),...
        '_etaPT_0-', num2str(etaPT_0),...
        't_eta_days-', num2str(t_eta_days),...
        '_ndays-', num2str(n_days),...
        '_notes-', notes,...
        '.mat');
    save(fname)
    fprintf('results saved to: \n %s \n', fname);
end

%---------------
% plot results
%---------------
T = T./60; % change to hours

%% All the days
T_all = []; Y_all = [];
for ii = 1:n_days
    temp = Tvals{ii}./60 + (24)*(ii - 1);
    T_all = [T_all; temp];
    Y_all = [Y_all; Yvals{ii}];
end
T_all = T_all./24;
figure(2)
clf;
nr = 2; nc = 2;
lw = 3; lwgray = 2; lsgray = '--';
cmap = parula(6);
c1 = cmap(1,:); c2 = cmap(2,:); c3 = cmap(3,:); c4 = cmap(4,:);
cgraymap = gray(5);
cgray = cgraymap(3,:);
subplot(nr,nc,1)
plot(T_all, Y_all(:,1), 'linewidth', lw, color = c1)
xlabel('Time (days)')
ylabel('Gut amount')
title('Gut amount')
grid on

subplot(nr,nc,2)
hold on
plot(T_all,Y_all(:,2)/pars.V_plasma, 'linewidth',lw,'color',c2)
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
xlabel('Time (days)')
ylabel('Plasma [K^+]')
title('Plasma [K^+]')
grid on

subplot(nr,nc,3)
hold on
plot(T_all,Y_all(:,3)/pars.V_interstitial,'linewidth',lw,'color',c3)
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
xlabel('Time (days)')
ylabel('Interstitial [K^+]')
title('Interstitial [K^+]')
grid on

subplot(nr,nc,4)
hold on
plot(T_all,Y_all(:,4)/pars.V_muscle,'linewidth',lw,'color',c4)
yline(120,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(140,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
xlabel('Time (days)')
ylabel('Intracellular [K^+]')
title('Intracellular [K^+]')
grid on


%-------------------------
% Functions used
%-------------------------
function [T, Y] = one_day_meals(IC0, Meal_Kamts, len_meal, MealTimes, params, opts)
    if MealTimes(1) > 0
    % start with fasting simulation if not starting at 0
        [tf, yf] = fast_sim(IC0, [0,MealTimes(1)], -60*6, params, opts);
        T = tf; 
        Y = yf;
    else
        error('MealTimes(1) = 0') % should not have meal at 0
    end
    % Do meal simulations per day
    for ii = 1:length(MealTimes)
        [tm, ym] = meal_sim(Y(end,:), MealTimes(ii), len_meal, Meal_Kamts(ii), ...
                                    params, opts);
        if ii < length(MealTimes)
            t_end = MealTimes(ii+1);
        else
            t_end = 24 * 60;
        end
        [tf, yf] = fast_sim(ym(end,:), [tm(end), t_end], tm(1), ...
                                    params, opts);
        T = [T; tm; tf]; Y = [Y; ym; yf];
    end
end

% Meal simulation
function [t, y] = meal_sim(IC, t0, len_meal, Kamt, params, opts)
    tf = t0 + len_meal; % meal length
    tspan = [t0, tf];
    Kintake = Kamt / (tf - t0);
    options = odeset('RelTol', 1.0e-6, 'AbsTol', 1e-9); % ode solver settings
    [t, y] = ode15s(@(t,y) kreg_eqns(t,y,params,...
                                'do_insulin', opts.do_insulin,...
                                'do_FF', opts.do_FF,...
                                'meal_time', t0,... % t0 is start of meal
                                'Kintake', Kintake,...
                                'TGF_eff', opts.do_TGFeff,...
                                'do_dyn_etaPT', opts.do_dyn_etaPT),...
                                tspan, IC, options);
end

% Fasting simulation
function [t, y] = fast_sim(IC, tspan, last_meal, params, opts)
    options = odeset('RelTol', 1.0e-6, 'AbsTol', 1e-9); % ode solver settings
    [t, y] = ode15s(@(t,y) kreg_eqns(t,y,params,...
                                'do_insulin', opts.do_insulin,...
                                'do_FF', opts.do_FF,...
                                'meal_time', last_meal,...
                                'Kintake', 0,... % fasting state
                                'do_dyn_etaPT', opts.do_dyn_etaPT,...
                                'TGF_eff', opts.do_TGFeff),... 
                                tspan, IC, options);
end