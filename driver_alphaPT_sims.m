% High K simulation experiment
% Options: PT impact of high K intake
clear all;

% initialize parameter values
fprintf('loading params \n')
pars = set_params();
[params, parnames] = pars2vector(pars, 0);

save_res = input('save results? (0/1) ');
notes = input('notes: ');

eta_ptKreab_vals = [0.1800,0.2400,0.3000,0.3600,0.4300,0.4900,0.5500,0.6100,0.6700];

alpha_TGF_vals = [0.5000, 0.6250, 0.7500, 0.8750, 1.0000, 1.1250, 1.2500, 1.3750, 1.5000] * pars.alpha_TGF;


for ii_etaPTKreab = 1:length(eta_ptKreab_vals)
    for jj_alphaTGF = 1:length(alpha_TGF_vals)

        % Vary over alpha_TGF and eta PT Kreab
        eta_ptKreab = eta_ptKreab_vals(ii_etaPTKreab);% get from vals
        alpha_TGF = alpha_TGF_vals(jj_alphaTGF); % get from vals

        fprintf('eta_ptKreab: %0.6f, alpha_TGF: %0.6f \n', eta_ptKreab, alpha_TGF)


        %------------- 
        % Settings
        %-------------
        MealInsulin = 1; % insulin + KCl meals
        len_meal = 30; % length of meal in minutes
        doFF = 1; % do FF effect on DT
        Kamt_high = 4 * 78 / 3; % high K intake, per meal
        Kamt_control = 78 / 3; % control K intake, per meal (seems to keep stable here...)
        Kamt_meal = Kamt_high; %Kamt_high; %Kamt_high; %Kamt_control; % control sim
        
        doMKX = 0; % do MKX in the simulation 1: DT K sec , 2: CDKsec, 3: CDKreab
        
        n_days = 50 % number of days for simullation
        
        TGF_eff = 1; %0; %3; %1 %1; % do TGF_eff (1: PT + GFR, 2: GFR only, 3: PT only)
        % TGF_eff parameters
        
        %------------------
        %------------------
        % sim settings
        opts.do_insulin = MealInsulin;
        opts.do_FF = doFF; 
        
        
        
        
        if doMKX > 0
            if doMKX == 1
                % slope tries (0.005, 0.01, 0.025, 0.05, 0.075, 0.1)
                MKXslope = 0.1; % dtKsec slope
            elseif doMKX == 2
                MKXslope = 0.1; % cdKsec slope
            elseif doMKX == 3
                MKXslope = -0.1; % cdKreab slope
            end
        else
            MKXslope = -1;
        end
        opts.do_MKX = [doMKX, MKXslope];
        opts.do_TGFeff = [TGF_eff, alpha_TGF, eta_ptKreab];
        
        %% set initial conditions
        temp = load('./SS/SS_4vars.mat');
        SS = temp.SS;
        [IC, ~, ~] = getSS(SS, params, ...
            'do_insulin', opts.do_insulin,...
            'do_FF',opts.do_FF,...
            'do_MKX', opts.do_MKX,...
            'do_figs', 0); % start at SS
        
        % No TGF effect because it is only under high intake (i.e., would need to 
        %    adjust Phi_Kin_ss as well....
        
        %% Intake for multiple day s
        
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
        if save_res
            if doMKX > 0
                fname = strcat('./MultiDaySim/', date, '_driver_multiday',...
                            '_insulin-', num2str(MealInsulin),...
                            '_Kamt_meal-', num2str(Kamt_meal),...
                            '_MKX-', num2str(doMKX),...
                            '_MKXSlope-', num2str(MKXslope),...
                            '_TGFeff-', num2str(TGF_eff),...
                            '_alphaTGF-', num2str(alpha_TGF),...
                            '_etaPTKreab-', num2str(eta_ptKreab),...
                            '_ndays-', num2str(n_days),...
                            '_notes-', notes,...
                            '.mat');
            else
                fname = strcat('./MultiDaySim/', date, '_driver_multiday',...
                    '_insulin-', num2str(MealInsulin),...
                    '_Kamt_meal-', num2str(Kamt_meal),...
                    '_TGFeff-', num2str(TGF_eff),...
                    '_alphaTGF-', num2str(alpha_TGF),...
                    '_etaPTKreab-', num2str(eta_ptKreab),...
                    '_ndays-', num2str(n_days),...
                    '_notes-', notes,...
                    '.mat');
            end
            save(fname)
            fprintf('results saved to: \n %s \n', fname);
        end
    end % for ii etaPTkrab
end % for jj alphaTGF


%%
%------------
% functions
%------------
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
                                'do_MKX', opts.do_MKX),...
                                tspan, IC, options);
%     vals = compute_vars(t,y,params,...
%                             'do_insulin', opts.do_insulin,...
%                             'do_FF', opts.do_FF, ...
%                             'meal_time', t0,...
%                             'TGF_eff', opts.do_TGFeff,...
%                             'do_MKX', opts.do_MKX);
end

% Fasting simulation
function [t, y] = fast_sim(IC, tspan, last_meal, params, opts)
    options = odeset('RelTol', 1.0e-6, 'AbsTol', 1e-9); % ode solver settings
    [t, y] = ode15s(@(t,y) kreg_eqns(t,y,params,...
                                'do_insulin', opts.do_insulin,...
                                'do_FF', opts.do_FF,...
                                'meal_time', last_meal,...
                                'Kintake', 0,... % fasting state
                                'do_MKX', opts.do_MKX,...
                                'TGF_eff', opts.do_TGFeff),... 
                                tspan, IC, options);
%     vals = compute_vars(t,y,params,...
%                             'do_insulin', opts.do_insulin,...
%                             'do_FF', opts.do_FF,...
%                             'meal_time', last_meal,...
%                             'Kintake', 0,...
%                             'do_MKX', opts.do_MKX,...
%                             'TGF_eff', opts.do_TGFeff); 
end