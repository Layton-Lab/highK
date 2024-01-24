% Make a table of all the model predictions

clear all;

%% load data
% simulation data files
f1 = './MultiDaySim/15-Jan-2024_driver_multiday_insulin-1_Kamt_meal-104_TGFeff-0_alphaTGF-0.11694_etaPTKreab-0.36_ndays-50_notes-noPTTGF.mat';  
f2 = './MultiDaySim/15-Jan-2024_driver_multiday_insulin-1_Kamt_meal-104_TGFeff-1_alphaTGF-0.11694_etaPTKreab-0.36_ndays-50_notes-PTTGF.mat'; 
f3 = './MultiDaySim/15-Jan-2024_driver_multiday_insulin-1_Kamt_meal-104_TGFeff-3_alphaTGF-0.11694_etaPTKreab-0.36_ndays-50_notes-PTonly.mat'; 
f4 = './MultiDaySim/15-Jan-2024_driver_multiday_insulin-1_Kamt_meal-26_TGFeff-0_alphaTGF-0.11694_etaPTKreab-0.36_ndays-50_notes-control.mat'; 

% load the files
dat1 = load(f1);
dat2 = load(f2);
dat3 = load(f3);
dat4 = load(f4);

% labels
lab1 = 'High K^+ - no PT/TGF effect';
lab2 = 'High K^+ - PT + TGF effects';
lab3 = 'High K^+ - only PT effect'; 
lab4 = 'Control K^+';

% Day 1
day = 1;
dat1_day1 = get_vals(dat1, day);
dat2_day1 = get_vals(dat2, day);
dat3_day1 = get_vals(dat3, day);
dat4_day1 = get_vals(dat4, day);

% Table for Day 1
T1 = struct2table([dat1_day1; dat2_day1; dat3_day1; dat4_day1],...
                    'RowNames', {lab1, lab2, lab3, lab4});

writetable(T1,'endofday1.csv',...
                    'WriteRowNames',true);

% Day 50
day = 50;
dat1_day50 = get_vals(dat1, day);
dat2_day50 = get_vals(dat2, day);
dat3_day50 = get_vals(dat3, day);
dat4_day50 = get_vals(dat4, day);

% Table for Day 50
T2 = struct2table([dat1_day50; dat2_day50; dat3_day50; dat4_day50],...
                    'RowNames', {lab1, lab2, lab3, lab4});

writetable(T1,'endofday50.csv',...
                    'WriteRowNames', true);






%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%
% Get values
function v_end = get_vals(dat, day)
    % Input: data, day
    % Output: values at the end of the day
    Y_day = dat.Yvals{day};
    Y_day_end = Y_day(end,:);

    v = compute_kidney_vars(Y_day, dat.params,...
                                'do_MKX', [dat.doMKX, dat.MKXslope], ...
                                'do_FF', dat.doFF, ...
                                'TGF_eff', [dat.TGF_eff, ...
                                        dat.alpha_TGF, ...
                                        dat.eta_ptKreab]);
    v_end.MKgut = Y_day_end(1);
    v_end.MKplas = Y_day_end(2);
    v_end.MKinter = Y_day_end(3);
    v_end.MKmuscle = Y_day_end(4);

    fields = fieldnames(v);
    for i = 1:numel(fields)
        fieldName = fields{i};
        temp = v.(fieldName);
        val_end = temp(end);
        v_end.(fieldName) = val_end;
    end
end


