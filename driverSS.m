% This script can be used for the "getSS" function

% clear 
clear all;

%%--------------
% User Input
%%-------------

%--------------
%--------------

%% initialize parameter values
fprintf('loading params \n')
pars = set_params();
[params, parnames] = pars2vector(pars,0);

%% set initial conditions
temp = load('./SS/SS_4vars.mat');
IC = temp.SS;
[SS, exitflag, residual] = getSS(IC, params,...
                                    'do_figs', true); % show figs?

% get kidney values from SS 
v = compute_kidney_vars(SS, params, ...
                        'SS', true);
printSS=true;
print_renal = true;
if printSS
   %% volumes
   pars = set_params();
   fprintf('final steady states \n')
   fprintf('               SS \n')
   fprintf('M_Kgut         %0.4f\n', SS(1))
   fprintf('M_Kplas        %0.4f\n', SS(2))
   fprintf('M_Kinter       %0.4f\n', SS(3))
   fprintf('M_Kmuscle      %0.4f\n', SS(4))
   fprintf('\n')
   fprintf('K_plas         %0.4f\n', SS(2)/pars.V_plasma)
   fprintf('K_inter        %0.4f\n', SS(3)/pars.V_interstitial)
   fprintf('K_muscle       %0.4f\n', SS(4)/pars.V_muscle)
end

if print_renal
    fprintf('\n')
    fprintf('Kidney steady states: \n')
    fprintf('Phi_{filK}        %0.4f\n', v.filK)
    fprintf('Phi_{psKreab}     %0.4f\n', v.psKreab)
    fprintf('Phi_{dtKsec}      %0.4f\n', v.dtKsec)
    fprintf('Phi_{cdKsec}      %0.4f\n', v.cdKsec)
    fprintf('Phi_{cdKreab}     %0.4f\n', v.cdKreab)
    fprintf('UrineK            %0.4f\n', v.UrineK)
end

