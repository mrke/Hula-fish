function [par, metapar, txt_par, chem] = pars_init_Trachinops_caudimaculatus(metadata)

% parameters: initial values at 
T_C = 273.15;     % K, temperature at 0 degrees C
T_ref = T_C + 20; % K, reference temperature
metapar.T_ref = T_ref;
metapar.model = 'abj'; % see online manual for explanation and alternatives

par.z = 1.95;          free.z     = 1;    units.z = '-';          label.z = 'z';          % zoom factor; for z = 1: L_m = 1 cm
par.F_m = 6.5;      free.F_m   = 0;    units.F_m = 'l/d.cm^2'; label.F_m = '{F_M}';    % max spec searching rate
par.kap_X = 0.8;    free.kap_X = 0;    units.kap_X = '-';      label.kap_X = 'kap_X';  % digestion efficiency of food to reserve
par.kap_P = 0.1;    free.kap_P = 0;    units.kap_P = '-';      label.kap_P = 'kap_P';  % faecation efficiency of food to faeces
par.v = 0.04;       free.v     = 1;    units.v = 'cm/d';       label.v = 'v';          % energy conductance
par.kap = 0.97;      free.kap   = 1;    units.kap = '-';        label.kap = 'kap';      % allocation fraction to soma
par.kap_R = 0.95;   free.kap_R = 0;    units.kap_R = '-';      label.kap_R = 'kap_R';  % reproduction efficiency
par.p_M = 50;       free.p_M   = 1;    units.p_M = 'J/d.cm^3'; label.p_M = '[p_M]';    % vol-spec somatic maint
par.p_T =  0;       free.p_T   = 0;    units.p_T = 'J/d.cm^2'; label.p_T = '{p_T}';    % surf-spec somatic maint
par.k_J = 0.002;    free.k_J   = 1;    units.k_J = '1/d';      label.k_J = 'k_J';      % maturity maint rate coefficient
par.E_G = 5400;     free.E_G   = 1;    units.E_G = 'J/cm^3';   label.E_G = '[E_G]';    % spec cost for structure
par.E_Hh = 3.617e-4;free.E_Hh  = 1; units.E_Hh = 'J';       label.E_Hh = 'E_Hh';   % maturity at hatch
par.E_Hb = 2.0e0;free.E_Hb  = 1; units.E_Hb = 'J';       label.E_Hb = 'E_Hb';   % maturity at birth
par.E_Hj = 2.5e0; free.E_Hj  = 1; units.E_Hj = 'J';       label.E_Hj = 'E_Hj';   % maturity at metam
par.E_Hp = 1.5e2; free.E_Hp  = 1; units.E_Hp = 'J';       label.E_Hp = 'E_Hp';   % maturity at puberty
par.h_a = 1.5e-7; free.h_a   = 1; units.h_a = '1/d^2';    label.h_a = 'h_a';     % Weibull aging acceleration
par.s_G = 1e-4;     free.s_G   = 0; units.s_G = '-';        label.s_G = 's_G';     % Gompertz stress coefficient

%% auxiliary parameters
par.T_A = 9800;     free.T_A   = 0; units.T_A = 'K';        label.T_A = 'T_A';     % Arrhenius temp
par.del_M = 0.1209; free.del_M = 1; units.del_M = '-'; label.del_M = 'del_M_SL'; % shape coefficient for standard length
par.del_M_TL = 0.1707; free.del_M_TL = 1; units.del_M_TL = '-'; label.del_M_TL = 'del_M_TL'; % shape coefficient for total length


% par.T_AL   = 50000;   free.T_AL   = 0;   units.T_AL = 'K';       label.T_AL = 'T_AL';      % Arrhenius temperature
% par.T_AH   = 90000;   free.T_AH   = 0;   units.T_AH = 'K';       label.T_AH = 'T_AH';      % Arrhenius temperature
% par.T_L   = T_C+0;   free.T_L   = 0;    units.T_L = 'K';        label.T_L = 'T_L';      % Arrhenius temperature
% par.T_H   = T_C+34;   free.T_H   = 0;    units.T_H = 'K';        label.T_H = 'T_H';      % Arrhenius temperature

% Environmental parameters (temperatures are in data)
par.f    = 1.0;     free.f     = 0; units.f    = '-';       label.f = 'f';         % scaled functional response for 0-var data
par.f_tL = 0.86;    free.f_tL  = 0; units.f_tL = '-';       label.f_tL = 'f_tL';         % scaled functional response for 0-var data

txt_par.units = units; txt_par.label = label; par.free = free; %pack units, label, free in  structure

%% set chemical parameters from Kooy2010 

% specific densities
d_V = get_d_V(metadata.phylum, metadata.class); d_E = d_V;

%
%       X     V     E     P       food, structure, reserve, product (=faeces)
d_O = [0.2;  d_V;  d_E;  0.2];    % g/cm^3, specific densities for organics

% chemical potentials from Kooy2010 Tab 4.2
%        X     V     E     P
mu_O = [525; 500;  550;  480] * 1000; % J/mol, chemical potentials for organics
%       C  H  O  N     CO2, H2O, O2, NH3
mu_M = [0; 0; 0; 0]; % chemical potential of minerals

% chemical indices from Kooy2010 Fig 4.15
%       X     V     E     P
n_O = [1.00, 1.00, 1.00, 1.00;  % C/C, equals 1 by definition
       1.80, 1.80, 1.80, 1.80;  % H/C  these values show that we consider dry-mass
       0.50, 0.50, 0.50, 0.50;  % O/C
       0.15, 0.15, 0.15, 0.15]; % N/C
%       C     H     O     N
n_M = [ 1     0     0     0;    % C/C, equals 0 or 1
        0     2     0     3;    % H/C
        2     1     2     0;    % O/C
        0     0     0     1];   % N/C

% molar volume of gas at 1 bar and 20 C is 24.4 L/mol
T = T_C + 20;            % K, temp of measurement equipment
X_gas = T_ref/ T/ 24.4;  % M, mol of gas per litre at T_ref (= 20 C) and 1 bar 

%% pack chem 
chem = struct(...
'T_ref', T_ref, ...
'mu_M', mu_M,  ...
'mu_X', mu_O(1),  ...
'mu_V', mu_O(2),  ...
'mu_E', mu_O(3),  ...
'mu_P', mu_O(4),  ...
'd_X',  d_O(1),   ...
'd_V',  d_O(2),   ...
'd_E',  d_O(3),   ...
'd_P',  d_O(4),   ...
'X_gas', X_gas, ...
'n_O',  n_O,   ...
'n_M',  n_M);

