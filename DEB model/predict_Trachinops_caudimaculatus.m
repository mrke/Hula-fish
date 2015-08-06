function [Prd_data, info] = predict_Engraulis_encrasicolus(par, chem, T_ref, data)

%% unpack par, chem, cpar and data
cpar = parscomp_st(par, chem);
v2struct(par); v2struct(chem); v2struct(cpar);
v2struct(data);
   
%% compute temperature correction factors
  TC_ab = tempcorr(temp.ab, T_ref, T_A);
  TC_aj = tempcorr(temp.aj, T_ref, T_A);
  TC_ap = tempcorr(temp.ap, T_ref, T_A);
  TC_am = tempcorr(temp.am, T_ref, T_A);
  TC_Ri = tempcorr(temp.Ri, T_ref, T_A);
  TC_tL = tempcorr(temp.tL, T_ref, T_A);

  % life cycle
  pars_tj = [g k l_T v_Hb v_Hj v_Hp];
  [t_j t_p t_b l_j l_p l_b l_i rho_j rho_B] = get_tj(pars_tj, f);

  % birth
  L_b = L_m * l_b;                  % cm, structural length at birth 
  Lw_b = L_b/ del_M;             % cm, standard length at birth
  Wd_b = L_b^3 * (1 + f * w);       % g, dry weight at birth
  aT_b = t_b/ k_M/ TC_ab;           % d, age at birth of foetus at f and T

  % metam
  L_j = L_m * l_j;                  % cm, structural length at metam
  Lw_j = L_j/ del_M;             % cm, standard length at metam
  Ww_j = L_j^3 * (1 + f * w);       % g, wet weight at metam
  aT_j = t_j/ k_M/ TC_aj;           % d, age at metam

  % puberty 
  L_p = L_m * l_p;                  % cm, structural length at puberty at f
  Lw_p = L_p/ del_M_TL;             % cm, total length at puberty at f
  Ww_p = L_p^3 * (1 + f * w);       % g, wet weight at puberty 
  aT_p = t_p/ k_M/ TC_ap;           % d, time since birth at puberty at f and T

  % ultimate
  L_i = L_m * l_i;                  % cm, ultimate structural length at f
  Lw_i = L_i/ del_M_TL;             % cm, ultimate total length at f
  Ww_i = L_i^3 * (1 + f * w);       % g, ultimate wet weight 

  % reproduction
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hj; U_Hp]; % compose parameter vector
  RT_i = TC_Ri * reprod_rate_j(L_i, f, pars_R);                 % ultimate reproduction rate

  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f, l_b);      % -, scaled mean life span at T_ref
  aT_m = t_m/ k_M/ TC_am;               % d, mean life span at T

  % pack to output
  Prd_data.ab = aT_b;
  Prd_data.aj = aT_j;
  Prd_data.ap = aT_p;
  Prd_data.am = aT_m;
  Prd_data.Lb = Lw_b;
  Prd_data.Lj = Lw_j;
  Prd_data.Lp = Lw_p;
  Prd_data.Li = Lw_i;
  Prd_data.Wb = Wd_b;
  Prd_data.Wj = Ww_j;
  Prd_data.Wp = Ww_p;
  Prd_data.Wi = Ww_i;
  Prd_data.Ri = RT_i;

  %% uni-variate data
  
  % 1st data set for late juveniles and adults
  [t_j t_p t_b l_j l_p l_b l_i rho_j rho_B info] = get_tj(pars_tj, f_tL);
  rT_B = TC_tL * rho_B * k_M; rT_j = TC_tL * rho_j * k_M;              % 1/d, von Bert, exponential growth rate
  aT_b = t_b/ k_M/ TC_tL; aT_j = t_j/ k_M/ TC_tL;
  t_j = aT_j - aT_b; % time since birth at metamorphosis
  t_bj = tL(tL(:,1) < t_j,1); % select times between birth & metamorphosis
  Lw_b = l_b * L_m/ del_M_TL; Lw_j = l_j * L_m/ del_M_TL; Lw_i = l_i * L_m/ del_M_TL;
  EL_bj = Lw_b * exp(t_bj * rT_j/3); % exponential growth as V1-morph
  t_ji = tL(tL(:,1) >= t_j,1); % selects times after metamorphosis
  EL_ji = Lw_i - (Lw_i - Lw_j) * exp( - rT_B * (t_ji - t_j)); % cm, expected length at time
  EL = [EL_bj; EL_ji]; % catenate lengths


  % length-weight
  EW = (LW(:,1) * del_M).^3 * (1 + f * w);                   % g, expected wet weight at time

  %% pack to output
  % the names of the fields in the structure must be the same as the data names in the mydata file
  Prd_data.tL = EL;
  Prd_data.LW = EW;