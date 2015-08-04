%% mydata_my_pet
% Sets referenced data

%%
function [data, txt_data, metadata] = mydata_Trachinops_caudimaculatus 
  % created by Starrlight Augustine, Bas Kooijman, Dina Lika, Goncalo Marques and Laure Pecquerie 2015/03/31
  % last modified: 07/07/2015 by starrlight
  
  %% Syntax
  % [data, txt_data, metadata] = <../mydata_my_pet.m *mydata_my_pet*>
  
  %% Description
  % Sets data, pseudodata, metadata, explanatory text, weight coefficients.
  % Meant to be a template in add_my_pet
  %
  % Output
  %
  % * data: structure with data
  % * txt_data: text vector for the presentation of results
  % * metadata: structure with info about this entry
  
  %% To do (remove these remarks after editing this file)
  % * copy this template; replace 'my_pet' by the name of your species
  % * fill in metadata fields with the proper information
  % * insert references for the data (an example is given), for multiple references, please use commas to separate references
  % * edit fact-list for your species, where you can add species and/or data properties
  % * edit real data; remove all data that to not belong to your pet
  % * complete reference list
  % * OPTIONAL : add discussion points / comments before the reference list

%% set metadata

T_C = 273.15; % K, temperature at 0 degrees C (used in T_typical)

metadata.phylum     = 'Chordata'; 
metadata.class      = 'Actinopterygii'; 
metadata.order      = 'Perciformes'; 
metadata.family     = 'Plesiopidae';
metadata.species    = 'Trachinops_caudimaculatus'; 
metadata.species_en = 'Southern Hulafish'; 
metadata.T_typical  = T_C + 20; % K
metadata.data_0     = {'ab'; 'ap'; 'am'; 'Lb'; 'Lp'; 'Li'; 'Wdb'; 'Wdp'; 'Wdi'; 'Ri'};  % tags for different types of zero-variate data
metadata.data_1     = {'L-W'}; % tags for different types of uni-variate data

metadata.COMPLETE = 2.5; % using criteria of LikaKear2011

metadata.author   = {'FirstName1 LastName1'};            % put names as authors as separate strings:  {'FirstName1 LastName2','FirstName2 LastName2'} , with corresponding author in first place 
metadata.date_subm = [2015 04 20];                       % [year month day], date at which the entry is submitted
metadata.email    = {'myname@myuniv.univ'};              % e-mail of corresponding author
metadata.address  = {'affiliation, zipcode, country'};   % affiliation, postcode, country of the corresponding author

% uncomment and fill in the following fields when the entry is updated:
% metadata.author_mod_1  = {'FirstName3 LastName3'};          % put names as authors as separate strings:  {'author1','author2'} , with corresponding author in first place 
% metadata.date_mod_1    = [2017 09 18];                      % [year month day], date modified entry is accepted into the collection
% metadata.email_mod_1   = {'myname@myuniv.univ'};            % e-mail of corresponding author
% metadata.address_mod_1 = {'affiliation, zipcode, country'}; % affiliation, postcode, country of the corresponding author

% for curators only ------------------------------
% metadata.curator     = {'FirstName LastName'};
% metadata.email_cur   = {'myname@myuniv.univ'}; 
% metadata.date_acc    = [2015 04 22]; 
%-------------------------------------------------

%% set data
% zero-variate data;
% typically depend on scaled functional response f.
% here assumed to be equal for all real data; the value of f is specified in pars_init_my_pet.

% age 0 is at onset of embryo development
data.ab = 15;      units.ab = 'd';    label.ab = 'age at birth';                bibkey.ab = 'MollCano2010';    
  temp.ab = T_C + 12;  % K, temperature 
  % observed age at birth is frequently larger than ab, because of diapauzes during incubation

% age 0 is at onset of embryo development
data.aj = 100;      units.aj = 'd';    label.aj = 'age at birth';                bibkey.aj = 'MollCano2010';    
  temp.aj = T_C + 13;  % K, temperature 
% observed age at birth is frequently larger than ab, because of diapauzes during incubation

data.ap = 200;     units.ap = 'd';    label.ap = 'age at puberty';              bibkey.ap = 'Anon2015';
  temp.ap = T_C + 15;  % K, temperature 
  % observed age at puberty is frequently larger than ap, 
  %   because allocation to reproduction starts before first eggs appear
data.am = 6*365;     units.am = 'd';    label.am = 'life span';                   bibkey.am = 'Wiki';   
  temp.am = T_C + 14;  % K, temperature 
% (accounting for aging only) 

% Please specify what type of length measurement is used for your species.
% We put here snout-to-vent length, but you should change this depending on your species:
% add an optional comment structure to give any additional explanations on
% how the value was chosen, see the last column of the Lb data set for an
% example
data.Lb  = 0.45;   units.Lb  = 'cm';   label.Lb  = 'snout to vent length at birth';    bibkey.Lb  = 'Anon2015'; comment.Lb  = 'mean value taken from several measurements';
data.Lp  = 3.8;   units.Lp  = 'cm';   label.Lp  = 'snout to vent length at puberty';  bibkey.Lp  = {'Anon2015','Wiki'}; % for multiple references, please use commas to separate references
data.Li  = 12.0;   units.Li  = 'cm';   label.Li  = 'ultimate snout to vent length';    bibkey.Li  = 'Wiki';
data.Wdb = 1.26e-5/5; units.Wdb = 'g';    label.Wdb = 'dry weight at birth';              bibkey.Wdb = 'Anon2015';
data.Wdp = .76/5;   units.Wdp = 'g';    label.Wdp = 'dry weight at puberty';            bibkey.Wdp = 'Anon2015';
data.Wdi = 15/5;   units.Wdi = 'g';    label.Wdi = 'ultimate dry weight';              bibkey.Wdi = 'Wiki';
data.Ri  = 14000/365;    units.Ri  = '#/d';  label.Ri  = 'maximum reprod rate';              bibkey.Ri  = 'Wiki';   
  % for an individual of ultimate length Li 
  temp.Ri = T_C + 14;  % K, temperature 
 
% uni-variate data

% uni-variate data at f = 0.8 (this value should be added in pars_init_my_pet as a parameter f_tL) 
% snout-to-vent length and wet weight were measured at the same time
% data.tL = [ 0     50  100 200 300 400 500 600;    % d, time since birth
%             0.45  1.1 1.7 2.7 3.4 4.0 4.5 4.9]';  % cm, snout-to-vent length at f and T
data.tL = csvread('uv_tL_Trachinops_caudimaculatus.csv');
data.tL = [data.tL(:,1), data.tL(:,2)/10];
units.tL = {'d', 'cm'};     label.tL = {'time since birth', 'snout to vent length'};  bibkey.tL = 'Anon2015';
  temp.tL = T_C + 14;  % K, temperature
comment.tL = 'Put here any remarks about the experimental protocol'; % optional field

%data.LW = [42.6	42.9	43.7	44.1	44.8	45	45.2	44.7	45.8	46.4	46.1	36.2	38.3	43.1	43.2	45.6	46.7	47.7	48.1	49.2	49.7	42.9	45.9	48.1	35.6	38.4	38.8	39.6	39.8	41.3	41.6	42.1	42.6	42.6	42.7	43	44.2	44.3	44.6	45	45.2	45.2	45.5	45.9	45.9	46.2	46.3	46.4	46.8	46.9	46.9	47.1	47.5	47.6	47.9	48.3	48.5	48.6	48.7	48.7	38.2	41.2	42.6	42.7	43.8	44.3	44.4	45	45.3	46.2	47.4	47.4	47.5	47.6	47.7	47.9	48.5	48.7	49.9	50	47.6	48.4	40.8	41.5	43	43.2	43.3	44	44.1	45.7	45.8	46	46.9	47.2	47.3	47.9	48.2	48.8	48.9	49.2	46	46.8	49.7	51.1	52.2	52.5	52.6	52.6	49.2	42.3	43.8	45.4	45.9	46.1	47	47.1	47.2	46.1	46.7	47.1	47.2	47.4	47.8	48	49.1	49.2	49.3	50.4	45.5	46.7	51.5	52.3	52.4	54.6	57.8	47.5	48.8	48.9	50.5	52.9	44.7	44.7	46.4	46.4	47.2	47.2	47.8	47.8	47.9	47.9	47.9	47.9	48.6	48.7	48.7	49	49	49.7	49.7	50.5	50.6	50.6	51.4	51.4	51.5	51.5	51.6	51.6	51.7	51.7	52.1	52.1	53.7	46	52.2	52.8	52.8	53.8	54	54.3	55.3	56.2	59.2	37.1	38.8	39.3	42.4	43.2	46	46.2	46.9	48.5	48.5	49.2	49.7	49.7	49.8	51.3	51.1	51.6	46.7	48.1	49	51.7	52.8	54.1	46.2	46.9	47	47.4	47.9	48.3	49	49.3	49.9	50.4	50.5	50.8	51.5	51.5	52.6	53.4	54.1	54.7	55.6	46.2	46.7	47	47.1	47.2	47.3	47.5	47.8	48.2	48.3	48.8	49.2	49.6	49.8	49.8	50.2	50.2	50.5	50.7	51.2	51.3	51.9	52.2	52.3	52.4	53	53.1	53.1	53.2	53.9	55.8	59.6	60.8	62.4	63.5	68.6	50.2	50.8	51.4	51.5	53	54.9	56.6	47.2	48.2	48.6	49.8	50.1	50.7	50.8	50.9	51	51	51.4	51.9	52	52	55.2	46.7	47.9	47	50.4	52.6	53.6	54.5	56	56.7	57.6	60.8	61.4	34.9	47.8	48	48.5	49.3	50.2	50.3	50.4	51.6	52.7	53	53.4	55.1	55.3	55.4	56.3	47.4	48.9	49.1	50.8	51.3	52.2	52.4	53.2	53.2	56.1	49	51.2	51.7	56.4	57	57.6	57.7	58.4	49	49.4	49.5	50.2	51.3	51.4	52	55.9	55.9	56.6	57.7	58.8	60.1	56.9	59.9	63	63.3	64.3	64.4	64.9	65.9	67.6	71.8	47.4	47.6	48.2	51	52.4	55.5	57.5	57.5	58.9	60.1	61	47.1	51.7	52.5	53.9	54.1	55.2	55.3	55.5	56.4	57.6	57.6	57.8	59.8	60.1	62	53.3	53.6	55.8	55.8	55.8	56.4	56.4	57.1	57.1	59.6	59.6	58.7	59	59	61.2	62	62.1	62.2	63.3	64	64.1	65	65.4	65.5	65.6	66.4	53.3	54.7	58.3	52.5	58.5	54.3	56	58.6	59	52.2	54.6	54.7	55.5	56	57	51.1	51.6	52.3	52.6	54	54.2	55	55.2	55.7	56.7	56.9	57.3	59	60.3	61	63.4	61.8	65	66.9	67.4	68.5	69.1	70.7	71.5	59.7	60.9	62.4	63.2	63.4	64.2	65.9	66.3	67.2	67.7	68.1	68.2	68.6	69.5	55	59.5	63.5	62.8	66.2	67.6	69.5	72.5	74.7	52.7	63.6	63.7	64.1	65.1	49.4	50.4	51.2	51.3	52.1	52.7	53.6	54	54.1	55	57.1	58.4	60.2	61.4	53.3	53.6	56.2	56.4	58.2	61.8	61	61.1	61.4	62.7	64.7	65.7	69.9	67.8	69.1	70.9	56.7	61.9	68.6	63.7	66.9	67.7	68.8	69.1	70.1	70.8	71.3	71.4	54	58.1	60.6	54.8	68.7	56.8	57.6	59.1	60.4	60.4	60.9	60.9	61.5	61.5	63.8	67	69.3	70	70.3	71.3	73.3	73.5	64.6	64.1	59	63.1	68.2	73	73.3	75.3	76.7	77.4	77.4	77.5	71.4	71.8	73.1	71.3	72.4	78.2	60.3	69.1	71.4	72.9	73.5	61.6	66.2	70.1	72	72.5	74.3	78	81.8	84.2	85.4	70.6	73.3	74.3	76	63.8	66.9	66.9	80.1	83.9	76.9	77.7	80	79	80.8	89.3	79.8;      % cm, snout-to-vent length at f
%            1.01	1.12	1.07	1.11	1.17	1.19	1.12	1.2	1.19	1.35	1.44	0.57	0.78	0.91	0.9	1.03	1.11	1.2	1.42	1.08	1.41	0.96	1.1	1.19	0.48	0.6	0.58	0.71	0.75	0.93	0.78	0.77	0.67	0.94	0.78	0.73	1.05	0.91	0.84	0.83	0.95	1.03	1.01	1.1	1.01	1.11	1.03	1.36	0.99	0.99	1	1.12	1.15	1.24	1.22	1.38	1.13	1.2	0.84	1.18	0.79	1.03	0.85	0.95	1.16	1.09	0.85	1.22	1.2	1.23	1.43	1.45	1.15	1.19	1.24	1.45	1.21	1.49	1.73	1.14	1.11	1.46	0.89	0.94	0.96	1.17	0.98	1.14	1.16	1.22	1.49	1.13	1.42	1.47	1.31	1.5	1.48	1.52	1.4	1.35	1.43	1.37	1.72	1.84	2.06	1.99	1.89	2.1	1.5	1.09	1.15	1.27	1.31	1.32	1.35	1.35	1.35	1.18	1.25	1.34	1.34	1.27	1.31	1.38	1.25	1.37	1.46	1.52	1.42	1.37	1.91	2.01	2.26	2.2	2.48	1.58	1.57	1.64	1.88	2.11	1.18	1.18	1.27	1.27	1.55	1.55	1.55	1.55	1.46	1.58	1.46	1.58	1.47	1.45	1.45	1.65	1.65	1.65	1.65	1.72	2.06	2.06	1.7	1.7	1.93	1.93	1.82	1.82	2.03	2.03	1.95	1.95	1.77	1.29	1.85	1.73	1.89	2.19	2	2.11	1.98	2.2	2.71	0.7	0.8	0.88	1.01	1.1	1.35	1.37	1.31	1.51	1.52	1.69	1.49	1.51	1.82	1.62	1.73	1.8	1.46	1.37	1.59	1.87	2.27	2.21	1.38	1.44	1.57	1.43	1.54	1.2	1.45	1.6	1.9	1.46	1.67	1.96	1.97	2.04	1.98	2.13	2.04	1.89	2.24	1.52	1.45	1.35	1.58	1.48	1.25	1.59	1.41	1.53	1.45	1.66	1.59	1.42	1.5	1.76	1.68	1.46	1.68	1.53	1.6	1.6	1.58	2.08	1.92	1.71	2.21	1.73	1.81	1.89	2.02	2.22	2.35	2.75	3.3	3.1	3.74	1.93	1.96	1.8	1.77	1.62	2.1	2.42	1.61	1.6	1.75	1.76	1.92	2.15	1.96	2.16	1.97	1.99	1.99	2.01	2.01	2.09	2.44	1	1.06	1.23	1.38	1.37	1.74	1.86	2.06	1.96	1.73	1.33	2.71	0.4	1.35	1.02	1.22	1.37	1.27	1.27	1.21	1.44	1.47	1.46	1.8	1.78	1.39	1.55	1.98	1.46	1.51	1.51	1.8	1.56	2.04	1.62	1.5	1.6	2.1	1.55	1.45	1.88	2.05	1.95	1.94	2.2	2.36	1.47	1.62	1.66	1.62	1.69	1.75	1.87	2.48	2.4	2.74	2.37	2.54	2.43	2.11	2.55	2.81	2.93	3.41	3.4	3.12	3.45	3.67	4.13	1.46	1.49	1.39	1.59	1.77	1.9	2.89	2.46	2.61	3.17	3.4	1.35	2	2.01	2.06	2.15	2.26	2.45	2.41	2.41	2.68	2.67	2.25	2.56	2.67	3.11	1.91	1.82	2.41	2.41	1.95	2.42	2.42	2.59	2.59	2.7	2.7	2.44	2.39	2.51	2.81	3.02	3.16	3.08	2.71	3.47	3.2	3.3	3.48	3.94	3.24	3.81	1.9	1.89	2.44	1.88	2.29	2.24	2.47	2.72	2.64	1.85	1.93	2.09	2.15	2.24	2.43	1.94	1.8	1.9	1.75	2.32	1.85	1.76	1.74	2.25	2.53	2.67	2.63	2.85	3.07	3.13	3.31	2.68	2.67	3.39	4.06	3.41	3.71	4.64	4.51	2.78	2.94	3.06	3.25	3.52	3.4	3.6	3.64	4	4.11	3.91	4.15	4.35	4.25	1.92	3.11	3.08	3.36	4.03	3.78	4.19	4.56	6.65	1.57	2.12	2.59	2.88	2.79	1.32	1.44	1.41	1.05	1.54	1.45	1.46	1.17	1.84	1.76	1.81	1.76	2.36	2.05	1.39	1.87	1.9	1.8	1.78	2.38	2.23	2.7	2.95	2.64	3.69	3.25	3.92	3.38	3.16	3.78	2.24	2.88	3.71	3.06	3.85	3.91	3.76	4.52	3.91	4.46	3.99	4.06	1.75	1.82	2.4	2.42	3.8	2.2	2	2.43	3.05	3.05	3.05	3.05	3.16	3.16	3.35	3.62	3.93	4	3.52	4.68	4.71	4.67	2.98	3.51	2.55	2.93	4.27	4.52	4.83	5.14	6.02	5.5	6.61	5.6	4.69	4.85	5.13	3.94	4.63	5.75	2.03	3.43	3.04	4.25	3.76	2.45	3.11	4.06	3.93	3.2	5.2	5.14	4.53	5.39	6.3	4.58	4.19	4.4	4.83	3.35	3.38	3.38	5.69	6.98	4.94	3.75	4.17	5.18	5.45	6.36	5.45]';   % g, wet weight at f and T
data.LW = csvread('uv_LW_Trachinops_caudimaculatus.csv');
data.LW = [data.LW(:,1)/10, data.LW(:,2)];
units.LW = {'cm', 'g'};     label.LW = {'snout to vent length', 'wet weight'};  bibkey.LW = 'Anon2015';
comment.LW = 'Put here any remarks about the experimental protocol'; % optional field

%% set weights for all real data
weight = setweights(data, []);

%% overwriting weights (remove these remarks after editing the file)
% the weights were set automatically with the function setweigths,
% if one wants to ovewrite one of the weights it should always present an explanation example:
%
% zero-variate data:
% weight.Wdi = 100 * weight.Wdi; % Much more confidence in the ultimate dry
%                                % weight than the other data points
% uni-variate data: 
% weight.tL = 2 * weight.tL;
%weight.Wdb = 0;
weight.tL = 0*weight.tL;
weight.LW = 0*weight.LW;

%% set pseudodata and respective weights
% (pseudo data are in data.psd and weights are in weight.psd)
[data, units, label, weight] = addpseudodata(data, units, label, weight);

%% overwriting pseudodata and respective weights (remove these remarks after editing the file)
% the pseudodata and respective weights were set automatically with the function setpseudodata
% if one wants to overwrite one of the values it should always present an explanation
% example:
% data.psd.p_M = 1000;                    % my_pet belongs to a group with high somatic maint 
% weight.psd.kap = 10 * weight.psd.kap;   % I need to give this pseudo data a higher weight

%% pack data and txt_data for output
data.weight = weight;
data.temp = temp;
txt_data.units = units;
txt_data.label = label;
txt_data.bibkey = bibkey;
if exist('comment','var')
  txt_data.comment = comment;
end

%% Facts
% list facts: F1, F2, etc.
% make sure each fact has a corresponding bib key
% do not put any DEB modelling assumptions here, only relevant information on
% biology and life-cycles etc.

F1 = 'The larval stage last 202 days and no feeding occurs';
metadata.bibkey.F1 = 'Wiki'; % optional bibkey
metadata.facts = struct('F1',F1);

%% Discussion points
D1 = 'Author_mod_1: I found information on the number of eggs per female as a function of length in Anon2013 that was much higher than in Anon2015 but chose to not include it as the temperature was not provided';
% optional bibkey: bibkey.D1 = 'Anon2013';
D2 = 'Author_mod_1: I was surprised to observe that the weight coefficient for ab changed so much the parameter values';     
% optional bibkey: bibkey.D2 = 'Kooy2010';
metadata.discussion = struct('D1', D1, 'D2', D2);

%% References
  % the following two references should be kept-----------------------------------------------------------
  bibkey = 'Kooy2010'; type = 'Book'; bib = [ ...  % used in setting of chemical parameters and pseudodata
  'author = {Kooijman, S.A.L.M.}, ' ...
  'year = {2010}, ' ...
  'title  = {Dynamic Energy Budget theory for metabolic organisation}, ' ...
  'publisher = {Cambridge Univ. Press, Cambridge}, ' ...
  'pages = {Table 4.2 (page 150), 8.1 (page 300)}, ' ...
  'howpublished = {\url{http://www.bio.vu.nl/thb/research/bib/Kooy2010.html}}'];
  eval(['metadata.biblist.' bibkey, '= ''@', type, '{', bibkey, ', ' bib, '}'';']);
  %
  bibkey = 'LikaKear2011'; type = 'Article'; bib = [ ...  % used for the estimation method
   'author = {Lika, K. and Kearney, M.R. and Freitas, V. and van der Veer, H.W. and van der Meer, J. and Wijsman, J.W.M. and Pecquerie, L. and Kooijman, S.A.L.M.},'...
   'year = {2011},'...
   'title = {The ''''covariation method'''' for estimating the parameters of the standard Dynamic Energy Budget model \textrm{I}: Philosophy and approach},'...
   'journal = {Journal of Sea Research},'...
   'volume = {66},'...
   'number = {4},'...
   'pages = {270-277},'...
   'DOI = {10.1016/j.seares.2011.07.010},'...
   'howpublished = {\url{http://www.sciencedirect.com/science/article/pii/S1385110111001055}}'];
  eval(['metadata.biblist.' bibkey, '= ''@', type, '{', bibkey, ', ' bib, '}'';']);
  %------------------------------------------------------------------------------------------------------

  % References for the data, following BibTex rules
 % author names : author = {Last Name, F. and Last Name2, F2. and Last Name 3, F3. and Last Name 4, F4.}
 % latin names in title e.g. \emph{Pleurobrachia pileus}

  bibkey = 'Wiki'; type = 'Misc'; bib = [...
  'howpublished = {\url{http://en.wikipedia.org/wiki/my_pet}},'...% replace my_pet by latin species name
  'note = {Accessed : 2015-04-30}'];
  eval(['metadata.biblist.' bibkey, '= ''@', type, '{', bibkey, ', ' bib, '}'';']);
  %
  bibkey = 'MollCano2010'; type = 'Article'; bib = [ ... % meant as example; replace this and further bib entries
  'author = {M{\o}ller, L. F. and Canon, J. M. and Tiselius, P.}, ' ... 
  'year = {2010}, ' ...
  'title = {Bioenergetics and growth in the ctenophore \emph{Pleurobrachia pileus}}, ' ...
  'journal = {Hydrobiologia}, ' ...
  'volume = {645}, ' ...
  'number = {4}, '...
  'pages = {167-178}'];
  eval(['metadata.biblist.' bibkey, '= ''@', type, '{', bibkey, ', ' bib, '}'';']);
  %
  bibkey = 'Anon2015'; type = 'Misc'; bib = [ ...
  'author = {Anonymous}, ' ...
  'year = {2015}, ' ...
  'howpublished = {\url{http://www.fishbase.org/summary/Rhincodon-typus.html}}'];
  eval(['metadata.biblist.' bibkey, '= ''@', type, '{', bibkey, ', ' bib, '}'';']);


