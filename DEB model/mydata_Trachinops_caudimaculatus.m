  function [data, txt_data, metadata] = mydata_Engraulis_encrasicolus

T_C = 273.15; % K, temperature at 0 degrees C

%% set other metadata

metadata.phylum     = 'Chordata'; 
metadata.class      = 'Actinopterygii'; 
metadata.order      = 'Perciformes'; 
metadata.family     = 'Plesiopidae';
metadata.species    = 'Trachinops_caudimaculatus'; 
metadata.species_en = 'Southern Hulafish';
metadata.T_typical  = T_C + 15; % K
metadata.data_0     = {'ab'; 'aj'; 'ap'; 'am'; 'Lb'; 'Lj'; 'Lp'; 'Li'; 'Wdb'; 'Wwj'; 'Wwp'; 'Wwi'; 'Ri'};  % tags for different types of zero-variate data
metadata.data_1     = {'t-L'; 'L-W'}; % tags for different types of uni-variate data

metadata.COMPLETE = 2.7; % using criteria of LikaKear2011

metadata.author   = {'Seann Chia'; 'Steve Swearer'; 'Mike Kearney'};        
metadata.date_acc = [2015 01 01];                           
metadata.email    = {'seann.chia@unimelb.edu.au'};                 
metadata.address  = {'University of Melbourne, Australia'}; 

%% set data
% zero-variate data
data.ab = 17;      units.ab = 'd';    label.ab = 'age at birth';           bibkey.ab = 'GarrSaiz2012'; 
  temp.ab = T_C + 12;  % K, temperature 
data.aj = 60;     units.aj = 'd';    label.aj = 'age at metam';           bibkey.aj = 'Re1996'; 
  temp.aj = T_C + 13;  % K, temperature 
data.ap = .8*365; units.ap = 'd';    label.ap = 'age at puberty';         bibkey.ap = 'PecqPeti2009';
  temp.ap = T_C + 15;  % K, temperature: average sea surface temperature off West Portuguese coast, Angélico, pers. comm.)\
data.am = 6*365;  units.am = 'd';    label.am = 'life span';              bibkey.am = 'Pecq2008';   
  temp.am = T_C + 14;  % K, temperature
  % average life span: 7-8 years, though max. age observed is 13 years old

data.Lb = 0.45;    units.Lb = 'cm';   label.Lb = 'standard length at birth';  bibkey.Lb = 'CataFolk2010';
% mid-range of smaller size class with gut content (p.309)
data.Lj = 2;   units.Lj = 'cm';   label.Lj = 'standard length at metam';bibkey.Lj = 'Re1996';
data.Lp = 3.8;      units.Lp = 'cm';   label.Lp = 'standard length at puberty';bibkey.Lp = 'PecqPeti2009';
data.Li = 12;     units.Li = 'cm';   label.Li = 'ultimate total length';  bibkey.Li = 'Pecq2008';
% Average length for Age 3 individuals in spring 2000, Ifremer PELGAS surveys (Pecq2008 p. 63)

data.Wb = (0.0009*(data.Lb*10)^2 - 0.021*(data.Lb*10) + 0.1405)/5; units.Wb = 'g';    label.Wb = 'dry weight at birth';    bibkey.Wb = 'CataFolk2010';
% (CataFolk2010) ln DW (\mug) = 0.867 + 0.5 SL (mm) - field data
data.Wj = (0.0014*(data.Lj*10)^2 - 0.048*(data.Lj*10) + 0.4953)/5;  units.Wj = 'g';    label.Wj = 'dry weight at metam';    bibkey.Wj = 'Pecq2008';
% from Ifremer Juvesu 1999 survey, 3 - 4.4 cm SL fish average
data.Wp = 0.0014*(data.Lp*10)^2 - 0.048*(data.Lp*10) + 0.4953; units.Wp = 'g';    label.Wp = 'wet weight at puberty';  bibkey.Wp = 'Pecq2008';
% From Ifremer PELGAS surveys (PecqPeti2009) - Total length (need checking, includes reprod buffer already) 
data.Wi = 0.0014*(data.Li*10)^2 - 0.048*(data.Li*10) + 0.4953; units.Wi = 'g';    label.Wi = 'ultimate wet weight';    bibkey.Wi = 'Pecq2008';
% From Ifremer PELGAS surveys (PecqPeti2009) - Total length (need checking)
data.Ri = (0.0014*(data.Li*10)^3.3714)/365;   units.Ri = '#/d';  label.Ri = 'maximum reprod rate';    bibkey.Ri = {'Pecq2008'};   
  temp.Ri = T_C + 14;  % K, temperature
% 20 batches of eggs ; between 200 and 500 eggs/g gonad-free wet
% weight/batch ; GSI  = 0.04 (PethRoos2013)
% Ri = 20 * 400 * 0.0041 * Li ^3.2099 * 0.96 / 365= 922
  
% uni-variate data
% Age-Weight
data.tL = csvread('uv_tL_Trachinops_caudimaculatus.csv');
data.tL = [data.tL(:,1), data.tL(:,2)/10];
units.tL = {'d', 'cm'};     label.tL = {'time since birth', 'standard length'};  bibkey.tL = 'Anon2015';
  temp.tL = T_C + 14;  % K, temperature
comment.tL = 'Put here any remarks about the experimental protocol'; % optional field

% Length-Weight
data.LW = csvread('uv_LW_Trachinops_caudimaculatus.csv');
data.LW = [data.LW(:,1), data.LW(:,2)];
units.LW = {'cm', 'g'};     label.LW = {'standard length', 'wet weight'};  bibkey.LW = 'Anon2015';
comment.LW = 'Put here any remarks about the experimental protocol'; % optional field

%% set weights for all real data
weight = setweights(data, []);

%% set pseudodata and respective weights
[data, units, label, weight] = addpseudodata(data, units, label, weight);

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


