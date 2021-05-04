function Pgen = Loadgen(casefile_dyn, output)

% Pgen = Loadgen(casefile_dyn)
% 
% Loads generator data
% 
% INPUTS
% casefile_dyn = m-file or struct with dynamic data
% output = 1: print progress info, 0: do not print progress info 
% 
% OUTPUTS
% Pgen = generator parameter matrix

% MatDyn
% Copyright (C) 2009 Stijn Cole
% Katholieke Universiteit Leuven
% Dept. Electrical Engineering (ESAT), Div. ELECTA
% Kasteelpark Arenberg 10
% 3001 Leuven-Heverlee, Belgium

%% Load data
if isstruct(casefile_dyn)
    gen = casefile_dyn.gen;
else
    [gen,exc,gov,mmc,DCVolCon,OutCurCon,CirCurCon,freq,stepsize,stoptime] = feval(casefile_dyn);
end

Pgen = gen;

genmodel = Pgen(:,1);


%% Define generator models
d=[1:length(genmodel)]';
type1 = d(genmodel==1);
type2 = d(genmodel==2);


%% Check transient saliency
xd_tr = Pgen(type2,8);
xq_tr = Pgen(type2,9);

if sum(xd_tr~=xq_tr)>=1
    if output; fprintf('> Warning: transient saliency not supported\n                                     '); end
    Pgen(type2,9) = Pgen(type2,8);
end

return;