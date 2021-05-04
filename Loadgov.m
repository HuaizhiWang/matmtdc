function Pgov = Loadgov(casefile_dyn)

% Pgov = Loadgov(casefile_dyn)
% 
% Loads governor data
% 
% INPUTS
% casefile_dyn = m-file or struct with dynamic data
% 
% OUTPUTS
% Pgov = governor parameter matrix

% MatDyn
% Copyright (C) 2009 Stijn Cole
% Katholieke Universiteit Leuven
% Dept. Electrical Engineering (ESAT), Div. ELECTA
% Kasteelpark Arenberg 10
% 3001 Leuven-Heverlee, Belgium

%% Load data
if isstruct(casefile_dyn)
	gov = casefile_dyn.gov;
else
    [gen,exc,gov,mmc,DCVolCon,OutCurCon,CirCurCon,freq,stepsize,stoptime] = feval(casefile_dyn);
end

%% Consecutive numbering or rows
Pgov = gov;

for i = 1:length(gov(1,:))
    Pgov(gov(:,1),i) = gov(:,i);
end

return;