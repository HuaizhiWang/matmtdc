function Pexc = Loadexc(casefile_dyn)

% Pexc = Loadexc(casefile_dyn)
% 
% Loads exciter data
% 
% INPUTS
% casefile_dyn = m-file or struct with dynamic data
% 
% OUTPUTS
% Pexc = exciter parameter matrix

% MatDyn
% Copyright (C) 2009 Stijn Cole
% Katholieke Universiteit Leuven
% Dept. Electrical Engineering (ESAT), Div. ELECTA
% Kasteelpark Arenberg 10
% 3001 Leuven-Heverlee, Belgium

%% Load data
if isstruct(casefile_dyn)
	exc = casefile_dyn.exc;
else
    [gen,exc,gov,mmc,DCVolCon,OutCurCon,CirCurCon,freq,stepsize,stoptime] = feval(casefile_dyn);
end

%% Consecutive numbering or rows
Pexc = exc;

for i = 1:length(exc(1,:))
    Pexc(exc(:,1),i) = exc(:,i);
end

return;