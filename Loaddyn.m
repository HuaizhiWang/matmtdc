function [freq,stepsize,stoptime]=Loaddyn(casefile_dyn)

% [freq,stepsize,stoptime]=Loaddyn(casefile_dyn)
% 
% Loads dynamic data
% 
% INPUTS
% casefile_dyn = m-file with dynamic data
% 
% OUTPUTS
% freq = system frequency
% stepsize = step size for fixed-step integration algorithms
% stoptime = stop time of integration algorithm

% MatDyn
% Copyright (C) 2009 Stijn Cole
% Katholieke Universiteit Leuven
% Dept. Electrical Engineering (ESAT), Div. ELECTA
% Kasteelpark Arenberg 10
% 3001 Leuven-Heverlee, Belgium

%%
if isstruct(casefile_dyn)
	freq = casefile_dyn.freq;
    stepsize = casefile_dyn.stepsize;
    stoptime = casefile_dyn.stoptime;
else
    [gen,exc,gov,mmc,DCVolCon,OutCurCon,CirCurCon,freq,stepsize,stoptime] = feval(casefile_dyn);
end

return;