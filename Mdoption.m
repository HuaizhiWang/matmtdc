function options = Mdoption

% options = Mdoption
% 
% Returns default option vector

% MatDyn
% Copyright (C) 2009 Stijn Cole
% Katholieke Universiteit Leuven
% Dept. Electrical Engineering (ESAT), Div. ELECTA
% Kasteelpark Arenberg 10
% 3001 Leuven-Heverlee, Belgium

%% Options vector
options = [ 
    1;              % method
    1e-4;           % tolerance
    1e-3;           % minstepsize
    1e2;            % maxstepsize
    
    1;              % output
    1;              % plots
];

return;