function [event,buschange,linechange] = fault

% fault
% MatDyn event data file
% 
% MatDyn
% Copyright (C) 2009 Stijn Cole
% Katholieke Universiteit Leuven
% Dept. Electrical Engineering (ESAT), Div. ELECTA
% Kasteelpark Arenberg 10
% 3001 Leuven-Heverlee, Belgium

%%
% event = [time type]
event=[1.1   1; 
       1.15    1];  

% buschange = [time bus(row)  attribute(col) new_value]
buschange   = [1.1  9   6   -1e10;
               1.15  9   6   0];

% linechange = [time  line(row)  attribute(col) new_value]
linechange = [];

return;