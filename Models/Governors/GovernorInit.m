function [Xgov0,Pgov0] = GovernorInit(Xgov, Pgov, Vgov, govtype)

% F = Governor(Xgov, Pgov, Vgov, govtype)
%
% Calculate initial conditions governors
 
% MatDyn
% Copyright (C) 2009 Stijn Cole
% Katholieke Universiteit Leuven
% Dept. Electrical Engineering (ESAT), Div. ELECTA
% Kasteelpark Arenberg 10
% 3001 Leuven-Heverlee, Belgium

%% Init
global freq;

[ngen,c] = size(Xgov);
Xgov0 = zeros(ngen,c);
[ngen,c] = size(Pgov);
Pgov0 = zeros(ngen,c+2);
d=[1:length(govtype)]';

%% Define types
type1 = d(govtype==1);
type2 = d(govtype==2);

%% Governor type 1: constant power
Pm0 = Xgov(type1,1);
Xgov0(type1,1) = Pm0;
  
%% Governor type 2: IEEE general speed-governing system
Pm0 = Xgov(type2,1);

K = Pgov(type2,2);
T1 = Pgov(type2,3);
T2 = Pgov(type2,4);
T3 = Pgov(type2,5);
Pup = Pgov(type2,6);
Pdown = Pgov(type2,7);
Pmax = Pgov(type2,8);
Pmin = Pgov(type2,9);

omega0 = Vgov(type2,1);
        
zz0 = Pm0;
PP0 = Pm0;


P0 = K.*(2*pi*freq - omega0);
xx0 = T1.*(1-T2./T1).*(2*pi*freq - omega0);

Xgov0(type2,1:4) = [Pm0 P0 xx0 zz0];
Pgov0(type2,1:10) = [Pgov(type2,1), K, T1, T2, T3, Pup, Pdown, Pmax, Pmin, PP0];  

%% Governor type 3: 

%% Governor type 4: 

return;