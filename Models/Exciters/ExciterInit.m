function [Xexc0,Pexc0] = ExciterInit(Xexc, Pexc, Vexc, exctype)

% [Xexc0,Pexc0] = ExciterInit(Xexc, Pexc, Vexc, exctype)
%
% Calculate initial conditions exciters
 
% MatDyn
% Copyright (C) 2009 Stijn Cole
% Katholieke Universiteit Leuven
% Dept. Electrical Engineering (ESAT), Div. ELECTA
% Kasteelpark Arenberg 10
% 3001 Leuven-Heverlee, Belgium

%% Init
[ngen,c] = size(Xexc);
Xexc0 = zeros(ngen,c);
[ngen,c] = size(Pexc);
Pexc0 = zeros(ngen,c+2);
d=[1:length(exctype)]';

%% Define types
type1 = d(exctype==1);
type2 = d(exctype==2);

%% Exciter type 1: constant excitation
Efd0 = Xexc(type1,1);
Xexc0(type1,1) = Efd0;
  
%% Exciter type 2: IEEE DC1A 
Efd0 = Xexc(type2,1);
Ka = Pexc(type2,2);
Ta = Pexc(type2,3);
Ke = Pexc(type2,4);
Te = Pexc(type2,5);
Kf = Pexc(type2,6);
Tf = Pexc(type2,7);
Aex = Pexc(type2,8);
Bex = Pexc(type2,9);
Ur_min = Pexc(type2,10);
Ur_max = Pexc(type2,11);
        
U = Vexc(type2,1);
        
Uf = zeros(length(type2),1);
Ux = Aex.*exp(Bex.*Efd0);
Ur = Ux+Ke.*Efd0;
Uref2 = U + (Ux+Ke.*Efd0)./Ka - U;
Uref = U;

Xexc0(type2,1:3) = [Efd0, Uf, Ur];   
Pexc0(type2,1:13) = [Pexc(type2,1), Ka, Ta, Ke, Te, Kf, Tf, Aex, Bex, Ur_min, Ur_max, Uref, Uref2];   

%% Exciter type 3

%% Exciter type 4

return;