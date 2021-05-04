function F = Generator(Xgen, Xexc, Xgov, Pgen, Vgen, gentype)

% F = Generator(Xgen, Xexc, Xgov, Pgen, Vgen, gentype)
%
% Generator model
 
% MatDyn
% Stijn Cole
% Katholieke Universiteit Leuven
% Dept. Electrical Engineering (ESAT), Div. ELECTA
% Kasteelpark Arenberg 10
% 3001 Leuven-Heverlee, Belgium

%% Init
global freq;
omegas=2*pi*freq;

[r,c] = size(Xgen);
F = zeros(r,c);
d=[1:length(gentype)]';

%% Define generator types
type1 = d(gentype==1);
type2 = d(gentype==2);

%% Generator type 1: classical model
omega = Xgen(type1,2);
Pm0 = Xgov(type1,1);

H = Pgen(type1,4);
D = Pgen(type1,5);

Pe = Vgen(type1,3);

ddelta = omega - omegas;
domega = pi .* freq ./ H .* (-D.*(omega - omegas) + Pm0 - Pe);
dEq = zeros(length(type1),1);

F(type1,1:3) = [ddelta domega dEq];

%% Generator type 2: 4th order model
omega = Xgen(type2,2);
Eq_tr = Xgen(type2,3);
Ed_tr = Xgen(type2,4);

H = Pgen(type2,4);
D = Pgen(type2,5);
xd = Pgen(type2,6);
xq = Pgen(type2,7);
xd_tr = Pgen(type2,8);
xq_tr = Pgen(type2,9);
Td0_tr = Pgen(type2,10);
Tq0_tr = Pgen(type2,11);

Id = Vgen(type2,1);
Iq = Vgen(type2,2);
Pe = Vgen(type2,3);

Efd = Xexc(type2,1);
Pm = Xgov(type2,1);

ddelta = omega - omegas;
domega = pi .* freq ./ H .* (-D.*(omega - omegas) + Pm - Pe);
dEq = 1./Td0_tr .* (Efd - Eq_tr + (xd - xd_tr).*Id);
dEd = 1./Tq0_tr .* (-Ed_tr - (xq - xq_tr).*Iq);

F(type2,1:4) = [ddelta domega dEq dEd];

%% Generator type 3:

%% Generator type 4:

return;