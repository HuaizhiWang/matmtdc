function F = Governor(Xgov, Pgov, Vgov, govtype)

% F = Governor(Xgov, Pgov, Vgov, govtype)
%
% Governor model
 
% MatDyn
% Copyright (C) 2009 Stijn Cole
% Katholieke Universiteit Leuven
% Dept. Electrical Engineering (ESAT), Div. ELECTA
% Kasteelpark Arenberg 10
% 3001 Leuven-Heverlee, Belgium

%% Init
global freq;
omegas=2*pi*freq;

[r,c] = size(Xgov);
F = zeros(r,c);
d=[1:length(govtype)]';

%% Define governor types
type1 = d(govtype==1);
type2 = d(govtype==2);

%% Governor type 1: constant power  
F(type1,1) = 0;

%% Governor type 2: IEEE general speed-governing system
Pm = Xgov(type2,1);
P = Xgov(type2,2);
x = Xgov(type2,3);
z = Xgov(type2,4);

K = Pgov(type2,2);
T1 = Pgov(type2,3);
T2 = Pgov(type2,4);
T3 = Pgov(type2,5);
Pup = Pgov(type2,6);
Pdown = Pgov(type2,7);
Pmax = Pgov(type2,8);
Pmin = Pgov(type2,9);
P0 = Pgov(type2,10);
        
omega = Vgov(type2,1);
        
dx = K.*(-1./T1.*x + (1 - T2./T1).*(omega - omegas));
dP = 1./T1.*x + T2./T1.*(omega - omegas);

y = 1./T3.*(P0 - P - Pm);

y2 = y;

if sum(y>Pup)>=1
    y2 = (1 - (y>Pup)).*y2 + (y>Pup).*Pup;
end
if sum(y<Pdown)>=1
    y2 = (1 - (y<Pdown)).*y2 + (y<Pdown).*Pdown;
end

dz = y2;

dPm = y2;

if sum(z>Pmax)>=1
    dPm = (1 - (z>Pmax)).*dPm + (z>Pmax).*0;
end
if sum(z<Pmin)>=1
    dPm = (1 - (z<Pmin)).*dPm + (z<Pmin).*0;
end
        
F(type2,1:4) = [dPm dP dx dz];

%% Governor type 3

%% Governor type 4

return;