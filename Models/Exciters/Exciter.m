function F = Exciter(Xexc, Pexc, Vexc, exctype)

% F = Exciter(Xexc, Pexc, Vexc, exctype)
%
% Exciter model
 
% MatDyn
% Copyright (C) 2009 Stijn Cole
% Katholieke Universiteit Leuven
% Dept. Electrical Engineering (ESAT), Div. ELECTA
% Kasteelpark Arenberg 10
% 3001 Leuven-Heverlee, Belgium

%% Init
[r,c] = size(Xexc);
F = zeros(r,c);
d=[1:length(exctype)]';

%% Define exciter types
type1 = d(exctype==1);
type2 = d(exctype==2);

%% Exciter type 1: constant excitation  
F(type1,1) = 0;

%% Exciter type 2: IEEE DC1A
Efd = Xexc(type2,1);
Uf = Xexc(type2,2);
Ur = Xexc(type2,3);
        
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
Uref = Pexc(type2,12);
Uref2 = Pexc(type2,13);

U = Vexc(type2,1);
        
Ux = Aex.*exp(Bex.*Efd);
dUr = 1./Ta .* (Ka.*(Uref - U + Uref2 - Uf) - Ur);
dUf = 1./Tf .* (Kf./Te .* (Ur - Ux - Ke.*Efd) - Uf);
if sum(Ur>Ur_max)>=1
	Ur2=Ur_max;
elseif sum(Ur<Ur_min)>=1
	Ur2=Ur_min;
else
	Ur2=Ur;
end
dEfd = 1./Te .* ( Ur2 - Ux - Ke.*Efd);
        
F(type2,1:3) = [dEfd dUf dUr];

%% Exciter type 3: 

%% Exciter type 4:

return;