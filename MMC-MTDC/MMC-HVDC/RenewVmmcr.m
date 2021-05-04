function [Vmmcr1, U2] = RenewVmmcr(Xmmci1, Xmmcr1, Pmmc0, Vmmcr0, t)

% [U] = SolveNetwork(Xgen, Pgen, Ly ,Uy ,Py ,gbus, gentype)
% 
% Solve the network
% 
% INPUTS
% Xgen = state variables of generators
% Pgen = parameters of generators
% Ly ,Uy ,Py = factorised augmented bus admittance matrix
% gbus = generator buses
% gentype = generator models
% 
% OUTPUTS
% U = bus voltages
 
% MatDyn
% Copyright (C) 2009 Stijn Cole
% Katholieke Universiteit Leuven
% Dept. Electrical Engineering (ESAT), Div. ELECTA
% Kasteelpark Arenberg 10
% 3001 Leuven-Heverlee, Belgium

%% Init
% Pmmc0=[N M Vd Vsmax Ismax phi Srated w1 L R C Carm Cd Rdc]
Vsmax = Pmmc0(1,4);
w1 = Pmmc0(1,8);
Rdc = Pmmc0(1,14);

% Xmmc0 = [ic0; vcu0; vcl0; Is_alpha0; Is_beta0; Vd0]
Is_alphar = Xmmcr1(4,1)*1.5*Vsmax/100e6;
Is_betar = Xmmcr1(5,1)*1.5*Vsmax/100e6;
Vdr = Xmmcr1(6,1);
Vdi = Xmmci1(6,1);

% Vmmc0 = [Va_alpha; Va_beta; is0; Ps; Qs; id0]
Va_alphar = Vmmcr0(1,1)/Vsmax;
Va_betar = Vmmcr0(2,1)/Vsmax;
% Power output
Pr = real((Va_alphar+1j*Va_betar)*conj(Is_alphar+1j*Is_betar)); % in pu
Qr = imag((Va_alphar+1j*Va_betar)*conj(Is_alphar+1j*Is_betar)); % in pu
% transform the power into equivalent impedance
Zeq = -(abs(Va_alphar+1j*Va_betar)^2)/(Pr-1j*Qr);

%% Generator type 2: 4th order model
% delta = Xgen(2,1);
% Eq_tr = Xgen(2,3);
% Ed_tr = Xgen(2,4);
% xd_tr = Pgen(2,8);

% impedance of the transformer
Zt = 1j*0.06;

%% Calculations
% transform the power absorbed by the rectifier station into equivalent
% impedances
% xd_tr---Z_t---Z_eq connected in series
% calculate the equivalent impedances of the rectifier station
% Calculate generator currents for generator 2
% E = (Eq_tr + 1j*Ed_tr)*exp(1j*delta);
U2 = -(Zeq+Zt)*(Is_alphar+1j*Is_betar);
U8 = -Zeq*(Is_alphar+1j*Is_betar);

% renew of Vmmcr1
% Vmmc0 = [Va_alpha; Va_beta; is0; Ps; Qs; id0]
Va_alpha = real(U8)*Vsmax;
Va_beta = imag(U8)*Vsmax;

Is_alpha = Xmmcr1(4,1);
Is_beta = Xmmcr1(5,1);

Ps = 1.5*(Va_alpha*Is_alpha+Va_beta*Is_beta); % real time values
Qs = 1.5*(Va_beta*Is_alpha-Va_alpha*Is_beta); % real time values
Is = sqrt(Is_alpha^2+Is_beta^2);
phi = atan(Qs/Ps);
is = Is*cos(w1*t-phi);

id = -(Vdr-Vdi)/Rdc;

Vmmcr1 = [Va_alpha; Va_beta; Is; Ps; Qs; id];
return;