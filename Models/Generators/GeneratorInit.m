function [Efd0,Xgen0] = GeneratorInit(Pgen, U0, gen, baseMVA, gentype)

% [Efd0,Xgen0] = GeneratorInit(Pgen, U0, gen, baseMVA, gentype)
%
% calculate initial conditions generators
 
% MatDyn
% Copyright (C) 2009 Stijn Cole
% Katholieke Universiteit Leuven
% Dept. Electrical Engineering (ESAT), Div. ELECTA
% Kasteelpark Arenberg 10
% 3001 Leuven-Heverlee, Belgium

%% Init
global freq;

[ngen,c] = size(Pgen);
Xgen0 = zeros(ngen,1);
Efd0 = zeros(ngen,1);
d=[1:length(gentype)]';

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% Define types
type1 = d(gentype==1);
type2 = d(gentype==2);

%% Generator type 1: classical model
x_tr = Pgen(type1,7);

omega0=ones(length(type1),1).*2.*pi.*freq;

% Initial machine armature currents
Ia0 = (gen(type1,PG) - j.*gen(type1,QG))./( conj(U0(type1,1)))./baseMVA;

% Initial Steady-state internal EMF

Eq_tr0 = U0(type1,1) + j.*x_tr.*Ia0;
delta0 = angle(Eq_tr0);
Eq_tr0 = abs(Eq_tr0);

Xgen0(type1,1:3) = [delta0 omega0 Eq_tr0];

%% Generator type 2: 
xd = Pgen(type2,6);
xq = Pgen(type2,7);
xd_tr = Pgen(type2,8);
xq_tr = Pgen(type2,9);

omega0=ones(length(type2),1).*2.*pi.*freq;

% Initial machine armature currents
Ia0=(gen(type2,PG) - j.*gen(type2,QG))./(conj(U0(type2,1)))./baseMVA;
phi0=angle(Ia0);

% Initial Steady-state internal EMF
Eq0 = U0(type2,1) + j.*xq.*Ia0;
delta0 = angle(Eq0);

% Machine currents in dq frame
Id0 = -abs(Ia0).*sin(delta0 - phi0);
Iq0 = abs(Ia0).*cos(delta0 - phi0);

% Field voltage
Efd0(type2) = abs(Eq0) - (xd - xq).*Id0;

% Initial Transient internal EMF
Eq_tr0 = Efd0(type2,1) + (xd - xd_tr).*Id0;
Ed_tr0 = -(xq - xq_tr).*Iq0;

Xgen0(type2,1:4) = [delta0, omega0, Eq_tr0, Ed_tr0];

%% Generator type 3: 

%% Generator type 4: 

return;