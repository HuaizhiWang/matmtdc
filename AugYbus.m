function [Ly, Uy, Py] = AugYbus(baseMVA, bus, branch, xd_tr, gbus, P, Q, U0)

% [Ly,Uy,Py] = AugYbus(baseMVA, bus, branch, xd_tr, gbus, P, Q, U0)
% 
% Constructs augmented bus admittance matrix Ybus
% 
% INPUTS
% baseMVA = power base
% bus = bus data
% branch = branch data
% xd_tr = d component of transient reactance
% gbus = generator buses
% P = load active power
% Q = load reactive power
% U0 = steady-state bus voltages
% 
% OUTPUTS
% [Ly,Uy,Py] = factorised augmented bus admittance matrix

% MatDyn
% Copyright (C) 2009 Stijn Cole
% Katholieke Universiteit Leuven
% Dept. Electrical Engineering (ESAT), Div. ELECTA
% Kasteelpark Arenberg 10
% 3001 Leuven-Heverlee, Belgium

%%

% Calculate bus admittance matrix
Ybus = makeYbus(baseMVA, bus, branch);

% Calculate equivalent load admittance
yload = (P - 1j.*Q)./(abs(U0).^2);

% Calculate equivalent generator admittance
ygen=zeros(size(Ybus,1),1);
ygen(gbus) = 1./(1j.*xd_tr);
% 2号发电机被替换为MMC-HVDC的inverter station
ygen(2) = 0;

% Add equivalent load and generator admittance to Ybus matrix
for i=1:size(Ybus,1)
    Ybus(i,i) = Ybus(i,i) + ygen(i) + yload(i);
end

Y = Ybus;

% Factorise
[Ly,Uy,Py] = lu(Y,'vector');

return;