function Pmmc0 = Loadmmc(casefile_dyn)

[gen,exc,gov,mmc,DCVolCon,OutCurCon,CirCurCon,freq,stepsize,stoptime] = feval(casefile_dyn);

% mmc=[N M Vd Vsmax Ismax phi];
N = mmc(1,1); % number of submodels in each arm
M = mmc(1,2); % number of phases
Vd = mmc(1,3); % direct voltage
Vsmax = mmc(1,4); % peak line-to-neutral voltage
Ismax = mmc(1,5); % peak output current
% phi = mmc(1,6); % input power angle in [-pi/2, pi/2] rad

% parameter calculation for mmc
Srated = M*Vsmax*Ismax/2; % VA rated power
w1 = 2*pi*freq;
L = 0.08*Vsmax/(w1*Ismax/2); % arm inductance H = 8% of base impedance
R = 0.1*L*w1; % arm resistance 10% of X
Wrated = 30e-3; % J/VA for m=1 Q=0
C = N*2*Wrated*Srated/6/Vd^2; %F
Carm = C/N; % F
Cd = 100e-6; % pole to pole capacitance P145
% Vsm = Vd/N; % average voltage of SM [V]
% Vs = Vsmax*sqrt(3)/sqrt(2); % grid line-line rms voltage at PCC

% resistance of the dc transmission line
Rdc = 1.39; % 100km dc cable: 13.9 mili-Ohm/km

% Pmmc0=[N M Vd Vsmax Ismax phi Srated w1 L R C Carm Cd Rdc];
Pmmc0 = [mmc Srated w1 L R C Carm Cd Rdc];
end