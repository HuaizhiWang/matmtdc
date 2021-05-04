function Vmmci1 = RenewVmmci(Xmmci1, Xmmcr1, Pmmc0, U1, t)
% Pmmc0=[N M Vd Vsmax Ismax phi Srated w1 L R C Carm Cd Rdc]
Vsmax = Pmmc0(1,4);
w1 = Pmmc0(1,8);
Rdc = Pmmc0(1,14);

% 2号发电机被MMC-HVDC的inverter station 替代
Uinv = U1(2,1);
Va_alphai = real(Uinv)*Vsmax;
Va_betai = imag(Uinv)*Vsmax;

% Xmmc0 = [ic0; vcu0; vcl0; Is_alpha0; Is_beta0; Vd0]
Is_alphai = Xmmci1(4,1);
Is_betai = Xmmci1(5,1);
Vdi = Xmmci1(6,1);
Vdr = Xmmcr1(6,1);
idi = (Vdr-Vdi)/Rdc;

% inverter
Psi = 1.5*(Va_alphai*Is_alphai+Va_betai*Is_betai); % real time values
Qsi = 1.5*(Va_betai*Is_alphai-Va_alphai*Is_betai); % real time values
Isi = sqrt(Is_alphai^2+Is_betai^2);
phi = atan(Qsi/Psi);
% isi = Isi*cos(w1*t-phi);

% Vmmc0 = [Va_alpha; Va_beta; is0; Ps; Qs; id0]
Vmmci1 = [Va_alphai; Va_betai; Isi; Psi; Qsi; idi];

end