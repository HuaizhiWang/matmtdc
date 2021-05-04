function [Vouci1, Voucr1] = RenewVouc(Xouci1, Xoucr1, Pouc0, Vouci0, Voucr0, Xmmci0, Xmmcr0, Pmmc0, Vmmci0, Addvai, Vmmcr0, Vdcci1, MMCflow0, t)

% Pmmc0=[N M Vd Vsmax Ismax phi Srated w1 L R C Carm Cd Rdc]
Vd_ref = Pmmc0(1,3);
Vsmax = Pmmc0(1,4);
Ismax = Pmmc0(1,5);
w1 = Pmmc0(1,8);
Rdc = Pmmc0(1,14);
Kpos = Ismax/(0.4*Vsmax);

% Vmmc0 = [Va_alpha; Va_beta; is0; Ps; Qs; id0]
% inverter
Va_alphai = Vmmci0(1,1);
Va_betai = Vmmci0(2,1);
Vai = sqrt(Va_alphai^2+Va_betai^2);
% Va = Addvai(2,1);
Vai = Vai + 1*Addvai(1,1);
% rectifier
Va_alphar = Vmmcr0(1,1);
Va_betar = Vmmcr0(2,1);
Var = sqrt(Va_alphar^2+Va_betar^2);

% MMCflow0 = [Prinv Qrinv; Prrec Qrrec] in p.u.
Psri = MMCflow0(1,1)*100e6;
Qsri = MMCflow0(1,2)*100e6;
Psrr = MMCflow0(2,1)*100e6;
id0 = Psrr/Vd_ref;
Psrr = Psrr - id0^2*Rdc;
Qsrr = MMCflow0(2,2)*100e6;

% Xouc0 = [VR_alpha; VR_beta]
% inverter
VR_alphai = Xouci1(1,1);
VR_betai = Xouci1(2,1);
% rectifier
VR_alphar = Xoucr1(1,1);
VR_betar = Xoucr1(2,1);

% Pouc0 = [Kp_oc K0]
Kp_oc = Pouc0(1,1);

% Xmmc0 = [ic0; vcu0; vcl0; Is_alpha0; Is_beta0; Vd0];
% inverter
Is_alphai = Xmmci0(4,1);
Is_betai = Xmmci0(5,1);
Vdi = Xmmci0(6,1);
% rectifier
Is_alphar = Xmmcr0(4,1);
Is_betar = Xmmcr0(5,1);
Vdr = Xmmcr0(6,1);

% inverter
Is_alphari = Vdcci1;
% rectifier
Is_alpharr = (Va_alphar*Psrr+Va_betar*Qsrr)/(1.5*(Va_alphar^2+Va_betar^2));

% with LVRT capability: reactive power injection
% inverter
if Vai >= 0.9*Vsmax
    Is_betari = (Va_betai*Psri-Va_alphai*Qsri)/(1.5*(Va_alphai^2+Va_betai^2));
else
    Is_betari = (0.9*Vsmax-Vai)*Kpos;
end
% if Is_betari > Ismax
%     Is_betari = Ismax;
% elseif Is_betari < -Ismax
%     Is_betari = -Ismax;
% end
Vs_alphai = VR_alphai+Va_alphai+Kp_oc*(Is_alphari-Is_alphai);
if Vs_alphai > Vdi
    Vs_alphai = Vdi;
elseif Vs_alphai < -Vdi
    Vs_alphai = -Vdi;
end
Vs_betai = VR_betai+Va_betai+Kp_oc*(Is_betari-Is_betai);
if Vs_betai > Vdi
    Vs_betai = Vdi;
elseif Vs_betai < -Vdi
    Vs_betai = -Vdi;
end
Vsi = sqrt(Vs_alphai^2+Vs_betai^2);
vs_refi = Vsi*cos(w1*t);
vsi=vs_refi;
% inverter
Vouci1 = [Vs_alphai; Vs_betai; vs_refi; vsi];



% rectifier
if Var >= 0.9*Vsmax
    Is_betarr =(Va_betar*Psrr-Va_alphar*Qsrr)/(1.5*(Va_alphar^2+Va_betar^2));
else
    Is_betarr = (0.9*Vsmax-Var)*Kpos;
end
% if Is_betarr > Ismax
%     Is_betarr = Ismax;
% elseif Is_betarr < -Ismax
%     Is_betarr = -Ismax;
% end
Vs_alphar = VR_alphar+Va_alphar+Kp_oc*(Is_alpharr-Is_alphar);
if Vs_alphar > Vdr
    Vs_alphar = Vdr;
elseif Vs_alphar < -Vdr
    Vs_alphar = -Vdr;
end
Vs_betar = VR_betar+Va_betar+Kp_oc*(Is_betarr-Is_betar);
if Vs_betar > Vdr
    Vs_betar = Vdr;
elseif Vs_betar < -Vdr
    Vs_betar = -Vdr;
end
Vsr = sqrt(Vs_alphar^2+Vs_betar^2);
vs_refr = Vsr*cos(w1*t);
vsr=vs_refr;
Voucr1 = [Vs_alphar; Vs_betar; vs_refr; vsr];

end