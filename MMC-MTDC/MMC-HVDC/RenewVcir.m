function [Vciri1, Vcirr1] = RenewVcir(Xciri1, Xcirr1, Pcir0, Vciri0, Vcirr0, Xmmci0, Xmmcr0, Pmmc0, Vouci1, Voucr1)
% Pmmc0=[N M Vd Vsmax Ismax phi Srated w1 L R C Carm Cd Rdc]
M = Pmmc0(1,2);
R = Pmmc0(1,10);
Rdc = Pmmc0(1,14);

% Pcir0 = [Ra alpha2]
Ra = Pcir0(1,1);

% Xmmc0 = [ic0; vcu0; vcl0; Is_alpha0; Is_beta0; Vd0];
% inverter
ici = Xmmci0(1,1);
Vdi = Xmmci0(6,1);
% rectifier
icr = Xmmcr0(1,1);
Vdr = Xmmcr0(6,1);
id = (Vdr-Vdi)/Rdc;   % ×¢Òâµã

% Xcir0 = [xc10; xc20]
% inverter
xc1i = Xciri1(1,1);
% rectifier
xc1r = Xcirr1(1,1);

% Vcir0 = [vc_ref; vc; nu; nl]
% inverter
Vciri1 = Vciri0;
vc_refi = -(R+Ra)*id/M-Ra*xc1i+Ra*ici+Vdi/2;
Vciri1(1,1) = vc_refi;
vci = vc_refi;
Vciri1(2,1) = vci;
% rectifier
Vcirr1 = Vcirr0;
vc_refr = -(R+Ra)*id/M-Ra*xc1r+Ra*icr+Vdr/2;
Vcirr1(1,1) = vc_refr;
vcr = vc_refr;
Vcirr1(2,1) = vcr;

% Vouc0 = [Vs_alpha; Vs_beta; vs_ref; vs]
vs_refi = Vouci1(3,1);
vs_refr = Voucr1(3,1);

% inverter
nui = (vc_refi-vs_refi)/Vdi;
if nui > 1
    nui = 1;
elseif nui < 0
    nui = 0;
end
Vciri1(3,1) = nui;

nli = (vc_refi+vs_refi)/Vdi;
if nli > 1
    nli = 1;
elseif nli < 0
    nli = 0;
end
Vciri1(4,1) = nli;

% rectifier
nur = (vc_refr-vs_refr)/Vdr;
if nur > 1
    nur = 1;
elseif nur < 0
    nur = 0;
end
Vcirr1(3,1) = nur;

nlr = (vc_refr+vs_refr)/Vdr;
if nlr > 1
    nlr = 1;
elseif nlr < 0
    nlr = 0;
end
Vcirr1(4,1) = nlr;

end