function [Xouc0, Vouc0] = OucInit(Xmmc0, Pmmc0, Vmmc0)
% Pmmc0=[N M Vd Vsmax Ismax phi Srated w1 L R C Carm Cd Rdc]
R = Pmmc0(1,10);

% Vs_alpha Vs_beta
% Terminal voltage initialization
% Xmmc0 = [ic0; vcu0; vcl0; Is_alpha0; Is_beta0; Vd0]
Is_alpha0 = Xmmc0(4,1);
Is_beta0 = Xmmc0(5,1);
Vd0 = Xmmc0(6,1);
% Vmmc0 = [Va_alpha; Va_beta; is0; Ps; Qs; id0] Ps and Qs are real time
% values
Va_alpha = Vmmc0(1,1);
Va_beta = Vmmc0(2,1);
Vs_alpha = Va_alpha + (R/2)*Is_alpha0;
Vs_beta = Va_beta + (R/2)*Is_beta0;

if Vs_alpha > Vd0
    Vs_alpha = Vd0;
elseif Vs_alpha < -Vd0
    Vs_alpha = -Vd0;
end

if Vs_beta > Vd0
    Vs_beta = Vd0;
elseif Vs_beta < -Vd0
    Vs_beta = -Vd0;
end

VR_alpha = Vs_alpha - Va_alpha;
VR_beta = Vs_beta - Va_beta;
% Xouc0 = [VR_alpha; VR_beta]
Xouc0 = [VR_alpha; VR_beta];

Vs = sqrt(Vs_alpha^2 + Vs_beta^2);
vs_ref = Vs;
vs = vs_ref;
% Vouc0 = [Vs_alpha; Vs_beta; vs_ref; vs]
Vouc0 = [Vs_alpha; Vs_beta; vs_ref; vs];

end