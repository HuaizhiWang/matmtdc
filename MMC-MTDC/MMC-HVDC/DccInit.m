function [Xdcc0, Vdcc0] = DccInit(Pdcc0, Xmmc0, Pmmc0, Vmmc0)
% Xmmc0 = [ic0; vcu0; vcl0; Is_alpha0; Is_beta0, Vd0]
Is_alpha0 = Xmmc0(4,1);
Vd0 = Xmmc0(6,1);
% Vmmc0 = [Va_alpha; Va_beta; is0; Ps; Qs; id0]
Va_alpha = Vmmc0(1,1);
Va_beta = Vmmc0(2,1);
Va = sqrt(Va_alpha^2+Va_beta^2);

xdc0 = 3*Va*Is_alpha0/2;
% Xdcc0 = xdc0
Xdcc0 = xdc0;

% Pmmc0=[N M Vd Vsmax Ismax phi Srated w1 L R C Carm Cd]
N = Pmmc0(1,1);
M = Pmmc0(1,2);
Vd_ref = Pmmc0(1,3);
Ismax = Pmmc0(1,5);
C = Pmmc0(1,11);
Cd = Pmmc0(1,13);
Cdp = Cd+2*M*C/N;
% Pdcc0 = [alpha_d alpha_id]
alpha_d = Pdcc0(1,1);
% alpha_id = Pdcc0(1,2);

Is_alpha1r = 2*(xdc0+alpha_d*Cdp*(Vd0^2-Vd_ref^2)/2)/(3*Va);
if Is_alpha1r > Ismax
    Is_alphar = Ismax;
elseif Is_alpha1r < -Ismax
    Is_alphar = -Ismax;
else
    Is_alphar = Is_alpha1r;
end
% Vdcc0 = Is_alphar
Vdcc0 = Is_alphar;
end