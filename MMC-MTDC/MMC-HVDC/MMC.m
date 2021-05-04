function dFmmcr0 = MMC(Xmmc0, Pmmc0, Vmmc0, Vouc0, Vcir0, t)

%% Parameter input
% Xmmc0 = [ic0; vcu0; vcl0; Is_alpha0; Is_beta0; Vd0];
ic = Xmmc0(1,1);
vcu = Xmmc0(2,1);
vcl = Xmmc0(3,1);
Is_alpha = Xmmc0(4,1);
Is_beta = Xmmc0(5,1);
Vd = Xmmc0(6,1);

% Pmmc0=[N M Vd Vsmax Ismax phi Srated w1 L R C Carm Cd]
N = Pmmc0(1,1);
M = Pmmc0(1,2);
% Vd = Pmmc0(1,3);
% Vsmax = Pmmc0(1,4);
% Ismax = Pmmc0(1,5);
% phi = Pmmc0(1,6);
% Srated = Pmmc0(1,7);
w1 = Pmmc0(1,8);
L = Pmmc0(1,9);
R = Pmmc0(1,10);
C = Pmmc0(1,11);
Carm = Pmmc0(1,12);
Cd = Pmmc0(1,13);
Cdp = Cd+2*M*C/N;

% Vmmc0 = [Va_alpha; Va_beta; is0; Ps; Qs; id0]
Va_alpha = Vmmc0(1,1);
Va_beta = Vmmc0(2,1);
is = Vmmc0(3,1);
Ps = Vmmc0(4,1);
Qs = Vmmc0(5,1);
id0 = Vmmc0(6,1);
% is = is * atan(Qs/Ps);
phi = atan(Qs/Ps);
is = is * cos(w1*t-phi);

% Vouc0 = [Vs_alpha; Vs_beta; vs_ref; vs]
Vs_alpha = Vouc0(1,1);
Vs_beta = Vouc0(2,1);
% vs_ref = Vouc0(3,1);
% vs = Vouc0(4,1);

% Vcir0 = [vc_ref; vc; nu; nl]
% vc_ref = Vcir0(1,1);
nu = Vcir0(3,1);
nl = Vcir0(4,1);

%% State equations
dic = (Vd-nu*vcu-nl*vcl)/(2*L)-R*ic/L;
dvcu = nu*(ic+0.5*is)/Carm;
dvcl = nl*(ic-0.5*is)/Carm;
dIs_alpha = (Vs_alpha-Va_alpha-R*Is_alpha/2)*2/L;
dIs_beta = (Vs_beta-Va_beta-R*Is_beta/2)*2/L;
if Vd < 1e-6 && Vd >= 0
    dVd = (id0-Ps/1e-6)/Cdp;
elseif Vd > -1e-6 && Vd < 0
    dVd = (id0+Ps/1e-6)/Cdp;
else
    dVd = (id0-Ps/Vd)/Cdp;
end

dFmmcr0 = [dic; dvcu; dvcl; dIs_alpha; dIs_beta; dVd];

end