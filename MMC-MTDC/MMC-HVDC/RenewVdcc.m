function Vdcc1 = RenewVdcc(Xdcc1, Pdcc0, Xmmc0, Pmmc0, Vmmc0, Addvai)
% Pmmc0=[N M Vd Vsmax Ismax phi Srated w1 L R C Carm Cd Rdc]
N = Pmmc0(1,1);
M = Pmmc0(1,2);
Vd_ref = Pmmc0(1,3);
C = Pmmc0(1,11);
% Carm = Pmmc0(1,12);
Cd = Pmmc0(1,13);
Cdp = Cd + 2*M*C/N;

% Pdcc0 = [alpha_d alpha_id]
alpha_d = Pdcc0(1,1);

% Xdcc0 = xdc0 for inverter
xdci = Xdcc1;

% Xmmc0 = [ic0; vcu0; vcl0; Is_alpha0; Is_beta0; Vd0];
% inverter
Vdi = Xmmc0(6,1);

% Vmmc0 = [Va_alpha; Va_beta; is0; Ps; Qs; id0]
% inverter
Va_alphai = Vmmc0(1,1);
Va_betai = Vmmc0(2,1);
Va = sqrt(Va_alphai^2+Va_betai^2);
% Vai = Addvai(2,1);
Va = Va + 1*Addvai(1,1);

% Vdcc0 = Is_alphari for inverter
Is_alphari = 2*(xdci+alpha_d*Cdp*(Vdi^2-Vd_ref^2)/2)/(3*Va);
% if Is_alphari > Ismax
%     Is_alphari = Ismax;
% elseif Is_alphari < -Ismax
%     Is_alphari = -Ismax;
% end
Vdcc1 = Is_alphari;

end