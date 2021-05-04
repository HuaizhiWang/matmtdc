function dFouc0 = OucCon(Xouc0, Pouc0, Vouc0, Xmmc0, Pmmc0, Vmmc0, Addvai, Vdcc0, MMCflow0, Type)
% Pmmc0=[N M Vd Vsmax Ismax phi Srated w1 L R C Carm Cd Rdc];
M = Pmmc0(1,2);
Vd_ref = Pmmc0(1,3);
Vsmax = Pmmc0(1,4);
Ismax = Pmmc0(1,5);
% phi = Pmmc0(1,6);
% Srated = Pmmc0(1,7);
Rdc = Pmmc0(1,14);

% Vmmc0 = [Va_alpha; Va_beta; is0; Ps; Qs; id0]
Va_alpha = Vmmc0(1,1);
Va_beta = Vmmc0(2,1);
Va = sqrt(Va_alpha^2+Va_beta^2);    
% Vai = Addvai(2,1);
Va = Va + 1*Addvai(1,1);

% Xmmc0 = [ic0; vcu0; vcl0; Is_alpha0; Is_beta0; Vd0]
Is_alpha = Xmmc0(4,1);
Is_beta = Xmmc0(5,1);

% Vouc0 = [Vs_alpha; Vs_beta; vs_ref; vs]
Vs_alpha = Vouc0(1,1);
Vs_beta = Vouc0(2,1);

% Xouc0 = [VR_alpha; VR_beta]
VR_alpha = Xouc0(1,1);
VR_beta = Xouc0(2,1);

% Pouc0 = [Kp_oc K0]
Kp_oc = Pouc0(1,1);
K0 = Pouc0(1,2);
Kpos = Ismax/(0.4*Vsmax);

if Type == 1 % inverter Vdc/Vs control
    Qsr = MMCflow0(1,2)*100e6;
    % Vdcc0 = Is_alphar
    Is_alphar = Vdcc0;
    % with LVRT capability: reactive power injection
    if Va >= 0.9*Vsmax
        Is_betar = -2*Qsr/(3*Va);
    else
        Is_betar = (0.9*Vsmax-Va)*Kpos;
    end
elseif Type ==2 % rectifier PQ control
    Ps0 = MMCflow0(2,1)*100e6;
    ic0 = Ps0/M/Vd_ref; % Initial value for ic (A)
    % Vd*id = P
    id0 = ic0*M;
    Psr = Ps0 - Rdc*id0^2; % reference value
    Qsr = MMCflow0(2,2)*100e6;
    Is_alphar = 2*Psr/(3*Va);
    Is_betar = -2*Qsr/(3*Va);
end

dVR_alpha = K0*((Is_alphar-Is_alpha)+(Vs_alpha-(VR_alpha+Va_alpha+Kp_oc*(Is_alphar-Is_alpha)))/Kp_oc);
dVR_beta = K0*((Is_betar-Is_beta)+(Vs_beta-(VR_beta+Va_beta+Kp_oc*(Is_betar-Is_beta)))/Kp_oc);

dFouc0 = [dVR_alpha; dVR_beta];
end