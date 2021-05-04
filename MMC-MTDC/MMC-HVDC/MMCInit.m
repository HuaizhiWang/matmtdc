function [Xmmc0, Vmmc0] = MMCInit(Pmmc0, MMCflow, U0, Type)
% Pmmc0=[N M Vd_ref Vsmax Ismax phi Srated w1 L R C Carm Cd Rdc]
% N = Pmmc0(1,1);
M = Pmmc0(1,2);
Vd_ref = Pmmc0(1,3);
Vsmax = Pmmc0(1,4);
% Ismax = Pmmc0(1,5);
% phi = Pmmc0(1,6);
% Srated = Pmmc0(1,7);
% w1 = Pmmc0(1,8);
% L = Pmmc0(1,9);
% R = Pmmc0(1,10);
% C = Pmmc0(1,11);
% Carm = Pmmc0(1,12);
% Cd = Pmmc0(1,13);
Rdc = Pmmc0(1,14);

% Power flow initialization
if Type == 1 % inverter, and the output power Ps and Qs are positive
    Ps0 = MMCflow(1,1)*100e6; % reference value
    Qs0 = MMCflow(1,2)*100e6; % reference value
    Vd0 = Vd_ref;
    ic0 = Ps0/M/Vd0; % Initial value for ic (A)
    id0 = ic0*M;
elseif Type == 2 % rectifier, and the output power Ps and Qs are negative
    Ps0 = MMCflow(2,1)*100e6;
    ic0 = Ps0/M/Vd_ref; % Initial value for ic (A)
    % Vd*id = P
    id0 = ic0*M;
    Vd0 = Vd_ref - id0*Rdc;
    Ps0 = Ps0 - Rdc*id0^2; % reference value
    Qs0 = MMCflow(2,2)*100e6; % reference value
    ic0 = -Ps0/M/Vd0; % Initial value for ic (A)
end

% the second and third state
vcu0 = Vd0;
vcl0 = Vd0;

% Terminal voltage initialization
% MMC-HVDC is connected between bus 2 and bus 8
if Type == 1
    Va_alpha = real(U0(2,1))*Vsmax;
    Va_beta = imag(U0(2,1))*Vsmax;
elseif Type == 2
    % 对于rectifier的情况，需要重新计算潮流，使得节点8连接rectifier后收到Ps的有功 不发出无功
    % 发电机节点2的电压任然采用U0(2,1), ie
    u2x = real(U0(2,1));
    u2y = imag(U0(2,1));
    u8y = (-(u2x^2*u2y+u2y^3+2*(Ps0/100e6)*0.06*u2x)-sqrt((u2x^2*u2y+u2y^3+2*(Ps0/100e6)*0.06*u2x)^2-4*(u2x*u2y+(Ps0/100e6)*0.06)*(0.06*Ps0/100e6)))/(-2);
    u8x = (-0.06*Ps0/100e6+u2x*u8y)/u2y;
    Va_alpha = u8x*Vsmax;
    Va_beta = u8y*Vsmax;
    % 测试节点2和8输出的有功无功功率
%     Pg2 = (u2y*u8x-u2x*u8y)/0.06;
%     Qg2 = (1-u2x*u8x-u2y*u8y)/0.06;
%     Qs0 = (u2x*u8x+u2y*u8y-u8x^2-u8y^2)/0.06;
end

% the fourth and fifth state
% Ps0 + j Qs0 = 1.5*(Va_alpha + j*Va_beta)*(is_alpha - j*is_beta)
Is_alpha0 = (Va_alpha*Ps0+Va_beta*Qs0)/(1.5*(Va_alpha^2+Va_beta^2));
Is_beta0 = (Va_beta*Ps0-Va_alpha*Qs0)/(1.5*(Va_alpha^2+Va_beta^2));
Is0 = sqrt(Is_alpha0^2+Is_beta0^2);

% vs_ref = Vs*cos(w1*t) --> t=0 -->vs_ref = Vs
% vc_ref = Vd/2
% is = Is*cos(w1*t-phi) --> t=0 -->is = Is*cos(phi)
% phi = atan(Qs0/Ps0)
phi = atan(Qs0/Ps0);
is0 = Is0;

Ps = Ps0; % real time value
Qs = Qs0; % real time value

% Xmmc0 = [ic0; vcu0; vcl0; Is_alpha0; Is_beta0; Vd0]
Xmmc0 = [ic0; vcu0; vcl0; Is_alpha0; Is_beta0; Vd0];

% Vmmc0 = [Va_alpha; Va_beta; is0; Ps; Qs; id0]
Vmmc0 = [Va_alpha; Va_beta; is0; Ps; Qs; id0];

end