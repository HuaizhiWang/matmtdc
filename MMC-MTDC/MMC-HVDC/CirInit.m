function [Xcir0, Vcir0] = CirInit(Xmmc0, Pmmc0, Vmmc0, Vouc0, Pcir0, Type)

% Xmmc0 = [ic0; vcu0; vcl0; Is_alpha0; Is_beta0; Vd0]
ic = Xmmc0(1,1);
% vcu = Xmmc0(2,1);
% vcl = Xmmc0(3,1);
% Is_alpha = Xmmc0(4,1);
% Is_beta = Xmmc0(5,1);
Vd0 = Xmmc0(6,1);

% Pmmc0=[N M Vd Vsmax Ismax phi Srated w1 L R C Carm Cd]
% N = Pmmc0(1,1);
M = Pmmc0(1,2);
% Vd = Pmmc0(1,3);
% Vsmax = Pmmc0(1,4);
% Ismax = Pmmc0(1,5);
% phi = Pmmc0(1,6);
% Srated = Pmmc0(1,7);
w1 = Pmmc0(1,8);
L = Pmmc0(1,9);
R = Pmmc0(1,10);
% C = Pmmc0(1,11);
% Carm = Pmmc0(1,12);
% Cd = Pmmc0(1,13);

% Vmmc0 = [Va_alpha; Va_beta; is0; Ps; Qs; id0]
% Va_alpha = Vmmc0(1,1);
% Va_beta = Vmmc0(2,1);
% is = Vmmc0(3,1);
if Type == 1
    id0 = Vmmc0(6,1);
elseif Type == 2
    id0 = -Vmmc0(6,1);
end

% Vouc0 = [Vs_alpha; Vs_beta; vs_ref; vs]
% Vs_alpha = Vouc0(1,1);
% Vs_beta = Vouc0(2,1);
vs_ref = Vouc0(3,1);
% vs = Vmmc0(4,1);

% Pcir0 = [Ra alpha2];
Ra = Pcir0(1,1);
alpha2 = Pcir0(1,2);

xc20 = 0;
xc10 = (-2*alpha2*(R+Ra)*id0/(L*M)+2*alpha2*(R+Ra)*ic/L)/(2*alpha2*Ra/L+4*w1^2);

vc_ref = -(R+Ra)*id0/M-Ra*xc10+Ra*ic+Vd0/2;
vc = vc_ref;

nu = (vc_ref-vs_ref)/Vd0;
if nu > 1
    nu = 1;
elseif nu < 0
    nu=0;
end
nl = (vc_ref+vs_ref)/Vd0;
if nl > 1
    nl = 1;
elseif nl < 0
    nl = 0;
end

% Xcir0 = [xc10; xc20]
Xcir0 = [xc10; xc20];
% Vcir0 = [vc_ref; vc; nu; nl]
Vcir0 = [vc_ref; vc; nu; nl];

end