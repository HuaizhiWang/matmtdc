function dCir0 = CirCon(Xcir0, Pcir0, Vcir0, Xmmc0, Pmmc0)
% Xcir0 = [xc10; xc20]
xc1 = Xcir0(1,1);
xc2 = Xcir0(2,1);

% Pcir0 = [Ra alpha2]
% Ra = Pcir0(1,1);
alpha2 = Pcir0(1,2);

% Vcir0 = [vc_ref; vc; nu; nl]
% vc_ref = Vcir0(1,1);
% vc = Vcir0(2,1);
nu = Vcir0(3,1);
nl = Vcir0(4,1);

% Xmmc0 = [ic0; vcu0; vcl0; Is_alpha0; Is_beta0; Vd0];
ic = Xmmc0(1,1);
vcu = Xmmc0(2,1);
vcl = Xmmc0(3,1);
% Is_alpha = Xmmc0(4,1);
% Is_beta = Xmmc0(5,1);
Vd = Xmmc0(6,1);

% Pmmc0=[N M Vd Vsmax Ismax phi Srated w1 L R C Carm Cd]
% N = Pmmc0(1,1);
% M = Pmmc0(1,2);
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

dxc1 = xc2;
dxc2 = -2*alpha2*((Vd-nu*vcu-nl*vcl)/(2*L)-R*ic/L)-4*w1^2*xc1;

dCir0 = [dxc1; dxc2];

end