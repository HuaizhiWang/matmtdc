function Pouc0 = Loadouc(casefile_dyn, Pmmc0)
% Pmmc0=[N M Vd Vsmax Ismax phi Srated w1 L R C Carm Cd]
L = Pmmc0(1,9);

[gen,exc,gov,mmc,DCVolCon,OutCurCon,CirCurCon,freq,stepsize,stoptime] = feval(casefile_dyn);
% OutCurCon = [alpha_c alpha_1]
alpha_c = OutCurCon(1,1);
alpha_1 = OutCurCon(1,2);

% Refering to P158 (3.108)
% Kp_oc = alpha_c*L/2; K0 = alpha_1*alpha_c*L alpha_1<<alpha_c
Kp_oc = alpha_c*L/2;
K0 = alpha_1*alpha_c*L;

% Pouc0 = [Kp_oc K0]
Pouc0 = [Kp_oc K0];

end