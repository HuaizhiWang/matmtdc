function [Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, Xmmci0, Xmmcr0, Pmmc0, Vmmci0, Vmmcr0, Xouci0, Xoucr0, Pouc0, Vouci0, Voucr0, Xciri0, Xcirr0, Pcir0, Vciri0, Vcirr0, Xdcci0, Pdcc0, Vdcci0, U0, t, stepsize] = ModifiedEuler(t, Xgen0, Pgen, Vgen0, Xexc0, Pexc, Vexc0, Xgov0, Pgov, Vgov0, Xmmci0, Xmmcr0, Pmmc0, Vmmci0, Vmmcr0, Xouci0, Xoucr0, Pouc0, Vouci0, Voucr0, Xciri0, Xcirr0, Pcir0, Vciri0, Vcirr0, Xdcci0, Pdcc0, Vdcci0, MMCflow0, Ly, Uy, Py, gbus, genmodel, excmodel, govmodel, stepsize, U0,W_sumi,W_deltai,W_sumr,W_deltar)

% Modified Euler ODE solver
 
% MatDyn
% Copyright (C) 2009 Stijn Cole
% Katholieke Universiteit Leuven
% Dept. Electrical Engineering (ESAT), Div. ELECTA
% Kasteelpark Arenberg 10
% 3001 Leuven-Heverlee, Belgium

% Xmmc0 = [ic0; vcu0; vcl0; Is_alpha0; Is_beta0; Vd0]
% Pmmc0=[N M Vd Vsmax Ismax phi Srated w1 L R C Carm Cd Rdc]
% Vmmc0 = [Va_alpha; Va_beta; is0; Ps; Qs; id0]
% Xouc0 = [VR_alpha; VR_beta]
% Pouc0 = [Kp_oc K0]
% Vouc0 = [Vs_alpha; Vs_beta; vs_ref; vs]
% Xcir0 = [xc10; xc20]
% Pcir0 = [Ra alpha2]
% Vcir0 = [vc_ref; vc; nu; nl]
% Xdcc0 = xdc0
% Pdcc0 = [alpha_d alpha_id]
% Vdcc0 = Is_alphar
% Pgrid0=[GridPhaseInit SCR Lg Rg Rdc]
Type1 = 1; % inverter
Type2 = 2; % rectifier

%% First Euler step
%% Synchronous generators    
% EXCITERS
dFexc0 = Exciter(Xexc0, Pexc, Vexc0, excmodel);
Xexc1 = Xexc0 + stepsize.*dFexc0;

% GOVERNORS
dFgov0 = Governor(Xgov0, Pgov, Vgov0, govmodel);
Xgov1 = Xgov0 + stepsize.*dFgov0;
        
% GENERATORS
dFgen0 = Generator(Xgen0, Xexc1, Xgov1, Pgen, Vgen0, genmodel);
Xgen1 = Xgen0 + stepsize.*dFgen0;

%% MMC-HVDC
% DC voltage controller for inverter
dDcci0 = DvotCon(Xdcci0, Pdcc0, Vdcci0, Vmmci0);
Xdcci1 = Xdcci0 + stepsize.*dDcci0;
Vdcci1 = RenewVdcc(Xdcci1, Pdcc0, Xmmci0, Pmmc0, Vmmci0);

% Output current controller for inverter
dFouci0 = OucCon(Xouci0, Pouc0, Vouci0, Xmmci0, Pmmc0, Vmmci0, Vdcci1, MMCflow0, Type1);
Xouci1 = Xouci0 + stepsize.*dFouci0;
% Output current controller for rectifier
Vdccr0 = 0;
dFoucr0 = OucCon(Xoucr0, Pouc0, Voucr0, Xmmcr0, Pmmc0, Vmmcr0, Vdccr0, MMCflow0, Type2);
Xoucr1 = Xoucr0 + stepsize.*dFoucr0;
[Vouci1, Voucr1] = RenewVouc(Xouci1, Xoucr1, Pouc0, Vouci0, Voucr0, Xmmci0, Xmmcr0, Pmmc0, Vmmci0, Vmmcr0, Vdcci1, MMCflow0, t);

% Circulating current controller for inverter
 % Sum-Capacitor-Voltage Ripples
dCiri0 = CirCon(Xciri0, Pcir0, Vciri0, Xmmci0, Pmmc0,MMCflow0,W_sumi,W_deltai,Type1,t);
Xciri1 = Xciri0 + stepsize.*dCiri0;
% Circulating current controller for rectifier
dCirr0 = CirCon(Xcirr0, Pcir0, Vcirr0, Xmmcr0, Pmmc0,MMCflow0,W_sumr,W_deltar,Type2,t);
Xcirr1 = Xcirr0 + stepsize.*dCirr0;
[Vciri1, Vcirr1] = RenewVcir(Xciri1, Xcirr1, Pcir0, Vciri0, Vcirr0, Xmmci0, Xmmcr0, Pmmc0, Vouci1, Voucr1);

% MMC inverter
dFmmci0 = MMC(Xmmci0, Pmmc0, Vmmci0, Vouci1, Vciri1);
Xmmci1 = Xmmci0 + stepsize.*dFmmci0;
% MMC rectifier
dFmmcr0 = MMC(Xmmcr0, Pmmc0, Vmmcr0, Voucr1, Vcirr1);
Xmmcr1 = Xmmcr0 + stepsize.*dFmmcr0;
       
%% Calculate system voltages
U1 = SolveNetwork(Xgen1, Pgen, Ly, Uy, Py, gbus, genmodel, Xmmci1, Pmmc0);

% re-new Vmmc for the inverter of the MMC-HVDC
% 2号节点成为inverter station的输出母线
Vmmci1 = RenewVmmci(Xmmci1, Xmmcr1, Pmmc0, U1, t);

% Calculate the real terminal voltage of generator 2
% 8号节点成为rectifier station的输出母线
[Vmmcr1, Ug2] = RenewVmmcr(Xmmci1, Xmmcr1, Pmmc0, Vmmcr0, t);

% Calculate machine currents and power
U1(2,1) = Ug2;
[Id1,Iq1,Pe1] = MachineCurrents(Xgen1, Pgen, U1(gbus), genmodel);    
   
% Update variables that have changed
Vexc1 = abs(U1(gbus));
Vgen1 = [Id1,Iq1,Pe1];   
Vgov1 = [Xgen1(:,2)];


%% Second Euler step    
%% Synchronous generators    
% EXCITERS
dFexc1 = Exciter(Xexc1, Pexc, Vexc1, excmodel);
Xexc2 = Xexc0 + stepsize/2 .* (dFexc0 + dFexc1);
    
% GOVERNORS
dFgov1 = Governor(Xgov1, Pgov, Vgov1, govmodel);
Xgov2 = Xgov0 + stepsize/2 .* (dFgov0 + dFgov1);
              
% GENERATORS
dFgen1 = Generator(Xgen1, Xexc2, Xgov2, Pgen, Vgen1, genmodel);
Xgen2 = Xgen0 + stepsize/2 .* (dFgen0 + dFgen1);

%% MMC-HVDC
% DC voltage controller for inverter
dDcci1 = DvotCon(Xdcci1, Pdcc0, Vdcci1, Vmmci1);
Xdcci2 = Xdcci0 + stepsize/2.*(dDcci0+dDcci1);
Vdcci2 = RenewVdcc(Xdcci2, Pdcc0, Xmmci1, Pmmc0, Vmmci1);

% Output current controller for inverter
dFouci1 = OucCon(Xouci1, Pouc0, Vouci1, Xmmci1, Pmmc0, Vmmci1, Vdcci2, MMCflow0, Type1);
Xouci2 = Xouci0 + stepsize/2.*(dFouci0+dFouci1);
% Output current controller for rectifier
Vdccr2 = 0;
dFoucr1 = OucCon(Xoucr1, Pouc0, Voucr1, Xmmcr1, Pmmc0, Vmmcr1, Vdccr2, MMCflow0, Type2);
Xoucr2 = Xoucr0 + stepsize/2.*(dFoucr0+dFoucr1);
[Vouci2, Voucr2] = RenewVouc(Xouci2, Xoucr2, Pouc0, Vouci1, Voucr1, Xmmci1, Xmmcr1, Pmmc0, Vmmci1, Vmmcr1, Vdcci2, MMCflow0, t);

% Circulating current controller for inverter
dCiri1 = CirCon(Xciri1, Pcir0, Vciri1, Xmmci1, Pmmc0,MMCflow0,W_sumi,W_deltai,Type1,t);
Xciri2 = Xciri0 + stepsize/2.*(dCiri0+dCiri1);
% Circulating current controller for rectifier
dCirr1 = CirCon(Xcirr1, Pcir0, Vcirr1, Xmmcr1, Pmmc0,MMCflow0,W_sumr,W_deltar,Type2,t);
Xcirr2 = Xcirr0 + stepsize/2.*(dCirr0+dCirr1);
[Vciri2, Vcirr2] = RenewVcir(Xciri2, Xcirr2, Pcir0, Vciri1, Vcirr1, Xmmci1, Xmmcr1, Pmmc0, Vouci2, Voucr2);

% MMC inverter
dFmmci1 = MMC(Xmmci1, Pmmc0, Vmmci1, Vouci2, Vciri2);
Xmmci2 = Xmmci0 + stepsize/2.*(dFmmci0+dFmmci1);
% MMC rectifier
dFmmcr1 = MMC(Xmmcr1, Pmmc0, Vmmcr1, Voucr2, Vcirr2);
Xmmcr2 = Xmmcr0 + stepsize/2.*(dFmmcr0+dFmmcr1);
 
% Calculate system voltages
U2 = SolveNetwork(Xgen2, Pgen, Ly, Uy, Py, gbus, genmodel, Xmmci2, Pmmc0);

% re-new Vmmc for the inverter of the MMC-HVDC
% 2号节点成为inverter station的输出母线
Vmmci2 = RenewVmmci(Xmmci2, Xmmcr2, Pmmc0, U2, t);

% Calculate the real terminal voltage of generator 2
% 8号节点成为rectifier station的输出母线
[Vmmcr2, Ug2] = RenewVmmcr(Xmmci2, Xmmcr2, Pmmc0, Vmmcr1, t);

% Calculate machine currents and power
U2(2,1) = Ug2;
[Id2,Iq2,Pe2] = MachineCurrents(Xgen1, Pgen, U2(gbus), genmodel); 

% Update variables that have changed
Vgen2 = [Id2,Iq2,Pe2];
Vexc2 = abs(U2(gbus));
Vgov2 = [Xgen2(:,2)];

%% Update

U0 = U2;
    
Vgen0 = Vgen2;
Vgov0 = Vgov2;
Vexc0 = Vexc2;
Vdcci0 = Vdcci2;
Vouci0 = Vouci2;
Voucr0 = Voucr2;
Vciri0 = Vciri2;
Vcirr0 = Vcirr2;
Vmmci0 = Vmmci2;
Vmmcr0 = Vmmcr2;
    
Xgen0 = Xgen2; 
Xexc0 = Xexc2;
Xgov0 = Xgov2;
Xdcci0 = Xdcci2;
Xouci0 = Xouci2;
Xoucr0 = Xoucr2;
Xciri0 = Xciri2;
Xcirr0 = Xcirr2;
Xmmci0 = Xmmci2;
Xmmcr0 = Xmmcr2;
    
Pgen0 = Pgen; 
Pexc0 = Pexc;
Pgov0 = Pgov;         
 
return;