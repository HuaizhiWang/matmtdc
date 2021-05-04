function [Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, U0, t, stepsize] = RungeKutta(t0, Xgen0, Pgen, Vgen0, Xexc0, Pexc, Vexc0, Xgov0, Pgov, Vgov0, Ly, Uy, Py, gbus, genmodel, excmodel, govmodel, stepsize)

% [Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, U0, t,
% stepsize] = RungeKutta(t0, Xgen0, Pgen, Vgen0, Xexc0, Pexc, Vexc0, Xgov0
% Pgov, Vgov0, Y, gbus, faulted_bus, theta0, genmodel, excmodel, govmodel
% stepsize)
%  
% Standard 4th order Runge-Kutta ODE solver
 
% MatDyn
% Copyright (C) 2009 Stijn Cole
% Katholieke Universiteit Leuven
% Dept. Electrical Engineering (ESAT), Div. ELECTA
% Kasteelpark Arenberg 10
% 3001 Leuven-Heverlee, Belgium

%% Runge-Kutta coefficients

% c = [0 1/2 1/2 1]; not used
a = [0 0 0 0
     1/2 0 0 0
     0 1/2 0 0 
     0 0 1 0 ];
b = [1/6 2/6 2/6 1/6];

%% K1

% EXCITERS
Kexc1 = Exciter(Xexc0, Pexc, Vexc0, excmodel);
Xexc1 = Xexc0 + stepsize.*a(2,1).*Kexc1;

% GOVERNORS
Kgov1 = Governor(Xgov0, Pgov, Vgov0, govmodel);
Xgov1 = Xgov0 + stepsize.*a(2,1).*Kgov1;
        
% GENERATORS
Kgen1 = Generator(Xgen0, Xexc1, Xgov1, Pgen, Vgen0, genmodel);
Xgen1 = Xgen0 + stepsize.*a(2,1).*Kgen1;
     
% Calculate system voltages
U1 = SolveNetwork(Xgen1, Pgen, Ly, Uy, Py, gbus, genmodel);

% Calculate machine currents and power
[Id1,Iq1,Pe1] = MachineCurrents(Xgen1, Pgen, U1(gbus), genmodel);    
   
% Update variables that have changed
Vexc1 = abs(U1(gbus));
Vgen1 = [Id1,Iq1,Pe1];
Vgov1 = [Xgen1(:,2)];

%% K2

% EXCITERS
Kexc2 = Exciter(Xexc1, Pexc, Vexc1, excmodel);
Xexc2 = Xexc0 + stepsize.* (a(3,1).*Kexc1 + a(3,2).*Kexc2 );

% GOVERNORS
Kgov2 = Governor(Xgov1, Pgov, Vgov1, govmodel);
Xgov2 = Xgov0 + stepsize.* (a(3,1).*Kgov1 + a(3,2).*Kgov2 );
        
% GENERATORS
Kgen2 = Generator(Xgen1, Xexc2, Xgov2, Pgen, Vgen1, genmodel);
Xgen2 = Xgen0 + stepsize.* (a(3,1).*Kgen1 + a(3,2).*Kgen2 );
 
% Calculate system voltages
U2 = SolveNetwork(Xgen2, Pgen, Ly, Uy, Py, gbus, genmodel);

% Calculate machine currents and power
[Id2,Iq2,Pe2] = MachineCurrents(Xgen2, Pgen, U2(gbus), genmodel);    
   
% Update variables that have changed
Vexc2 = abs(U2(gbus));
Vgen2 = [Id2,Iq2,Pe2];  
Vgov2 = [Xgen2(:,2)];

%% K3

% EXCITERS
Kexc3 = Exciter(Xexc2, Pexc, Vexc2, excmodel);
Xexc3 = Xexc0 + stepsize.* (a(4,1).*Kexc1 + a(4,2).*Kexc2 + a(4,3).*Kexc3);

% GOVERNORS
Kgov3 = Governor(Xgov2, Pgov, Vgov2, govmodel);
Xgov3 = Xgov0 + stepsize.* (a(4,1).*Kgov1 + a(4,2).*Kgov2 + a(4,3).*Kgov3);
        
% GENERATORS
Kgen3 = Generator(Xgen2, Xexc3, Xgov3, Pgen, Vgen2, genmodel);
Xgen3 = Xgen0 + stepsize.* (a(4,1).*Kgen1 + a(4,2).*Kgen2 + a(4,3).*Kgen3);
 
% Calculate system voltages
U3 = SolveNetwork(Xgen3, Pgen, Ly, Uy, Py, gbus, genmodel);

% Calculate machine currents and power
[Id3,Iq3,Pe3] = MachineCurrents(Xgen3, Pgen, U3(gbus), genmodel);    
   
% Update variables that have changed
Vexc3 = abs(U3(gbus));
Vgen3 = [Id3,Iq3,Pe3];  
Vgov3 = [Xgen3(:,2)];


%% K4

% EXCITERS
Kexc4 = Exciter(Xexc3, Pexc, Vexc3, excmodel);
Xexc4 = Xexc0 + stepsize.* (b(1).*Kexc1 + b(2).*Kexc2 + b(3).*Kexc3 + b(4).*Kexc4);

% GOVERNORS
Kgov4 = Governor(Xgov3, Pgov, Vgov3, govmodel);
Xgov4 = Xgov0 + stepsize.* (b(1).*Kgov1 + b(2).*Kgov2 + b(3).*Kgov3 + b(4).*Kgov4);
        
% GENERATORS
Kgen4 = Generator(Xgen3, Xexc4, Xgov4, Pgen, Vgen3, genmodel);
Xgen4 = Xgen0 + stepsize.* (b(1).*Kgen1 + b(2).*Kgen2 + b(3).*Kgen3 + b(4).*Kgen4);
   
% Calculate system voltages
U4 = SolveNetwork(Xgen4, Pgen, Ly, Uy, Py, gbus, genmodel);

% Calculate machine currents and power
[Id4,Iq4,Pe4] = MachineCurrents(Xgen4, Pgen, U4(gbus), genmodel);    
   
% Update variables that have changed
Vexc4 = abs(U4(gbus));
Vgen4 = [Id4,Iq4,Pe4];  
Vgov4 = [Xgen4(:,2)];

%% Update
U0 = U4;
    
Vgen0 = Vgen4;
Vgov0 = Vgov4;
Vexc0 = Vexc4;
    
Xgen0 = Xgen4; 
Xexc0 = Xexc4;
Xgov0 = Xgov4;      
    
Pgen0 = Pgen; 
Pexc0 = Pexc;
Pgov0 = Pgov;       

t = t0;

return;