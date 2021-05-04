function [Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, U0, errest, failed, t, stepsize] = RungeKuttaHighamHall(t0, Xgen0, Pgen, Vgen0, Xexc0, Pexc, Vexc0, Xgov0, Pgov, Vgov0, U0, Ly, Uy, Py, gbus, genmodel, excmodel, govmodel, tol, maxstepsize, stepsize)

% [Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, U0
% errest, failed, t, stepsize] = RungeKuttaHighamHall(t0, Xgen0, Pgen
% Vgen0, Xexc0, Pexc, Vexc0, Xgov0, Pgov, Vgov0, U0, Y, gbus, faulted_bus
% theta0, genmodel, excmodel, govmodel, tol, maxstepsize, stepsize)
%  
% Runge-Kutta Higham and Hall ODE solver
 
% MatDyn
% Copyright (C) 2009 Stijn Cole
% Katholieke Universiteit Leuven
% Dept. Electrical Engineering (ESAT), Div. ELECTA
% Kasteelpark Arenberg 10
% 3001 Leuven-Heverlee, Belgium

%%
% Init
accept = false;
facmax = 4;
failed=0;


%% Runge-Kutta coefficients

% c = [0 2/9 1/3 1/2 3/5 1 1]; not used
a = [0 0 0 0 0 0
     2/9 0 0 0 0 0
     1/12 1/4 0 0 0 0
     1/8 0 3/8 0 0 0
     91/500 -27/100 78/125 8/125 0 0
     -11/20 27/20 12/5 -36/5 5 0
     1/12 0 27/32 -4/3 125/96 5/48];
b1 = [1/12 0 27/32 -4/3 125/96 5/48 0];
b2 = [2/15 0 27/80 -2/15 25/48 1/24 1/10];

%%
i=0;
while accept == false
    i=i+1;

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
Xexc4 = Xexc0 + stepsize.* (a(5,1).*Kexc1 + a(5,2).*Kexc2 + a(5,3).*Kexc3 + a(5,4).*Kexc4);

% GOVERNORS
Kgov4 = Governor(Xgov3, Pgov, Vgov3, govmodel);
Xgov4 = Xgov0 + stepsize.* (a(5,1).*Kgov1 + a(5,2).*Kgov2 + a(5,3).*Kgov3 + a(5,4).*Kgov4);
        
% GENERATORS
Kgen4 = Generator(Xgen3, Xexc4, Xgov4, Pgen, Vgen3, genmodel);
Xgen4 = Xgen0 + stepsize.* (a(5,1).*Kgen1 + a(5,2).*Kgen2 + a(5,3).*Kgen3 + a(5,4).*Kgen4);
   
% Calculate system voltages
U4 = SolveNetwork(Xgen4, Pgen, Ly, Uy, Py, gbus, genmodel);

% Calculate machine currents and power
[Id4,Iq4,Pe4] = MachineCurrents(Xgen4, Pgen, U4(gbus), genmodel);    
   
% Update variables that have changed
Vexc4 = abs(U4(gbus));
Vgen4 = [Id4,Iq4,Pe4];  
Vgov4 = [Xgen4(:,2)];


%% K5

% EXCITERS
Kexc5 = Exciter(Xexc4, Pexc, Vexc4, excmodel);
Xexc5 = Xexc0 + stepsize.*(a(6,1).*Kexc1 + a(6,2).*Kexc2 + a(6,3).*Kexc3 + a(6,4).*Kexc4 + a(6,5).*Kexc5);

% GOVERNORS
Kgov5 = Governor(Xgov4, Pgov, Vgov4, govmodel);
Xgov5 = Xgov0 + stepsize.*(a(6,1).*Kgov1 + a(6,2).*Kgov2 + a(6,3).*Kgov3 + a(6,4).*Kgov4 + a(6,5).*Kgov5);

% GENERATORS
Kgen5 = Generator(Xgen4, Xexc5, Xgov5, Pgen, Vgen4, genmodel);
Xgen5 = Xgen0 + stepsize.*(a(6,1).*Kgen1 + a(6,2).*Kgen2 + a(6,3).*Kgen3 + a(6,4).*Kgen4 + a(6,5).*Kgen5);

    
% Calculate system voltages
U5 = SolveNetwork(Xgen5, Pgen, Ly, Uy, Py, gbus, genmodel);

% Calculate machine currents and power
[Id5,Iq5,Pe5] = MachineCurrents(Xgen5, Pgen, U5(gbus), genmodel);    
   
% Update variables that have changed
Vexc5 = abs(U5(gbus));
Vgen5 = [Id5,Iq5,Pe5];  
Vgov5 = [Xgen5(:,2)];

%% K6

% EXCITERS
Kexc6 = Exciter(Xexc5, Pexc, Vexc5, excmodel);
Xexc6 = Xexc0 + stepsize.*(a(7,1).*Kexc1 + a(7,2).*Kexc2 + a(7,3).*Kexc3 + a(7,4).*Kexc4 + a(7,5).*Kexc5 + a(7,6).*Kexc6);

% GOVERNORS
Kgov6 = Governor(Xgov5, Pgov, Vgov5, govmodel);
Xgov6 = Xgov0 + stepsize.*(a(7,1).*Kgov1 + a(7,2).*Kgov2 + a(7,3).*Kgov3 + a(7,4).*Kgov4 + a(7,5).*Kgov5 + a(7,6).*Kgov6);

% GENERATORS
Kgen6 = Generator(Xgen5, Xexc6, Xgov6, Pgen, Vgen5, genmodel);
Xgen6 = Xgen0 + stepsize.*(a(7,1).*Kgen1 + a(7,2).*Kgen2 + a(7,3).*Kgen3 + a(7,4).*Kgen4 + a(7,5).*Kgen5 + a(7,6).*Kgen6);

    
% Calculate system voltages
U6 = SolveNetwork(Xgen6, Pgen, Ly, Uy, Py, gbus, genmodel);

% Calculate machine currents and power
[Id6,Iq6,Pe6] = MachineCurrents(Xgen6, Pgen, U6(gbus), genmodel);    
   
% Update variables that have changed
Vexc6 = abs(U6(gbus));
Vgen6 = [Id6,Iq6,Pe6];  
Vgov6 = [Xgen6(:,2)];


%% K7

% EXCITERS
Kexc7 = Exciter(Xexc6, Pexc, Vexc6, excmodel);
Xexc7 = Xexc0 + stepsize.*(b1(1).*Kexc1 + b1(2).*Kexc2 + b1(3).*Kexc3 + b1(4).*Kexc4 + b1(5).*Kexc5 + b1(6).*Kexc6 + b1(7).*Kexc7);

% GOVERNORS
Kgov7 = Governor(Xgov6, Pgov, Vgov6, govmodel);
Xgov7 = Xgov0 + stepsize.*(b1(1).*Kgov1 + b1(2).*Kgov2 + b1(3).*Kgov3 + b1(4).*Kgov4 + b1(5).*Kgov5 + b1(6).*Kgov6 + b1(7).*Kgov7);

% GENERATORS
Kgen7 = Generator(Xgen6, Xexc7, Xgov7, Pgen, Vgen6, genmodel);
Xgen7 = Xgen0 + stepsize.*(b1(1).*Kgen1 + b1(2).*Kgen2 + b1(3).*Kgen3 + b1(4).*Kgen4 + b1(5).*Kgen5 + b1(6).*Kgen6 + b1(7).*Kgen7);

    
% Calculate system voltages
U7 = SolveNetwork(Xgen7, Pgen, Ly, Uy, Py, gbus, genmodel);

% Calculate machine currents and power
[Id7,Iq7,Pe7] = MachineCurrents(Xgen7, Pgen, U7(gbus), genmodel);    
   
% Update variables that have changed
Vexc7 = abs(U7(gbus));
Vgen7 = [Id7,Iq7,Pe7];  
Vgov7 = [Xgen7(:,2)];

     
%% Second, higher order solution

Xexc72 = Xexc0 + stepsize.*(b2(1).*Kexc1 + b2(2).*Kexc2 + b2(3).*Kexc3 + b2(4).*Kexc4 + b2(5).*Kexc5 + b2(6).*Kexc6 + b2(7).*Kexc7);
Xgov72 = Xgov0 + stepsize.*(b2(1).*Kgov1 + b2(2).*Kgov2 + b2(3).*Kgov3 + b2(4).*Kgov4 + b2(5).*Kgov5 + b2(6).*Kgov6 + b2(7).*Kgov7);
Xgen72 = Xgen0 + stepsize.*(b2(1).*Kgen1 + b2(2).*Kgen2 + b2(3).*Kgen3 + b2(4).*Kgen4 + b2(5).*Kgen5 + b2(6).*Kgen6 + b2(7).*Kgen7);

%% Error estimate

Xexc = abs((Xexc72-Xexc7)');
Xgov = abs((Xgov72-Xgov7)');
Xgen = abs((Xgen72-Xgen7)');
errest = max( [max(max(Xexc)) max(max(Xgov)) max(max(Xgen)) ]);

if errest < eps
    errest = eps;
end
q = 0.84*(tol/errest)^(1/4);

if (errest < tol)
	accept = true;
    
	U0 = U7;
    
    Vgen0 = Vgen7;
    Vgov0 = Vgov7;
    Vexc0 = Vexc7;
    
    Xgen0 = Xgen7; 
    Xexc0 = Xexc7;
    Xgov0 = Xgov7;      
    
    Pgen0 = Pgen; 
    Pexc0 = Pexc;
    Pgov0 = Pgov;      

    t=t0;
else   
    failed=failed+1;    
    facmax = 1;
    
    t=t0;
    Pgen0 = Pgen; 
    Pexc0 = Pexc;
    Pgov0 = Pgov;  
    
    stepsize = min(max(q, 0.1), facmax)*stepsize;
    return
end

stepsize = min(max(q, 0.1), facmax)*stepsize;
    
if (stepsize > maxstepsize)
	stepsize = maxstepsize;
end

    
end
 
return;