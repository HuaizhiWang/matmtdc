function [Angles,Speeds,Eq_tr,Ed_tr,Efd,PM,Voltages,Stepsize,Errest,Time] = rundyn(casefile_pf, casefile_dyn, casefile_ev, mdopt)

% [Angles,Speeds,Eq_tr,Ed_tr,Efd,PM,Voltages,Stepsize,Errest,Time] =
% rundyn(casefile_pf, casefile_dyn, casefile_ev, mdopt)
% 
% Runs dynamic simulation
% 
% INPUTS
% casefile_pf = m-file with power flow data
% casefile_dyn = m-file with dynamic data
% casefile_ev = m-file with event data
% mdopt = options vector
% 
% OUTPUTS
% Angles = generator angles
% Speeds = generator speeds
% Eq_tr = q component of transient voltage behind reactance
% Ed_tr = d component of transient voltage behind reactance
% Efd = Excitation voltage
% PM = mechanical power
% Voltages = bus voltages
% Stepsize = step size integration method
% Errest = estimation of integration error
% Failed = failed steps
% Time = time points
 
% MatDyn
% Copyright (C) 2009 Stijn Cole
% Katholieke Universiteit Leuven
% Dept. Electrical Engineering (ESAT), Div. ELECTA
% Kasteelpark Arenberg 10
% 3001 Leuven-Heverlee, Belgium

%% Begin timing
tic;
%% Add subdirectories to path

addpath([cd '/Solvers/']);
addpath([cd '/Models/Generators/']);
addpath([cd '/Models/Exciters/']);
addpath([cd '/Models/Governors/']);
addpath([cd '/Models/MMC-HVDC/']);
addpath([cd '/Cases/Powerflow/']);
addpath([cd '/Cases/Dynamic/']);
addpath([cd '/Cases/Events/']);
addpath([cd '/Matpower/']);


%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

%% Options
if nargin < 4
	mdopt = Mdoption;   
end
method = mdopt(1);
tol = mdopt(2);
minstepsize = mdopt(3);
maxstepsize = mdopt(4);
output = mdopt(5);
plots = mdopt(6);


%% Load all data

% Load dynamic simulation data
if output; disp('> Loading dynamic simulation data...'); end
global freq
[freq,stepsize,stoptime] = Loaddyn(casefile_dyn);

% Load generator data
Pgen0 = Loadgen(casefile_dyn, output);

% Load exciter data
Pexc0 = Loadexc(casefile_dyn);

% Load governor data
Pgov0 = Loadgov(casefile_dyn);

% Load mmc data
Pmmc0 = Loadmmc(casefile_dyn);

% Load circulating current controller parameters
Pcir0 = Loadcir(casefile_dyn);

% Load ouput current controller parameters
Pouc0 = Loadouc(casefile_dyn, Pmmc0);

% Load DC voltage controller parameters
Pdcc0 = Loaddcc(casefile_dyn);

% Load event data
if ~isempty(casefile_ev)
    [event,buschange,linechange] = Loadevents(casefile_ev);
else
    event=[];
end

genmodel = Pgen0(:,1);
excmodel = Pgen0(:,2);
govmodel = Pgen0(:,3);

%% Initialization: Power Flow 
Type1 = 1; % inverter of MMC-HVDC
Type2 = 2; % rectifier of MMC-HVDC
% Power flow options
mpopt=mpoption;
mpopt(31)=0;
mpopt(32)=0;
% Run power flow
[baseMVA, bus, gen, branch, success] = runpf(casefile_pf,mpopt);
if ~success
    fprintf('> Error: Power flow did not converge. Exiting...\n')
    return;
else
    if output; fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b> Power flow converged\n'); end
end

U0=bus(:,VM).*(cos(bus(:,VA)*pi/180) + 1j.*sin(bus(:,VA)*pi/180));
U00=U0;
% Get generator info
on = find(gen(:, GEN_STATUS) > 0);     %% which generators are on?
gbus = gen(on, GEN_BUS);               %% what buses are they at?
ngen = length(gbus);

% Calculate steady-state power flow for MMC-HVDC
MMCflow0 = BranchFlow(baseMVA, bus, branch, U00);

%% Construct augmented Ybus 
if output; disp('> Constructing augmented admittance matrix...'); end
Pl=bus(:,PD)./baseMVA;                  %% load power
Ql=bus(:,QD)./baseMVA;

xd_tr = zeros(ngen,1);
xd_tr(genmodel==2) = Pgen0(genmodel==2,8); % 4th order model: xd_tr column 8
xd_tr(genmodel==1) = Pgen0(genmodel==1,7); % classical model: xd_tr column 7

[Ly, Uy, Py] = AugYbus(baseMVA, bus, branch, xd_tr, gbus, Pl, Ql, U0);

%% MMC initial conditions
% Xmmci0, Vmmci0 for inverter
[Xmmci0, Vmmci0] = MMCInit(Pmmc0, MMCflow0, U00, Type1);  %第一处
% Xmmcr0, Vmmcr0 for rectifier
[Xmmcr0, Vmmcr0] = MMCInit(Pmmc0, MMCflow0, U00, Type2);
Prec = Vmmcr0(4,1)/1e6;

%% 更新2号发电机输出的有功无功
u8x = Vmmcr0(1,1)/Pmmc0(1,4);   %第二处
u8y = Vmmcr0(2,1)/Pmmc0(1,4);
u2x = real(U0(2,1));
u2y = imag(U0(2,1));
gen(2,2) = -Prec;
% gen(2,2) = (MMCflow0(1,1)*100e6+((MMCflow0(1,1)*100e6)/Pmmc0(1,3))^2*Pmmc0(1,14))/1e6; %Ps0-Ps0/Vd_ref*Rdc
gen(2,3) = (1-u2x*u8x-u2y*u8y)/0.06;
%% Output current controller initialization
% Xouci0, Vouci0 for inverter
[Xouci0, Vouci0] = OucInit(Xmmci0, Pmmc0, Vmmci0);
% Xoucr0, Voucr0 for inverter
[Xoucr0, Voucr0] = OucInit(Xmmcr0, Pmmc0, Vmmcr0);

%% Circulating current controller initialization
% Xciri0, Vciri0 for inverter
[Xciri0, Vciri0] = CirInit(Xmmci0, Pmmc0, Vmmci0, Vouci0, Pcir0, Type1);
% Xcirr0, Vcirr0 for rectifier
[Xcirr0, Vcirr0] = CirInit(Xmmcr0, Pmmc0, Vmmcr0, Voucr0, Pcir0, Type2);

% Xdcci0, Vdcci0 for inverter % 第三处
[Xdcci0, Vdcci0] = DccInit(Pdcc0, Xmmci0, Pmmc0, Vmmci0);
Vdccr0 = 0;

%% Calculate Initial machine state
if output; fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b> Calculating initial state...\n'); end
[Efd0, Xgen0] = GeneratorInit(Pgen0, U0(gbus), gen, baseMVA, genmodel);

omega0 = Xgen0(:,2);
U3 = SolveNetwork(Xgen0, Pgen0, Ly, Uy, Py, gbus, genmodel, Xmmci0, Pmmc0);

[Id0,Iq0,Pe0] = MachineCurrents(Xgen0, Pgen0, U0(gbus), genmodel);
Vgen0 = [Id0, Iq0, Pe0];

%% Exciter initial conditions

Vexc0 = [abs(U0(gbus))];
[Xexc0,Pexc0] = ExciterInit(Efd0, Pexc0, Vexc0, excmodel);


%% Governor initial conditions

Pm0 = Pe0;
[Xgov0, Pgov0] = GovernorInit(Pm0, Pgov0, omega0, govmodel);
Vgov0 = omega0;

%% Check Steady-state

Fexc0 = Exciter(Xexc0, Pexc0, Vexc0, excmodel);
Fgov0 = Governor(Xgov0, Pgov0, Vgov0, govmodel);
Fgen0 = Generator(Xgen0, Xexc0, Xgov0, Pgen0, Vgen0, genmodel);
Fmmci0 = MMC(Xmmci0, Pmmc0, Vmmci0, Vouci0, Vciri0);
Fmmcr0 = MMC(Xmmcr0, Pmmc0, Vmmcr0, Voucr0, Vcirr0);
Fouci0 = OucCon(Xouci0, Pouc0, Vouci0, Xmmci0, Pmmc0, Vmmci0, Vdcci0, MMCflow0, Type1);
Foucr0 = OucCon(Xoucr0, Pouc0, Voucr0, Xmmcr0, Pmmc0, Vmmcr0, Vdccr0, MMCflow0, Type2);
% Ciri0 = CirCon(Xciri0, Pcir0, Vciri0, Xmmci0, Pmmc0,MMCflow0,Wff_sumi,Wff_deltai,Type1,t);
% Cirr0 = CirCon(Xcirr0, Pcir0, Vcirr0, Xmmcr0, Pmmc0,MMCflow0,Wff_sumr,Wff_deltar,Type2,t);
Dcc0 = DvotCon(Xdcci0, Pdcc0, Vdcci0, Vmmci0);
% Check Generator Steady-state
if sum(sum(abs(Fgen0))) > 1e-6
    fprintf('> Error: Generator not in steady-state\n> Exiting...\n')
    return;
end
% Check Exciter Steady-state
if sum(sum(abs(Fexc0))) > 1e-6
	fprintf('> Error: Exciter not in steady-state\n> Exiting...\n')
	return;
end
% Check Governor Steady-state
if sum(sum(abs(Fgov0))) > 1e-6
    fprintf('> Error: Governor not in steady-state\n> Exiting...\n')
    return;
end

if output; fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b> System in steady-state\n'); end


%% Initialization of main stability loop
t=0;
errest=0;
failed=0;
eulerfailed = 0;

if method==3 || method==4
    stepsize = minstepsize;
end

if ~output
    fprintf('                   ')
end

ev=1;
eventhappened = false;
i=0;

%% Allocate memory for variables

if output; fprintf('> Allocate memory..'); end
chunk = 5000;

Time = zeros(chunk,1); Time(1,:) = t;
Errest = zeros(chunk,1); Errest(1,:) = errest;
Stepsize = zeros(chunk,1); Stepsize(1,:) = stepsize;

% System variables
Voltages = zeros(chunk, length(U0)); Voltages(1,:) = U0.';

% Generator
Angles = zeros(chunk,ngen); Angles(1,:) = Xgen0(:,1).*180./pi;
Speeds = zeros(chunk,ngen); Speeds(1,:) = Xgen0(:,2)./(2.*pi.*freq);
Eq_tr = zeros(chunk,ngen); Eq_tr(1,:) = Xgen0(:,3);
Ed_tr = zeros(chunk,ngen); Ed_tr(1,:) = Xgen0(:,4);

% Exciter and governor
Efd = zeros(chunk,ngen); Efd(1,:) = Efd0(:,1);  
PM = zeros(chunk,ngen); PM(1,:) = Pm0(:,1);

% MMC-HVDC
% Xmmc0 = [ic0; vcu0; vcl0; Is_alpha0; Is_beta0; Vd0]
Icsi = zeros(chunk,1); Icsi(1,1) = Xmmci0(1,1); Icsr = zeros(chunk,1); Icsr(1,1) = Xmmcr0(1,1);
Vcusi = zeros(chunk,1); Vcusi(1,1) = Xmmci0(2,1); Vcusr = zeros(chunk,1); Vcusr(1,1) = Xmmcr0(2,1);
Vclsi = zeros(chunk,1); Vclsi(1,1) = Xmmci0(3,1); Vclsr = zeros(chunk,1); Vclsr(1,1) = Xmmcr0(3,1);
Is_alphasi = zeros(chunk,1); Is_alphasi(1,1) = Xmmci0(4,1); Is_alphasr = zeros(chunk,1); Is_alphasr(1,1) = Xmmcr0(4,1);
Is_betasi = zeros(chunk,1); Is_betasi(1,1) = Xmmci0(5,1); Is_betasr = zeros(chunk,1); Is_betasr(1,1) = Xmmcr0(5,1);
Vdsi = zeros(chunk,1); Vdsi(1,1) = Xmmci0(6,1); Vdsr = zeros(chunk,1); Vdsr(1,1) = Xmmcr0(6,1);
% Vmmc0 = [Va_alpha; Va_beta; is0; Ps; Qs; id0]
Vai = zeros(chunk,1); Vai(1,1) = sqrt(Vmmci0(1,1)^2+Vmmci0(2,1)^2);
Var = zeros(chunk,1); Var(1,1) = sqrt(Vmmcr0(1,1)^2+Vmmcr0(2,1)^2);
Issi = zeros(chunk,1); Issi(1,1) = Vmmci0(3,1); Issr = zeros(chunk,1); Issr(1,1) = Vmmcr0(3,1);
Pssi = zeros(chunk,1); Pssi(1,1) = Vmmci0(4,1); Pssr = zeros(chunk,1); Pssr(1,1) = Vmmcr0(4,1);
Qssi = zeros(chunk,1); Qssi(1,1) = Vmmci0(5,1); Qssr = zeros(chunk,1); Qssr(1,1) = Vmmcr0(5,1);
% Xouc0 = [VR_alpha; VR_beta]
VR_alphasi = zeros(chunk,1); VR_alphasi(1,1) = Xouci0(1,1); VR_alphasr = zeros(chunk,1); VR_alphasr(1,1) = Xoucr0(1,1);
VR_betasi = zeros(chunk,1); VR_betasi(1,1) = Xouci0(2,1); VR_betasr = zeros(chunk,1); VR_betasr(1,1) = Xoucr0(2,1);
% Vouc0 = [Vs_alpha; Vs_beta; vs_ref; vs]
vs_refsi = zeros(chunk,1); vs_refsi(1,1) = Vouci0(3,1); vs_refsr = zeros(chunk,1); vs_refsr(1,1) = Voucr0(3,1);
Vssi = zeros(chunk,1); Vssi(1,1) = Vouci0(4,1); Vssr = zeros(chunk,1); Vssr(1,1) = Voucr0(4,1);
% Xcir0 = [xc10; xc20]
Xc1si = zeros(chunk,1); Xc1si(1,1) = Xciri0(1,1); Xc1sr = zeros(chunk,1); Xc1sr(1,1) = Xcirr0(1,1);
Xc2si = zeros(chunk,1); Xc2si(1,1) = Xciri0(2,1); Xc2sr = zeros(chunk,1); Xc2sr(1,1) = Xcirr0(2,1);
% Vcir0 = [vc_ref; vc; nu; nl];
vc_refsi = zeros(chunk,1); vc_refsi(1,1) = Vciri0(1,1); vc_refsr = zeros(chunk,1); vc_refsr(1,1) = Vcirr0(1,1);
Vcsi = zeros(chunk,1); Vcsi(1,1) = Vciri0(2,1); Vcsr = zeros(chunk,1); Vcsr(1,1) = Vcirr0(2,1);
nusi = zeros(chunk,1); nusi(1,1) = Vciri0(3,1); nusr = zeros(chunk,1); nusr(1,1) = Vcirr0(3,1);
nlsi = zeros(chunk,1); nlsi(1,1) = Vciri0(4,1); nlsr = zeros(chunk,1); nlsr(1,1) = Vcirr0(4,1);
W_sumi = zeros(chunk,1); W_deltai = zeros(chunk,1);
W_sumr = zeros(chunk,1); W_deltar = zeros(chunk,1);
Wf_sumi = zeros(chunk,1); Wf_deltai = zeros(chunk,1);
Wf_sumr = zeros(chunk,1); Wf_deltar = zeros(chunk,1);
Wff_sumi = zeros(chunk,1); Wff_deltai = zeros(chunk,1);
Wff_sumr = zeros(chunk,1); Wff_deltar = zeros(chunk,1);
Wsumi = 0;Wdeltai = 0;Wsumr = 0;Wdeltar = 0;
% first filter
wl = 2*pi*50;
num1 = [1 0 wl^2];den1 = [1 2*wl wl^2];
sys1 = tf(num1,den1);dsys1 =c2d(sys1,stepsize,'zoh');
b1= cell2mat(dsys1.Numerator);a1 = cell2mat(dsys1.Denominator);
% second filter
num2 = [1 0 (2*wl)^2];den2 = [1 2*2*wl (2*wl)^2];
sys2 = tf(num2,den2);dsys2 =c2d(sys2,stepsize,'zoh');
b2= cell2mat(dsys2.Numerator);a2 = cell2mat(dsys2.Denominator);
% Xdcc0 = xdc0
Xdcsi = zeros(chunk,1); Xdcsi(1,1) = Xdcci0;

%% Main stability loop
while t < stoptime + stepsize

    %% Output    
    i=i+1;
    if mod(i,45)==0 && output
       fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b> %6.2f%% completed', t/stoptime*100)
    end  
    
    %% Numerical Method
    switch method
        case 1    % 第四处
            [Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, Xmmci0, Xmmcr0, Pmmc0, Vmmci0, Vmmcr0, Xouci0, Xoucr0, Pouc0, Vouci0, Voucr0, Xciri0, Xcirr0, Pcir0, Vciri0, Vcirr0, Xdcci0, Pdcc0, Vdcci0, U0, t, newstepsize] = ModifiedEuler(t, Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, Xmmci0, Xmmcr0, Pmmc0, Vmmci0, Vmmcr0, Xouci0, Xoucr0, Pouc0, Vouci0, Voucr0, Xciri0, Xcirr0, Pcir0, Vciri0, Vcirr0, Xdcci0, Pdcc0, Vdcci0, MMCflow0, Ly, Uy, Py, gbus, genmodel, excmodel, govmodel, stepsize, U0,Wsumi,Wdeltai,Wsumr,Wdeltar);            
        case 2
            [Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, U0, t, newstepsize] = RungeKutta(t, Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, Ly, Uy, Py, gbus, genmodel, excmodel, govmodel, stepsize);
        case 3
            [Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, U0, errest, failed, t, newstepsize] = RungeKuttaFehlberg(t, Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, U0, Ly, Uy, Py, gbus, genmodel, excmodel, govmodel, tol, maxstepsize, stepsize);
        case 4          
            [Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, U0, errest, failed, t, newstepsize] = RungeKuttaHighamHall(t, Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, U0, Ly, Uy, Py, gbus, genmodel, excmodel, govmodel, tol, maxstepsize, stepsize);                 
        case 5
            [Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, U0, t, eulerfailed, newstepsize] = ModifiedEuler2(t, Xgen0, Pgen0, Vgen0, Xexc0, Pexc0, Vexc0, Xgov0, Pgov0, Vgov0, Ly, Uy, Py, gbus, genmodel, excmodel, govmodel, stepsize);
    end
        
    if eulerfailed
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b> Error: No solution found. Try lowering tolerance or increasing maximum number of iterations in ModifiedEuler2. Exiting... \n')
    	return;
    end
    
    if failed
        t = t-stepsize;
    end
    
    % End exactly at stop time
    if t + newstepsize > stoptime
        newstepsize = stoptime - t;
    end;
    
    
    %% Allocate new memory chunk if matrices are full
    if i>size(Time,1)
        Stepsize = [Stepsize; zeros(chunk,1)];Errest = [Errest; zeros(chunk,1)];Time = [Time; zeros(chunk,1)];
        Voltages = [Voltages; zeros(chunk,length(U0))];Efd = [Efd; zeros(chunk,ngen)];PM = [PM; zeros(chunk,ngen)];
        Angles=[Angles;zeros(chunk,ngen)];Speeds=[Speeds;zeros(chunk,ngen)];Eq_tr=[Eq_tr;zeros(chunk,ngen)];
        Ed_tr=[Ed_tr;zeros(chunk,ngen)]; Icsi = [Icsi; zeros(chunk,1)]; Icsr = [Icsr; zeros(chunk,1)];
        Vcusi=[Vcusi; zeros(chunk,1)]; Vcusr=[Vcusr; zeros(chunk,1)]; Vclsi=[Vclsi; zeros(chunk,1)];
        Vclsr=[Vclsr; zeros(chunk,1)]; Is_alphasi=[Is_alphasi;zeros(chunk,1)]; Is_alphasr=[Is_alphasr;zeros(chunk,1)];
        Is_betasi=[Is_betasi;zeros(chunk,1)]; Is_betasr=[Is_betasr;zeros(chunk,1)];
        Vdsi=[Vdsi;zeros(chunk,1)]; Vdsr=[Vdsr;zeros(chunk,1)]; Vai=[Vai;zeros(chunk,1)];
        Var=[Var;zeros(chunk,1)]; Issi=[Issi;zeros(chunk,1)]; Issr=[Issr;zeros(chunk,1)];
        Pssi=[Pssi;zeros(chunk,1)]; Pssr=[Pssr;zeros(chunk,1)]; Qssi=[Qssi;zeros(chunk,1)];
        Qssr=[Qssr;zeros(chunk,1)]; VR_alphasi=[VR_alphasi;zeros(chunk,1)]; VR_alphasr=[VR_alphasr;zeros(chunk,1)];
        VR_betasi=[VR_betasi;zeros(chunk,1)]; VR_betasr=[VR_betasr;zeros(chunk,1)];
        vs_refsi=[vs_refsi;zeros(chunk,1)]; vs_refsr=[vs_refsr;zeros(chunk,1)];
        Vssi=[Vssi;zeros(chunk,1)]; Vssr=[Vssr;zeros(chunk,1)]; Xc1si=[Xc1si;zeros(chunk,1)];
        Xc1sr=[Xc1sr;zeros(chunk,1)]; Xc2si=[Xc2si;zeros(chunk,1)]; Xc2sr=[Xc2sr;zeros(chunk,1)];
        vc_refsi=[vc_refsi;zeros(chunk,1)]; vc_refsr=[vc_refsr;zeros(chunk,1)];
        Vcsi=[Vcsi;zeros(chunk,1)]; Vcsr=[Vcsr;zeros(chunk,1)]; nusi=[nusi;zeros(chunk,1)];
        nusr=[nusr;zeros(chunk,1)]; nlsi=[nlsi;zeros(chunk,1)]; nlsr=[nlsr;zeros(chunk,1)];
        Xdcsi=[Xdcsi;zeros(chunk,1)];W_sumi = [W_sumi;zeros(chunk,1)]; W_deltai = [W_deltai;zeros(chunk,1)];
        W_sumr =[W_sumr;zeros(chunk,1)];W_deltar =[W_deltar;zeros(chunk,1)];
        Wf_sumi =[Wf_sumi;zeros(chunk,1)]; Wf_deltai =[Wf_deltai;zeros(chunk,1)];
        Wf_sumr =[Wf_sumr;zeros(chunk,1)]; Wf_deltar =[Wf_deltar;zeros(chunk,1)];
        Wff_sumi =[Wff_sumi;zeros(chunk,1)]; Wff_deltai =[Wff_deltai;zeros(chunk,1)];
        Wff_sumr =[Wff_sumr;zeros(chunk,1)]; Wff_deltar =[Wff_deltar;zeros(chunk,1)];
    end
    
    %% Save values
    Stepsize(i,:) = stepsize.';
    Errest(i,:) = errest.';
    Time(i,:) = t;    
    Voltages(i,:) = U0.';
    % exc
    Efd(i,:) = Xexc0(:,1).*(genmodel>1); % Set Efd to zero when using classical generator model    
    % gov
    PM(i,:) = Xgov0(:,1);   
    % gen
	Angles(i,:) = Xgen0(:,1).*180./pi;
    Speeds(i,:) = Xgen0(:,2)./(2.*pi.*freq);
    Eq_tr(i,:) = Xgen0(:,3);
    Ed_tr(i,:) = Xgen0(:,4);
    % MMC-HVDC
    Icsi(i,1) = Xmmci0(1,1); Icsr(i,1) = Xmmcr0(1,1);
    Vcusi(i,1) = Xmmci0(2,1); Vcusr(i,1) = Xmmcr0(2,1);
    Vclsi(i,1) = Xmmci0(3,1); Vclsr(i,1) = Xmmcr0(3,1);
    Is_alphasi(i,1) = Xmmci0(4,1); Is_alphasr(i,1) = Xmmcr0(4,1);
    Is_betasi(i,1) = Xmmci0(5,1); Is_betasr(i,1) = Xmmcr0(5,1);
    Vdsi(i,1) = Xmmci0(6,1); Vdsr(i,1) = Xmmcr0(6,1);
    Vai(i,1) = sqrt(Vmmci0(1,1)^2+Vmmci0(2,1)^2);
    Var(i,1) = sqrt(Vmmcr0(1,1)^2+Vmmcr0(2,1)^2);
    Issi(i,1) = Vmmci0(3,1); Issr(i,1) = Vmmcr0(3,1);
    Pssi(i,1) = Vmmci0(4,1); Pssr(i,1) = Vmmcr0(4,1);
    Qssi(i,1) = Vmmci0(5,1); Qssr(i,1) = Vmmcr0(5,1);
    VR_alphasi(i,1) = Xouci0(1,1); VR_alphasr(i,1) = Xoucr0(1,1);
    VR_betasi(i,1) = Xouci0(2,1); VR_betasr(i,1) = Xoucr0(2,1);
    vs_refsi(i,1) = Vouci0(3,1); vs_refsr(i,1) = Voucr0(3,1);
    Vssi(i,1) = Vouci0(4,1); Vssr(i,1) = Voucr0(4,1);
    Xc1si(i,1) = Xciri0(1,1); Xc1sr(i,1) = Xcirr0(1,1);
    Xc2si(i,1) = Xciri0(2,1); Xc2sr(i,1) = Xcirr0(2,1);
    Xc3si(i,1) = Xciri0(3,1); Xc3sr(i,1) = Xcirr0(3,1);
    Xc4si(i,1) = Xciri0(4,1); Xc4sr(i,1) = Xcirr0(4,1);
    vc_refsi(i,1) = Vciri0(1,1); vc_refsr(i,1) = Vcirr0(1,1);
    Vcsi(i,1) = Vciri0(2,1); Vcsr(i,1) = Vcirr0(2,1);
    nusi(i,1) = Vciri0(3,1); nusr(i,1) = Vcirr0(3,1);
    nlsi(i,1) = Vciri0(4,1); nlsr(i,1) = Vcirr0(4,1);
    Xdcsi(i,1) = Xdcci0;
    %% Adapt step size if event will occur in next step
    if ~isempty(event) && ev <= size(event,1) && (method == 3 || method == 4)      
        if t + newstepsize >= event(ev,1)
            if event(ev,1) - t < newstepsize
                newstepsize = event(ev,1) - t;
            end
        end
    end    
    %% Check for events
    if ~isempty(event) && ev <= size(event,1)           
        for k=ev:size(event,1) % cycle through all events ..   
            if abs(t-event(ev,1))>1e-10 ||  ev > size(event,1) %.. that happen on time t               
                break;
            else
                eventhappened = true;
            end
                switch event(ev,2)
                    case 1
                        bus(buschange(ev,2),buschange(ev,3)) = buschange(ev,4);
                    case 2
                        branch(linechange(ev,2),linechange(ev,3)) = linechange(ev,4);                        
                end 
                ev=ev+1;
        end          
            
        if eventhappened
            % Refactorise
            [Ly, Uy, Py] = AugYbus(baseMVA, bus, branch, xd_tr, gbus, bus(:,PD)./baseMVA, bus(:,QD)./baseMVA, U00);
            U0 = SolveNetwork(Xgen0, Pgen0, Ly, Uy, Py, gbus, genmodel, Xmmci0, Pmmc0);
            % MMC-HVDC
            Vmmci0 = RenewVmmci(Xmmci0, Xmmcr0, Pmmc0, U0, t);  %第五处
            [Vmmcr0, Ug2] = RenewVmmcr(Xmmci0, Xmmcr0, Pmmc0, Vmmcr0, t);  %第六处
            U0(2,1) = Ug2;
            [Id0,Iq0,Pe0] = MachineCurrents(Xgen0, Pgen0, U0(gbus), genmodel);
            Vgen0 = [Id0,Iq0,Pe0];
            Vexc0 = abs(U0(gbus));            
            % decrease stepsize after event occured
            if method==3 || method==4
                newstepsize = minstepsize;
            end            
            i=i+1; % if event occurs, save values at t- and t+
           %% Save values
            Stepsize(i,:) = stepsize.';
            Errest(i,:) = errest.';            
            Time(i,:) = t;    
            Voltages(i,:) = U0.';
            % exc
            Efd(i,:) = Xexc0(:,1).*(genmodel>1); % Set Efd to zero when using classical generator model    
            % gov
            PM(i,:) = Xgov0(:,1);   
            % gen
            Angles(i,:) = Xgen0(:,1).*180./pi;
            Speeds(i,:) = Xgen0(:,2)./(2.*pi.*freq);
            Eq_tr(i,:) = Xgen0(:,3);
            Ed_tr(i,:) = Xgen0(:,4);
            % MMC-HVDC
            Icsi(i,1) = Xmmci0(1,1); Icsr(i,1) = Xmmcr0(1,1);
            Vcusi(i,1) = Xmmci0(2,1); Vcusr(i,1) = Xmmcr0(2,1);
            Vclsi(i,1) = Xmmci0(3,1); Vclsr(i,1) = Xmmcr0(3,1);
            Is_alphasi(i,1) = Xmmci0(4,1); Is_alphasr(i,1) = Xmmcr0(4,1);
            Is_betasi(i,1) = Xmmci0(5,1); Is_betasr(i,1) = Xmmcr0(5,1);
            Vdsi(i,1) = Xmmci0(6,1); Vdsr(i,1) = Xmmcr0(6,1);
            Vai(i,1) = sqrt(Vmmci0(1,1)^2+Vmmci0(2,1)^2);
            Var(i,1) = sqrt(Vmmcr0(1,1)^2+Vmmcr0(2,1)^2);
            Issi(i,1) = Vmmci0(3,1); Issr(i,1) = Vmmcr0(3,1);
            Pssi(i,1) = Vmmci0(4,1); Pssr(i,1) = Vmmcr0(4,1);
            Qssi(i,1) = Vmmci0(5,1); Qssr(i,1) = Vmmcr0(5,1);
            VR_alphasi(i,1) = Xouci0(1,1); VR_alphasr(i,1) = Xoucr0(1,1);
            VR_betasi(i,1) = Xouci0(2,1); VR_betasr(i,1) = Xoucr0(2,1);
            vs_refsi(i,1) = Vouci0(3,1); vs_refsr(i,1) = Voucr0(3,1);
            Vssi(i,1) = Vouci0(4,1); Vssr(i,1) = Voucr0(4,1);
            Xc1si(i,1) = Xciri0(1,1); Xc1sr(i,1) = Xcirr0(1,1);
            Xc2si(i,1) = Xciri0(2,1); Xc2sr(i,1) = Xcirr0(2,1);
            Xc3si(i,1) = Xciri0(3,1); Xc3sr(i,1) = Xcirr0(3,1);
            Xc4si(i,1) = Xciri0(4,1); Xc4sr(i,1) = Xcirr0(4,1);
            vc_refsi(i,1) = Vciri0(1,1); vc_refsr(i,1) = Vcirr0(1,1);
            Vcsi(i,1) = Vciri0(2,1); Vcsr(i,1) = Vcirr0(2,1);
            nusi(i,1) = Vciri0(3,1); nusr(i,1) = Vcirr0(3,1);
            nlsi(i,1) = Vciri0(4,1); nlsr(i,1) = Vcirr0(4,1);
            Xdcsi(i,1) = Xdcci0;
            
            eventhappened = false;
        end
    end   
    % get sum-capacitor-voltage ripples
        W_sumi(i) = (Vcusi(i)^2 + Vclsi(i)^2); W_deltai(i) = (Vcusi(i)^2 - Vclsi(i)^2);
        W_sumr(i) = (Vcusr(i)^2 + Vclsr(i)^2); W_deltar(i) = (Vcusr(i)^2 - Vclsr(i)^2);

        % first filter
        if i <= 2
        Wf_sumi(1) = W_sumi(1);Wf_sumi(2) = -a1(2)*Wf_sumi(1)+b1(1)*W_sumi(2)+b1(2)*W_sumi(1);
        Wf_deltai(1) = W_deltai(1);Wf_deltai(2) = -a1(2)*Wf_deltai(1)+b1(1)*W_deltai(2)+b1(2)*W_deltai(1);
        Wf_sumr(1) = W_sumr(1);Wf_sumr(2) = -a1(2)*Wf_sumr(1)+b1(1)*W_sumr(2)+b1(2)*W_sumr(1);
        Wf_deltar(1) = W_deltar(1);Wf_deltar(2) = -a1(2)*Wf_deltar(1)+b1(1)*W_deltar(2)+b1(2)*W_deltar(1);
        end
        if i >= 3 
        Wf_sumi(i) = self_filter(Wf_sumi(i-1),Wf_sumi(i-2),W_sumi(i),W_sumi(i-1),W_sumi(i-2),a1,b1);
        Wf_deltai(i) = self_filter(Wf_deltai(i-1),Wf_deltai(i-2),W_deltai(i),W_deltai(i-1),W_deltai(i-2),a1,b1);
        Wf_sumr(i) = self_filter(Wf_sumr(i-1),Wf_sumr(i-2),W_sumr(i),W_sumr(i-1),W_sumr(i-2),a1,b1);
        Wf_deltar(i) = self_filter(Wf_deltar(i-1),Wf_deltar(i-2),W_deltar(i),W_deltar(i-1),W_deltar(i-2),a1,b1);
        end

        % second filter
        if i <= 2
        Wff_sumi(1) = Wf_sumi(1);Wff_sumi(2) = -a2(2)*Wff_sumi(1)+b2(1)*Wf_sumi(2)+b2(2)*Wf_sumi(1);
        Wff_deltai(1) = Wf_deltai(1);Wff_deltai(2) = -a2(2)*Wff_deltai(1)+b2(1)*Wf_deltai(2)+b2(2)*Wf_deltai(1);
        Wff_sumr(1) = Wf_sumr(1);Wf_sumr(2) = -a2(2)*Wff_sumr(1)+b2(1)*Wf_sumr(2)+b2(2)*Wf_sumr(1);
        Wff_deltar(1) = Wf_deltar(1);Wff_deltar(2) = -a2(2)*Wff_deltar(1)+b2(1)*Wf_deltar(2)+b2(2)*Wf_deltar(1);
        end
        if i >= 3 
        Wff_sumi(i) = self_filter(Wff_sumi(i-1),Wff_sumi(i-2),Wf_sumi(i),Wf_sumi(i-1),Wf_sumi(i-2),a2,b2);
        Wff_deltai(i) = self_filter(Wff_deltai(i-1),Wff_deltai(i-2),Wf_deltai(i),Wf_deltai(i-1),Wf_deltai(i-2),a2,b2);
        Wff_sumr(i) = self_filter(Wff_sumr(i-1),Wff_sumr(i-2),Wf_sumr(i),Wf_sumr(i-1),Wf_sumr(i-2),a2,b2);
        Wff_deltar(i) = self_filter(Wff_deltar(i-1),Wff_deltar(i-2),Wf_deltar(i),Wf_deltar(i-1),Wf_deltar(i-2),a2,b2);
        end
        Wsumi = Wff_sumi(i);Wdeltai = Wff_deltai(i);Wsumr = Wff_sumr(i);Wdeltar = Wff_deltar(i);
    %% Advance time    
    stepsize = newstepsize;
    t = t + stepsize;    
end % end of main stability loop

%% Output
if output
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b> 100%% completed')
else
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
end
simulationtime=toc;
if output; fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b> Simulation completed in %5.2f seconds\n', simulationtime); end

%% Clean up
rmpath([cd '/Solvers/']);
rmpath([cd '/Models/Generators/']);
rmpath([cd '/Models/Exciters/']);
rmpath([cd '/Models/MMC-HVDC/']);
rmpath([cd '/Models/Governors/']);
rmpath([cd '/Cases/Powerflow/']);
rmpath([cd '/Cases/Dynamic/']);
rmpath([cd '/Cases/Events/']);

%% Plot
close all;
PlotMMCHVDC(i, Time, Icsi, Icsr, Vcusi, Vcusr, Vclsi, Vclsr, Is_alphasi, Is_alphasr, Is_betasi, Is_betasr, Vdsi, Vdsr, Vai, Var, Issi, Issr, Pssi, Pssr, Qssi, Qssr, VR_alphasi, VR_alphasr, VR_betasi, VR_betasr, vs_refsi, vs_refsr, Vssi, Vssr, Xc1si, Xc1sr, Xc2si, Xc2sr, vc_refsi, vc_refsr, Vcsi, Vcsr, nusi, nusr, nlsi, nlsr, Xdcsi);
PlotSG(i, Time, Angles, Speeds, Eq_tr, Ed_tr, Efd, PM, Voltages);
return;