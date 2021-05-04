function [gen,exc,gov,mmc,DCVolCon,OutCurCon,CirCurCon,freq,stepsize,stoptime] = case9dyn

% case9dyn
% MatDyn dynamic data file
% 
% MatDyn
% Copyright (C) 2009 Stijn Cole
% Katholieke Universiteit Leuven
% Dept. Electrical Engineering (ESAT), Div. ELECTA
% Kasteelpark Arenberg 10
% 3001 Leuven-Heverlee, Belgium

%% General data

freq = 50;
stepsize = 0.0001;
stoptime = 5;
%% Generator data

% [genmodel excmodel govmodel H D xd xq xd_tr xq_tr Td_tr Tq_tr]
% [ genmodel    excmodel    govmodel    H       D       x       x'     ]

gen=[2  2   1   1.20   0.02   0.14    0.14   0.25    0.25    5.2 0.81;
     2  2   1   2.40   0.01   0.14    0.14   0.25    0.25    5.2 0.81;
     2  2   1   5.74   0.02   1.93    1.77    0.25    0.25    5.2 0.81;
];

%% Exciter data

% [gen Ka  Ta  Ke  Te  Kf  Tf  Aex  Bex  Ur_min  Ur_max]
exc=[1 50 0.05    -0.17   0.95    0.04    1   0.014   1.55    -1.7    1.7;
     2 50 0.05    -0.17   0.95    0.04    1   0.014   1.55    -1.7    1.7;
     3 50 0.05    -0.17   0.95    0.04    1   0.014   1.55    -1.7    1.7;
];


%% Governor data
% [gen K  T1  T2  T3  Pup  Pdown  Pmax  Pmin]
gov=[1  0   0   0   0   0   0   0   0;
     2  0   0   0   0   0   0   0   0;
     3  0   0   0   0   0   0   0   0;
];

%% MMC data: MMC-HVDC容量设计为210MVA
% mmc=[N M Vd Vsmax Ismax phi];
mmc = [12 3 140e3 70e3 2e3 0];

%% Circulating current controller
% parameters of circulating current controller
Ra = 20; % P167
alpha2 = 200; % P166
CirCurCon = [Ra alpha2];

%% Output current controller
% parameters of the output current controller of a MMC
% Refering to P 158 (3.108)
% Kp_oc = alpha_c*L/2; K0 = alpha_1*alpha_c*L alpha_1<<alpha_c
% alpha_c = 1000 rad/s and alpha_1 = 100 rad/s
alpha_c = 1000;
alpha_1 = 100;
OutCurCon = [alpha_c alpha_1];

%% DC voltage controller
% refer to P201
alpha_d = 50;
alpha_id = 25;
DCVolCon = [alpha_d alpha_id];

return;