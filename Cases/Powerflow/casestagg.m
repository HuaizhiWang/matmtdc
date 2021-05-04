function [baseMVA, bus, gen, branch, area, gencost] = casestagg

%%-----  Power Flow Data  -----%%
%% system MVA base
baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
bus = [
	1	3	0	0	0	0	1	1	0	345	1	1.1	0.9;
	2	1	20	10	0	0	1	1	0	345	1	1.1	0.9;
	3	1	45	15	0	0	1	1	0	345	1	1.1	0.9;
	4	1	40	5	0	0	1	1	0	345	1	1.1	0.9;
	5	1	60	10	0	0	1	1	0	345	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin
gen = [
	1	0	0	300	-300	1.06	100	1	250	10;
	2	40	30	300	-300	1.06	100	1	300	10;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status
branch = [
	1	2	0.02	0.06	0.030*2	250	250	250	0	0	1;
	1	3	0.08	0.24	0.025*2	250	250	250	0	0	1;
	2	3	0.06	0.18	0.020*2	250	250	250	0	0	1;
	2	4	0.06	0.18	0.020*2	250	250	250	0	0	1;
	2	5	0.04	0.12	0.015*2	250	250	250	0	0	1;
	3	4	0.01	0.03	0.010*2	250	250	250	0	0	1;
	4	5	0.08	0.24	0.025*2	250	250	250	0	0	1;
];

area=[];
gencost=[];

return;
