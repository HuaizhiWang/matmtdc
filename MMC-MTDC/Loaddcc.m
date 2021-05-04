function Pdcc0 = Loaddcc(casefile_dyn)
[gen,exc,gov,mmc,DCVolCon,OutCurCon,CirCurCon,freq,stepsize,stoptime] = feval(casefile_dyn);
% DCVolCon = [alpha_d alpha_id]
Pdcc0 = DCVolCon;
end