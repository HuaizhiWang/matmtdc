function Pcir0 = Loadcir(casefile_dyn)

[gen,exc,gov,mmc,DCVolCon,OutCurCon,CirCurCon,freq,stepsize,stoptime] = feval(casefile_dyn);


Pcir0 = CirCurCon;

end