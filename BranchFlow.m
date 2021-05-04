function MMCflow = BranchFlow(baseMVA, bus, branch, U0)

% Calculate bus admittance matrix
Ybus = makeYbus(baseMVA, bus, branch);

nb = size(Ybus,1);

Bpflow = zeros(100,4);
count = 0;

for m = 1:nb
    for n = m+1:nb
        if Ybus(m,n) ~= 0
            count = count +1;
            Ps = -real(U0(m,1)*conj(Ybus(m,n)*(U0(m,1)-U0(n,1))));
            Qs = -imag(U0(m,1)*conj(Ybus(m,n)*(U0(m,1)-U0(n,1))));
            Bpflow(count,:) = [m n Ps Qs];
        end
    end
end

Bpflow = Bpflow(1:count,:);

%% MMC is connected between bus 2 and bus 8
% output power for inverter
Pinv = Bpflow(2,3);
Qinv = Bpflow(2,4);
% output power for rectifier
Prec = -Pinv; %-real(U0(8,1)*conj(Ybus(2,8)*(U0(8,1)-U0(2,1))));
Qrec = 0;

MMCflow = [Pinv Qinv;
    Prec Qrec]; 
end