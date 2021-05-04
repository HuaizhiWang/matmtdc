function dDcc0 = DvotCon(Xdcc0, Pdcc0, Vdcc0, Vmmc0, Addvai)
% Vmmc0 = [Va_alpha; Va_beta; is0; Ps; Qs; id0]
Va_alpha = Vmmc0(1,1);
Va_beta = Vmmc0(2,1);
Va = sqrt(Va_alpha^2+Va_beta^2);
% Vai = Addvai(2,1);
Va = Va + 1*Addvai(1,1);
% Xdcc0 = xdc0
xdc = Xdcc0;
% Pdcc0 = [alpha_d alpha_id]
alpha_id = Pdcc0(1,2);
% Vdcc0 = Is_alphar
Is_alphar = Vdcc0;

dxdc = alpha_id*(3*Va*Is_alphar/2-xdc);

dDcc0 = dxdc;

end