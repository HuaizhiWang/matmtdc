function PlotMMCHVDC(i, Time, Icsi, Icsr, Vcusi, Vcusr, Vclsi, Vclsr, Is_alphasi, Is_alphasr, Is_betasi, Is_betasr, Vdsi, Vdsr, Vai, Var, Issi, Issr, Pssi, Pssr, Qssi, Qssr, VR_alphasi, VR_alphasr, VR_betasi, VR_betasr, vs_refsi, vs_refsr, Vssi, Vssr, Xc1si, Xc1sr, Xc2si, Xc2sr, vc_refsi, vc_refsr, Vcsi, Vcsr, nusi, nusr, nlsi, nlsr, Xdcsi)

Time = Time(1:i,:);
Icsi = Icsi(1:i,:); Icsr = Icsr(1:i,:);
Vcusi = Vcusi(1:i,:); Vcusr = Vcusr(1:i,:);
Vclsi = Vclsi(1:i,:); Vclsr = Vclsr(1:i,:);
Is_alphasi = Is_alphasi(1:i,:); Is_alphasr = Is_alphasr(1:i,:);
Is_betasi = Is_betasi(1:i,:); Is_betasr = Is_betasr(1:i,:);
Vdsi = Vdsi(1:i,:); Vdsr = Vdsr(1:i,:);
Vai = Vai(1:i,:); Var = Var(1:i,:);
Issi = Issi(1:i,:); Issr = Issr(1:i,:);
Pssi = Pssi(1:i,:); Pssr = Pssr(1:i,:);
Qssi = Qssi(1:i,:); Qssr = Qssr(1:i,:);
VR_alphasi = VR_alphasi(1:i,:); VR_alphasr = VR_alphasr(1:i,:);
VR_betasi = VR_betasi(1:i,:); VR_betasr = VR_betasr(1:i,:);
vs_refsi = vs_refsi(1:i,:); vs_refsr = vs_refsr(1:i,:);
Vssi = Vssi(1:i,:); Vssr = Vssr(1:i,:);
Xc1si = Xc1si(1:i,:); Xc1sr = Xc1sr(1:i,:);
Xc2si = Xc2si(1:i,:); Xc2sr = Xc2sr(1:i,:);
% vc_refsi = vc_refsi(1:i,:); vc_refsr = vc_refsr(1:i,:);
Vcsi = Vcsi(1:i,:); Vcsr = Vcsr(1:i,:);
nusi = nusi(1:i,:); nusr = nusr(1:i,:);
nlsi = nlsi(1:i,:); nlsr = nlsr(1:i,:);
Xdcsi = Xdcsi(1:i,:);

figure
xlabel('Time (s)');
ylabel('$i_{\mathrm{ci}}$','interpreter','latex');
hold on;
plot(Time, Icsi, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$v_{\mathrm{cui}}$','interpreter','latex');
hold on;
plot(Time, Vcusi, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$v_{\mathrm{cli}}$','interpreter','latex');
hold on;
plot(Time, Vclsi, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$I_{\mathrm{s}\alpha\mathrm{i}}$','interpreter','latex');
hold on;
plot(Time, Is_alphasi, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$I_{\mathrm{s}\beta\mathrm{i}}$','interpreter','latex');
hold on;
plot(Time, Is_betasi, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$V_{\mathrm{di}}$','interpreter','latex');
hold on;
plot(Time, Vdsi, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$V_{\mathrm{ai}}$','interpreter','latex');
hold on;
plot(Time, Vai, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$i_{\mathrm{si}}$','interpreter','latex');
hold on;
plot(Time, Issi, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$P_{\mathrm{si}}$','interpreter','latex');
hold on;
plot(Time, Pssi, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$Q_{\mathrm{si}}$','interpreter','latex');
hold on;
plot(Time, Qssi, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$V_{\mathrm{R}\alpha\mathrm{i}}$','interpreter','latex');
hold on;
plot(Time, VR_alphasi, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$V_{\mathrm{R}\beta\mathrm{i}}$','interpreter','latex');
hold on;
plot(Time, VR_betasi, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$v_{\mathrm{si}}$','interpreter','latex');
hold on;
plot(Time, Vssi, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$v_{\mathrm{refi}}$','interpreter','latex');
hold on;
plot(Time, vs_refsi, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$x_{\mathrm{c1i}}$','interpreter','latex');
hold on;
plot(Time, Xc1si, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$x_{\mathrm{c2i}}$','interpreter','latex');
hold on;
plot(Time, Xc2si, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$v_{\mathrm{ci}}$','interpreter','latex');
hold on;
plot(Time, Vcsi, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$n_{\mathrm{ui}}$','interpreter','latex');
hold on;
plot(Time, nusi, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$n_{\mathrm{li}}$','interpreter','latex');
hold on;
plot(Time, nlsi, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$x_{\mathrm{dci}}$','interpreter','latex');
hold on;
plot(Time, Xdcsi, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);






figure
xlabel('Time (s)');
ylabel('$i_{\mathrm{cr}}$','interpreter','latex');
hold on;
plot(Time, Icsr, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$v_{\mathrm{cur}}$','interpreter','latex');
hold on;
plot(Time, Vcusr, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$v_{\mathrm{clr}}$','interpreter','latex');
hold on;
plot(Time, Vclsr, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$I_{\mathrm{s}\alpha\mathrm{r}}$','interpreter','latex');
hold on;
plot(Time, Is_alphasr, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$I_{\mathrm{s}\beta\mathrm{r}}$','interpreter','latex');
hold on;
plot(Time, Is_betasr, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$V_{\mathrm{dr}}$','interpreter','latex');
hold on;
plot(Time, Vdsr, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$V_{\mathrm{ar}}$','interpreter','latex');
hold on;
plot(Time, Var, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$i_{\mathrm{sr}}$','interpreter','latex');
hold on;
plot(Time, Issr, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$P_{\mathrm{sr}}$','interpreter','latex');
hold on;
plot(Time, Pssr, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$Q_{\mathrm{sr}}$','interpreter','latex');
hold on;
plot(Time, Qssr, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$V_{\mathrm{R}\alpha\mathrm{r}}$','interpreter','latex');
hold on;
plot(Time, VR_alphasr, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$V_{\mathrm{R}\beta\mathrm{r}}$','interpreter','latex');
hold on;
plot(Time, VR_betasr, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$v_{\mathrm{sr}}$','interpreter','latex');
hold on;
plot(Time, Vssr, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$v_{\mathrm{refr}}$','interpreter','latex');
hold on;
plot(Time, vs_refsr, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$x_{\mathrm{c1r}}$','interpreter','latex');
hold on;
plot(Time, Xc1sr, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$x_{\mathrm{c2r}}$','interpreter','latex');
hold on;
plot(Time, Xc2sr, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$v_{\mathrm{cr}}$','interpreter','latex');
hold on;
plot(Time, Vcsr, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$n_{\mathrm{ur}}$','interpreter','latex');
hold on;
plot(Time, nusr, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time (s)');
ylabel('$n_{\mathrm{lr}}$','interpreter','latex');
hold on;
plot(Time, nlsr, 'k-', 'LineWidth', 2);
hold on;
% axis([-0.01 Times(end) -0.1 0]);
set(gca,'FontName','Times New Roman','FontSize',24);

end