function PlotSG(i, Time, Angles, Speeds, Eq_tr, Ed_tr, Efd, PM, Voltages)
% Save only the first i elements
Time = Time(1:i,:);
Angles = Angles(1:i,:);
Speeds = Speeds(1:i,:);
Eq_tr = Eq_tr(1:i,:);
Ed_tr = Ed_tr(1:i,:);
Efd = Efd(1:i,:);
PM = PM(1:i,:);
Voltages = Voltages(1:i,:);

figure
xlabel('Time [s]');
ylabel('$\delta$ [deg]','interpreter','latex');
hold on;
plot(Time, Angles, 'k-', 'LineWidth', 2);
hold on;
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time [s]');
ylabel('$\omega_{\mathrm{r}}$ [p.u.]','interpreter','latex');
hold on;
plot(Time, Speeds, 'k-', 'LineWidth', 2);
hold on;
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time [s]');
ylabel('$E_{q}$ [p.u.]','interpreter','latex');
hold on;
plot(Time, Eq_tr, 'k-', 'LineWidth', 2);
hold on;
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time [s]');
ylabel('$E_{d}$ [p.u.]','interpreter','latex');
hold on;
plot(Time, Ed_tr, 'k-', 'LineWidth', 2);
hold on;
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time [s]');
ylabel('$V_{\mathrm{bus}}$ [p.u.]','interpreter','latex');
hold on;
plot(Time, abs(Voltages), 'k-', 'LineWidth', 2);
hold on;
set(gca,'FontName','Times New Roman','FontSize',24);

figure
xlabel('Time [s]');
ylabel('$E_{\mathrm{f}}$ [p.u.]','interpreter','latex');
hold on;
plot(Time, Efd, 'k-', 'LineWidth', 2);
hold on;
set(gca,'FontName','Times New Roman','FontSize',24);

end