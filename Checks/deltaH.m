function deltaH

lecture = fopen('DH.txt');
vector = textscan(lecture, '%f %f %f');
fclose(lecture);

nstep = vector{1};
deltaH = vector{2};
error = vector{3};

deltaT = 1./nstep;
deltaT = deltaT.^2;

e = errorbar(deltaT, deltaH, error);
set(e, 'Linestyle', 'none', 'Marker', '*','Linewidth', 1.5);
hold on
Ofitline = polyfit(deltaT, deltaH, 1);
plot(deltaT, polyval(Ofitline,deltaT),'--b', 'Linewidth', 2);

set(gca, 'fontsize', 25);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.8, 0.6, 0.68]);
xlabel('\delta^2\tau');
ylabel('|\DeltaH|');
title('\fontsize{19}\lambda=1.3282, k=0.18169, D=3, L/a=4');
legend('data', 'fit: y=(17.3\pm0.3)x-(0.008\pm0.004)');
set(legend,'location','southeast');
print('deltaH', '-dpng');