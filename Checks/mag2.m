function mag2

lecture = fopen('magSquared.txt');
vector = textscan(lecture, '%f %f %f');
fclose(lecture);

lecture1 = fopen('squaredMagnetization.txt');
vector1 = textscan(lecture1, '%f %f %f');
fclose(lecture1);

nstep = vector{1};
magnet = vector{2};
error = vector{3};

magnet1 = vector1{2};
error1 = vector1{3};

deltaT = 1./nstep;
deltaT = deltaT.^2;

e1 = errorbar(deltaT, magnet, error);
set(e1, 'Linestyle', 'none', 'Marker', '*','Linewidth', 1.5);
hold on
e2 = errorbar(deltaT, magnet1, error1);
set(e2, 'Linestyle', 'none', 'Marker', '*', 'Linewidth', 1.5);

Ofitline = polyfit(deltaT, magnet1, 1);
plot(deltaT, polyval(Ofitline,deltaT),'--r', 'Linewidth', 2);


Ofitline1 = polyfit(deltaT, magnet, 1);
plot(deltaT, polyval(Ofitline1,deltaT),'--b','Linewidth', 2);


legend('y=(-2.203\pm0.002)x+(0.224\pm0.001)', 'HMC: y=(-0.006\pm0.004)x+(0.225\pm0.001)');
set(legend,'location','southeast');
set(gca, 'fontsize', 25);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.8, 0.6, 0.68]);
xlabel('\delta^2 \tau');
ylabel('\langlem^2\rangle/V^2');
title('\fontsize{19} \lambda=1.3282, k=0.18169, D=3, L/a=4');
print('magnetization2', '-dpng');