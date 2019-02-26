function binderCumulant

%L=4
lecture4 = fopen('binderCumulantL4.txt');
vector4 = textscan(lecture4, '%f %f %f');
fclose(lecture4);

kappa4 = vector4{1};
binder4 = vector4{2};
error4 = vector4{3};
%--------------------------------------------------

%L=6
lecture6 = fopen('binderCumulantL6.txt');
vector6 = textscan(lecture6, '%f %f %f');
fclose(lecture6);

kappa6 = vector6{1};
binder6 = vector6{2};
error6 = vector6{3};
%--------------------------------------------------

%L=8
lecture8 = fopen('binderCumulantL8.txt');
vector8 = textscan(lecture8, '%f %f %f');
fclose(lecture8);

kappa8 = vector8{1};
binder8 = vector8{2};
error8 = vector8{3};
%--------------------------------------------------

%L=10
lecture10 = fopen('binderCumulantL10.txt');
vector10 = textscan(lecture10, '%f %f %f');
fclose(lecture10);

kappa10 = vector10{1};
binder10 = vector10{2};
error10 = vector10{3};
%--------------------------------------------------

%L=12
lecture12 = fopen('binderCumulantL12.txt');
vector12 = textscan(lecture12, '%f %f %f');
fclose(lecture12);

kappa12 = vector12{1};
binder12 = vector12{2};
error12 = vector12{3};
%--------------------------------------------------

%L=14
lecture14 = fopen('binderCumulantL14.txt');
vector14 = textscan(lecture14, '%f %f %f');
fclose(lecture14);

kappa14 = vector14{1};
binder14 = vector14{2};
error14 = vector14{3};
%--------------------------------------------------

%L=16
lecture16 = fopen('binderCumulantL16.txt');
vector16 = textscan(lecture16, '%f %f %f');
fclose(lecture16);

kappa16 = vector16{1};
binder16 = vector16{2};
error16 = vector16{3};
%--------------------------------------------------

errorbar(kappa4, binder4, error4, 'Linewidth', 2,'Marker', '.', 'MarkerSize', 20, 'Linestyle', 'none');
hold on
errorbar(kappa6, binder6, error6,'Linewidth', 2,'Marker', '.', 'MarkerSize', 20, 'Linestyle', 'none');
errorbar(kappa8, binder8, error8,'Linewidth', 2,'Marker', '.', 'MarkerSize', 20, 'Linestyle', 'none');
errorbar(kappa10, binder10, error10,'Linewidth', 2,'Marker', '.', 'MarkerSize', 20, 'Linestyle', 'none');
errorbar(kappa12, binder12, error12,'Linewidth', 2,'Marker', '.', 'MarkerSize', 20, 'Linestyle', 'none');
errorbar(kappa14, binder14, error14,'Linewidth', 2,'Marker', '.', 'MarkerSize', 20, 'Linestyle', 'none');
errorbar(kappa16, binder16, error16,'Linewidth', 2,'Marker', '.', 'MarkerSize', 20, 'Linestyle', 'none');

line([0.18644 0.18644], [0.9 3], 'linestyle', '--','Linewidth', 2);

axis([0.15 0.23 0.9 3]);

set(gca, 'fontsize', 25);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 1, 0.75, 0.68]);
xlabel('k');
ylabel('U=\langlem^4\rangle/(\langlem^2\rangle)^2');
title('\fontsize{19} Binder Cumulant, \lambda=1.145, D=3');
legend('L/a = 4', 'L/a = 6', 'L/a = 8', 'L/a = 10', 'L/a = 12', 'L/a = 14', 'L/a = 16')
print('binder', '-dpng');

end
