function [] = PlotTemp()
global nElectrons T L W
global x y Vx Vy time Temp dt

plot(time,Temp, 'ro', 'markers',1,'MarkerFaceColor', 'b');
hold on
axis([0 time+dt 250 1050]);
xlabel('Time(s)');
ylabel('Temp (K)');
grid on
title('Temperature vs. Time');


end