function [ output_args ] = PlotTemp(Limits)
global nElectrons T L W MarkerSize
global x y Vx Vy time Temp dt

plot(time,Temp, 'ro', 'markers',MarkerSize,'MarkerFaceColor', 'b');
hold on
axis([0 time+dt 400 500]);
xlabel('Time(s)');
ylabel('Temp (K)');
grid on
title('Temperature vs. Time');


end