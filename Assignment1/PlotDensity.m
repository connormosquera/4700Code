function [] = PlotDensity()
global nElectrons T L W MarkerSize
global x y Vx Vy

n=hist3([x',y'],[50 50]);

pcolor(n');
colorbar;


end