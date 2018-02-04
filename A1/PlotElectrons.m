function [ output_args ] = PlotElectrons()
global nElectrons nPlot T L W
global x y Vx Vy

plot(x(1:nPlot), y(1:nPlot),'bo','markers',1,'MarkerFaceColor','b');
%hold on
axis([0 L 0 W]);
%xlabel('X');
%ylabel('Y');
%title('Electron Position');


color code

cc = jet(nPlot);
for i=1:nPlot
    plot(x(i), y(i), 'o','markers', 1, 'Color', cc(i,:));
end

axis([0 L 0 W]);


end
