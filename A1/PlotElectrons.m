function [ output_args ] = PlotElectrons(cc)
global nElectrons T L W MarkerSize
global x y Vx Vy

plot(x, y, 'bo', 'markers',MarkerSize,'MarkerFaceColor', 'b');
%hold on
axis([0 L 0 W]);
%xlabel('X');
%ylabel('Y');
%title('Electron Position');


% color code

% for i=1:nElectrons
%     plot(x(i), y(i), 'o','markers', 1, 'Color', cc(i,:));
% end
% 
% axis([0 L 0 W]);


end
