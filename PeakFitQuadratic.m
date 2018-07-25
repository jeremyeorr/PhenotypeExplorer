function [x_peak,y_peak] = PeakFitQuadratic(x,y);

%x = [aD bD cD];
%y = [AD BD CD];

a = [y(3)*x(2)-y(3)*x(1)-x(3)*y(2)+x(3)*y(1)-y(1)*x(2)+x(1)*y(2)]/[(x(2)-x(1))*(-x(2)*x(3)+x(2)*x(1)+x(3)^2-x(1)*x(3))];
b = [y(2)-a*x(2)^2-y(1)+a*x(1)^2]/[x(2)-x(1)];
c = y(1)-a*x(1)^2-b*x(1);

%x1 = 2:0.001:2.02;
%y1 = a*x1.^2+b*x1+c;

x_peak = -b/(2*a);
y_peak = a*x_peak.^2+b*x_peak+c;

%figure(200), plot(x,y,'.',x1,y1,x_peak,y_peak,'o');
