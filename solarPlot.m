function [visiblePlot] = solarPlot()
%SOLARPLOT Plots the sun and visible planets
%   Plots the sun and visible planets (including Earth) relative 
%   to the vernal equinox.

%constants from JPL
% venus = [1.147,.7213];
% mars = [.3839,1.418];
% earth = [1.9323,.9834];
% jupiter = [3.3415,5.4562];
% saturn = [4.5305,10.0488];
% sun = [0,0];

%transform polar coordinates to cartesian
theta = [0 1.147 1.9323 .3839 3.3415 4.5305];
rho = [0 .7213 .9834 1.418 5.4562 10.0488];
[x,y] = pol2cart(theta,rho);

hold on 

%create axis
plot([-12 12],[0 0],'b');
plot([0 0], [-12 12],'b');

%plot visivle bodies and sun
for i = 1:6

    plot(x(i),y(i),'r*');
    plot([0 x(i)],[0 y(i)],'r');
    axis([-10 5 -12 5]);

end
    
end

