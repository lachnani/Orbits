function [ellipse] = makeellipse(a,e)
%MAKEELLIPSE Draws an ellipse or circle
%   Takes in the semimajor axis (or radius in the case of a circle) and the
%   eccentricity. 

theta = linspace(0,360,1000);

if e == 0
  r = a;
else
  r = (a*(1-e^2))./(1+(e.*cosd(theta))); 
end

x = r.*cosd(theta);
y = r.*sind(theta);

plot(x,y);


end

