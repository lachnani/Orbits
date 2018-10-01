function [periods,periodd,velocity] = orbitcalc(a)
%ORBITCALC Calculates Period and Velocity from radius and mu value
%   Assumes a circular orbit around the sun. Takes in the semi major axis
%   and interprets it as the radius. Outputs the period in seconds and
%   Sidereal days, and the velocity in km/s. 

%check input
if ~isscalar(a)
    error('Input must be a scalar');
end
if a <= 0
    error('Input must be positive');
end

%set solar mu
mu = 1.327e+11;
%calculate velocity
velocity = sqrt(mu/a);
%calculate periods
periods = (2*pi*a)/velocity;
periodd = periods/86164;


end

