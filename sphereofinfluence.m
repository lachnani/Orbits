function [rsoi] = sphereofinfluence(m1,m2,d)
%SPHEREOFINFLUENCE Calculates radius of sphere of influence from body m1 to
%to body m2, centered at m1. Takes the masses and distance between bodies
rsoi = d/(((m1/m2)^(-2/5))+1);
end

