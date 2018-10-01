function r = conicRadius( rp, ecc, theta )
%CONICRADIUS Answer the radius of a conic section
%with periapsis radius rp, eccentricity ecc, 
%at true anomaly theta (degrees).

r = rp*( 1 + ecc) ./ ( 1 + ecc*cosd(theta) ); % Polar equation
end

