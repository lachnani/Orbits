function [] = energyPlot(rp,ra)
%ENERGYPLOT Plots Energy, KE, and PE as function of theta
%   Takes radius of periapsis and radius of apoapsis in km and outputs plot

a = (ra+rp)/2; %semimajor axis
E = -1*398600/(2*a); %total energy
E(1:1000) = E;
theta = linspace(0,360,1000); %Theta array
e = (ra/a)-1; %eccentricity
p = a*(1-e^2);
r = p./(1+(e*cosd(theta))); %Radius array
v2 = 398600.*((2./r).*(1/a)); %Velocity squared array
KE = v2./2; %Kinetic Energy
PE = -1.*398600./r;
figure(1)
plot(theta,E)
xlabel('True Anomaly (Degrees)')
ylabel('Energy (10^6 Joules)')
title('Total Energy vs True Anomaly')
axis tight
figure(2)
plotyy(theta,KE,theta,PE)
xlabel('True Anomaly (Degrees)')
ylabel('Energy (10^6 Joules)')
title('Kinetic and Potential Energy vs True Anomaly')
legend({'KE','PE'})
axis tight

end

