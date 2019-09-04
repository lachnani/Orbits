function [vinj] = inject2planet(rpark,rplanet)
% INJECT2PLANET Calculates injection delta V to planet
%       Input parking orbit radius and mean planet-sun distance
re = 1.496e8; %Radius of Earth Orbit
mus = 1.327e11; %Gm of sun
at = (re+rplanet)/2; %Transfer orbit to Planet Semimajor axis
vp = sqrt(mus*((2/re)-(1/at))); %V+ 
ve = 29.785; %Velocity of Earth Relative to Sun
vinf = abs(ve-vp); %V infinity
mue = 398600; %Gm of earth
vpe = sqrt((vinf^2)+(2*mue/rpark)); %Perigee velocity
vinj = vpe - sqrt(mue/rpark); %Delta V injection
fprintf('atransfer = %1.3f km\n',at)
fprintf('v+ = %1.3f km/sec\n',vp)
fprintf('vinfinity/earth = %1.3f km/sec\n',vinf)
fprintf('vperigee/earth = %1.3f km/sec\n',vpe)
fprintf('DeltaV Injection = %1.3f km/sec\n',vinj)
disp('-----------------------------')
e = 1 + (rpark*(vinf^2)/mue);
thetainf = acosd(-1/e);
do2 = thetainf-90;
b = mue*sqrt((e^2)-1)/(vinf^2);
fprintf('eccentricity = %1.3f\n',e)
fprintf('theta infinity = %1.3f degrees\n',thetainf)
fprintf('Delta/2 = %1.3f degrees\n',do2)
fprintf('B = %1.3f km\n',b)
end

