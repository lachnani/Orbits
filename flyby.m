function [] = flyby(palt,Gmp,rplanet,flyrad,sundark)
%FLYBY Computes interplanetary flyby parameters
%   Inputs:
%   parking altitude in km
%   Gm of target planet in km^3/s^2
%   radius of planet orbit in km
%   radius of flyby closest approach in km
%   0 if dark passage. 1 if sunlit passage
rpark = palt+6378.14; %Radius of parking orbit in km
re = 1.496e8; %Radius of Earth Orbit in  km
mus = 1.327e11; %Gm of sun in km^3/s^2
at = (re+rplanet)/2; %Transfer orbit to Planet Semimajor axis in km
vp = sqrt(mus*((2/re)-(1/at))); %V+ in km/s
ve = 29.785; %Velocity of Earth Relative to Sun in km/s
vinfp = abs(ve-vp); %V infinity in km/s
mue = 398600; %Gm of earth in km^3/s^2
vpe = sqrt((vinfp^2)+(2*mue/rpark)); %Perigee velocity in km/s
vinj = vpe - sqrt(mue/rpark); %Delta V injection in km/s
% %Hyperbolic exit
% e = 1 + (rpark*(vinfp^2)/mue); %eccentricity
% thetainf = acosd(-1/e); %theta infinity from earth in degrees
% do2 = thetainf-90; %delta/2 in degrees
% b = mue*sqrt((e^2)-1)/(vinfp^2); %impact parameter in km
% fprintf('eccentricity = %1.3f\n',e)
% fprintf('theta infinity = %1.3f degrees\n',thetainf)
% fprintf('Delta/2 = %1.3f degrees\n',do2)
% fprintf('B = %1.3f km\n',b)
%Hyperbolic approach
vplanet = sqrt(mus/rplanet); %Velocity of target in km/s
vmin = re*vp/rplanet; %Vminus in km/s
vinfm = abs(vmin-vplanet); %Vinfinityminus in km/s
e2 = 1 + (flyrad*(vinfm^2)/Gmp); %Shape of approach hyperbola
do3 = 2*asind(1/e2);%Delta deflection in degrees
vp2 = sqrt((vinfm^2)+(vplanet^2)-(2*vinfm*vplanet*cosd(180-do3))); %Vplus at exit
%Post Venus Orbit
if sundark == 0 %Check side of passage
    sgnb = -1;
elseif sundark == 1
    sgnb = 1;
end
beta = sgnb*acosd(((vplanet^2)+(vp2^2)-(vinfm^2))/(2*vplanet*vp2)); %Flight path angle in degrees
X0 = rplanet*(vp2^2)/mus; %Calculation parameter
e3 = sqrt(((X0-1)^2)*(cosd(beta)^2)+(sind(beta)^2)); %Eccentricity of sun orbit
theta = atand(X0*sind(beta)*cosd(beta)/(X0*(cosd(beta)^2)-1)); %True anomaly of sun orbit in degrees
if theta <=0 %Make theta positive
    theta = theta + 360;
end
if theta >= 180 %Theta quadrant check
    if beta >= 0 
        theta = theta-180;
    end
elseif theta <= 180
    if beta <= 0
        theta = theta+180;
    end
end
h = rplanet*vp2*cosd(beta); %angular momentum in km^2/s
rp = ((h^2)/mus)/(1+e3); %radius of perihelion in km
a = (h^2)/(mus*(1-e3^2)); %semimajor axis in km
tau = 2*pi*sqrt((a^3)/mus); %period in km
disp('-----------------------------')
fprintf('deltaV Injection = %6.3f km/s\n',vinj)
fprintf('radius           = %6.3f km\n',rplanet)
fprintf('velocity         = %6.3f km/s\n',vp2)
fprintf('beta             = %6.3f degrees\n',beta)
fprintf('eccentricity     = %6.3f\n',e3)
fprintf('theta            = %6.3f degrees\n',theta)
fprintf('angular momentum = %6.3f km^2/s\n',h)
fprintf('r of periapsis   = %6.3f km\n',rp)
fprintf('semimajor axis   = %6.3f km\n',a)
fprintf('period           = %6.3f sec\n',tau)
end

