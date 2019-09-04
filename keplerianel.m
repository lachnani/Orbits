function [rt,vt] = keplerianel(r0,v0)
%KEPLERIANEL Determine Keplerian orbits and trajectory
%   Takes r0 in km and v0 in km/s vectors at time t0, delta t for propagation,and Gm and
%   Re of central body. Outputs orbit type, Keplerian elements, ballistic
%   analysis, and propagates r and v.

%Part 1
%Compute eccentricity using problem 2.4
r0mag = norm(r0);%Magnitude of r0 in km
v0mag = norm(v0);%Magnitude of v0 in km/s
Gm = input('Enter the Gm of the central body in km^3/s^2: ');
evec = ((v0mag^2 - (Gm/r0mag)).*r0 - (dot(r0,v0)).*v0)./Gm;%Eccentricity vector by problem 2.4
e = norm(evec);%Magnitude of e
%Determine orbit type from e
if e <= 1e-3 
    disp('The orbit is circular')
elseif e < 1+(1e-3) && e > 1-(1e-3)
    disp('The orbit is parabolic')
elseif e > 1+(1e-3)
    disp('The orbit is hyperbolic')
else
    disp('The orbit is elliptical')
end
disp('-----------------------------------')

%Part 2
h = cross(r0,v0); %Angular momentum by definition in km^2/s
energy = ((v0mag^2)/2)-(Gm/r0mag); %Energy by 2.20
p = dot(h,h)/Gm; %Semilatus rectum by definition in km
a = -1*Gm/(2*energy); %Semimajor axis by 2.20 in km
%Cartesian to Kepler
i = acosd(h(3)/norm(h)); %inclination by 4.40 in degrees
K = [0,0,1]; %K by definition
N = cross(K,h); %N by definition
Nhat = N./norm(N); %Nhat by 4.41
omega = acosd(Nhat(1)); %RAAN by 4.42 in degrees
if Nhat(2) < 0 %Quadrant check
    omega = 360 - omega;
end
argument = acosd(dot(evec,Nhat)/e); %Argument of Periapsis by 4.43 in degrees
if evec(3) < 0 %Quadrant check
    argument = 360 - argument;
end
Theta = acosd((p-r0mag)/(e*r0mag)); %True Anomaly by conic equation in degrees
ir = r0./norm(r0); %Radius direction
if sign(dot(v0,ir)) == -1 %Check approach or departure
    Theta = 360 - Theta; %Adjust Theta accordingly
end
%Conic Section cases
if e < 1-(1e-3) %elliptic case
    period = 2*pi*sqrt((a^3)/Gm); %Period by 2.34 in seconds
    rp = a*(1-e); %Radius of periapsis by exercise 2.10 in km
    ra = a*(1+e); %Radius of apoapsis by exercise 2.10 in km
    %To find True and Eccentric Anomaly
    ip = evec./e; %Periapsis direction by definition of eccentricity vector
    E = acosd((1 - (r0mag/a))/e); %Eccentric Anomaly by 4.12 in degrees
    %Quadrant check using True Anomaly (must be in same hemisphere)
    if Theta > 180
        E = 360 - E;
    end
    n = 360/period; %Mean motion in deg/sec
    M = E - e*sind(E); %Mean anomaly in degrees
    tp = M/n; %Time since periapsis in seconds
    fprintf('eccentricity          = %6.3f\n',e)
    fprintf('semilatus rectum      = %6.3f km\n',p)
    fprintf('semimajor axis        = %6.3f km\n',a)
    fprintf('period                = %6.3f sec\n',period)
    fprintf('radius of periapsis   = %6.3f km\n',rp)
    fprintf('radius of apoapsis    = %6.3f km\n',ra)
    fprintf('inclination           = %6.3f degrees\n',i)
    fprintf('RAAN                  = %6.3f degrees\n',omega)
    fprintf('argument of periapsis = %6.3f degrees\n',argument)
    fprintf('true anomaly          = %6.3f degrees\n',Theta)
    fprintf('eccemtric anomaly     = %6.3f degrees\n',E)
    fprintf('mean anomaly          = %6.3f degrees\n',M)
    fprintf('time since periapsis  = %6.3f sec\n',tp)
    disp('-----------------------------------')
elseif e < 1+(1e-3) && e > 1-(1e-3) %Parabolic case
    rp = p/2; %Radius of Periapsis by Conic Equation in km
    fprintf('semilatus rectum      = %6.3f km\n',p)
    fprintf('radius of periapsis   = %6.3f km\n',rp)
    fprintf('inclination           = %6.3f degrees\n',i)
    fprintf('RAAN                  = %6.3f degrees\n',omega)
    fprintf('argument of periapsis = %6.3f degrees\n',argument)
    fprintf('true anomaly          = %6.3f degrees\n',Theta)
    disp('-----------------------------------')
elseif e> 1+(1e-3) %Hyperbolic case
    rp = a*(1-e); %Radius of periapsis by exercise 3.26 in km
    vinf = sqrt(-1*Gm/a); %Vinfinity by 3.25 on km/s
    fprintf('eccentricity          = %6.3f\n',e)
    fprintf('semilatus rectum      = %6.3f km\n',p)
    fprintf('semimajor axis        = %6.3f km\n',a)
    fprintf('V infinity            = %6.3f km/s\n',vinf)
    fprintf('radius of periapsis   = %6.3f km\n',rp)
    fprintf('inclination           = %6.3f degrees\n',i)
    fprintf('RAAN                  = %6.3f degrees\n',omega)
    fprintf('argument of periapsis = %6.3f degrees\n',argument)
    fprintf('true anomaly          = %6.3f degrees\n',Theta)
    disp('-----------------------------------')
end

%Part 3
Re = input('Enter the radius of the central body in km: ');
if e >= 1 && Theta < 180 %Spacecraft has passed and will not return
    disp('The orbit will not impact the central body')    
elseif rp <= Re %Spacecraft will impact
    disp('The orbit is ballistic and the spacecraft will impact the central body')
else
    disp('The orbit will not impact the central body')
end
disp('-----------------------------------')

%Part 4
%Define Alpha0
if e == 1
    alpha0 = 0;
else
    alpha0 = (2/r0mag) - (dot(v0,v0)/Gm);
end
%Will solve for x via Newton's method. Stumpff functions defined below
x1 = 0;%Dummy variable
x2 = 1;%Initial guess
dt = input('Enter the delta t in seconds: ');
K = dt*sqrt(Gm);%LHS of Universal Kepler's Equation
while abs(x2-x1) >= .001 %Tolerance
    x1 = x2;
    Fx = ((dot(r0,v0)/sqrt(Gm))*(x1^2)*stumpffc(alpha0*(x1^2)))+...
        ((1-r0mag*alpha0)*(x1^3)*stumpffs(alpha0*(x1^2)))+(r0mag*x1)-K; %F(x)
    Fxp = ((dot(r0,v0)/sqrt(Gm))*(x1-(alpha0*(x1^3)*stumpffs(alpha0*(x1^2)))))+...
        ((1-r0mag*alpha0)*(x1^2)*stumpffc(alpha0*(x1^2)))+r0mag; %F'(x)
    x2 = x1 - (Fx/Fxp);
end
sf = stumpffs(alpha0*(x2^2));%final stumpff s value
cf = stumpffc(alpha0*(x2^2));%final stumpff c value
%Evaluate propagated vector from 4.35 and 4.36
rt = (1-((x2^2)*cf/r0mag)).*r0 + (dt - ((x2^3)*sf/sqrt(Gm))).*v0;%r(t) in km 
vt = ((sqrt(Gm)/(r0mag*norm(rt))).*((alpha0*(x2^3)*sf) - x2).*r0) +...
    (1 - (x2^2)*cf/norm(rt)).*v0;%v(t) in km/s
end
%Stumpff functions
function [sx] = stumpffs(x)
%STUMPFFS Computes S Stumpff function
if x > 10e-7
    sx = (sqrt(x)-sin(sqrt(x)))/(x^(3/2));
elseif x<-10e-7
    sx = (sinh(sqrt(-x))-sqrt(-x))/((-x)^(3/2));
else
    sx = 1/6;
end
end
function [cx] = stumpffc(x)
%STUMPFFC Computes C Stumpff function
if x > 10e-7
    cx = (1-cos(sqrt(x)))/x;
elseif x<-10e-7
    cx = (cosh(sqrt(-x))-1)/(-x);
else
    cx = 1/2;
end
end

