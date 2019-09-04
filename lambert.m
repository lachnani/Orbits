function [time, p] = lambert(r1,r2,c,a)
%LAMBERT Takes radius vectors r1 and r2 from the Sun to P and Q
%respectively with the assumption r1<r2. Also takes chord length between P
%and Q and semimajor axis. Outputs time and semilatus rectum.

mu = 132712000000;
r1mag = norm(r1); r2mag = norm(r2);
nu = acos(dot(r1,r2)/(r1mag*r2mag)); %heliocentric true anomaly
tau = 2*pi*sqrt((a^3)/mu); %orbit period
otype = input('Enter Orbit Type (e, p, or h)');
s = (r1mag + r2mag + c)/2; %semiperimeter of the triangle
%Minimums
if otype == 'e'
    cm = sqrt(r1mag^2 + r2mag^2 - 2*r1mag*r2mag*cos(nu)); %chord between points
    sm = (r1mag + r2mag + cm)/2; %semiperimeter of the triangle
    am = sm/2; %minimum semimajor axis of transfer ellipse
    pm = (-2/cm)*(sm-r1mag)*(sm-r2mag); %parameter of miminum energy ellipse
    Em = -mu(2*am); %minimum energy of transfer ellipse
    em = sqrt(1 - 2*pm/sm); %eccentricity of minimum energy ellipse
elseif otype == 'h'
    am = 0; %minimum semimajor axis of transfer ellipse
end
% Solving for semilatus rectum
alpha = 2*asin(sqrt(s/(2*a)));
beta = 2*asin(sqrt((s-c)/(2*a)));
gamma = 2*asinh(sqrt(s/(-2*a)));
delta = 2*asinh(sqrt((s-c)/(-2*a)));
if otype == 'e'
    p1 = (4*a/c^2)*(s-r1mag)*(s-r2mag)*(sin((alpha+beta)/2))^2; %parameter 1
    p2 = (4*a/c^2)*(s-r1mag)*(s-r2mag)*(sin((alpha-beta)/2))^2; %parameter 1
elseif otype == 'h'
    p1 = (4*a/c^2)*(s-r1mag)*(s-r2mag)*(sinh((gamma+delta)/2))^2; %parameter 1
    p2 = (4*a/c^2)*(s-r1mag)*(s-r2mag)*(sinh((gamma-delta)/2))^2; %parameter 1
elseif otype == 'p'
    p1 = (4/c^2)*(s-r1mag)*(s-r2mag)*(sqrt(s/2)+sqrt((s-c)/2)); %parameter 1
    p2 = (4/c^2)*(s-r1mag)*(s-r2mag)*(sqrt(s/2)-sqrt((s-c)/2)); %parameter 2
end
% Solving for time (Langrange's  solutions)
%The set of permissible F={F|r1+r2 = constant} is an ellipse with foci at P
%and Q and 2aF = r1+r2
%The set {F*}={F*|PF*+QF*=4a-(r1+r2) = constant} is an ellipse with foci at
%P and Q and 2aF*=4a-(r1+r2)
if otype == 'e'
    fstar = input('Enter branch of the vacant focus on hyperbola of {F*} (lower or upper):');
    if fstar == 'lower'
        time = tau * [(alpha - sin(alpha)) - (beta - sin(beta))]/(2*pi);
        p = p1;
        if nu > pi
            time = tau - time;
        end
    elseif fstar == 'upper'
        time = tau - (tau * ((alpha - sin(alpha)) + (beta - sin(beta)))/(2*pi));
        p = p2;
        if nu > pi
            time = tau - time;
        end
    else
        error('Invalid branch')
    end
elseif otype == 'h'
    fstar = input('Enter branch of the vacant focus on hyperbola of {F*} (lower or upper):');
    if fstar == 'upper'
        time = sqrt((-a^3)/mu)*((sinh(gamma)-gamma) - (sinh(delta)-delta));
        p = p1;
    elseif fstar == 'lower'
        time = sqrt((-a^3)/mu)*((sinh(gamma)-gamma) - (sinh(delta)-delta));
        p = p2;
    else
        error('Invalid branch')
    end
elseif otype == 'p'
    fstar = input('Enter branch of the vacant focus on hyperbola of {F*} (lower or upper):');
    if fstar == 'upper'
        time = sqrt(2/mu)*(s^(3/2)-(s-c)^(3/2))/3;
        p = p1;
    elseif fstar == 'lower'
        time = sqrt(2/mu)*(s^(3/2)+(s-c)^(3/2))/3;
        p = p2;
    else
        error('Invalid branch')
    end
end    
end

