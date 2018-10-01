% FUNCTION plotConic
% Plots a conic section given inputs.
% INPUTS
% rp: (float) periapsis radius
% ecc: (float) eccentricity
% style: (string) accepts values per MATLAB's LineSpec
% type: (string) in case you want to use this for a transfer for other application

function plotConic( rp, ecc, style, type )

% set limits of conic section dependent on ecc
if ecc < 0
    disp( 'Illegal value of eccentricity! can''t be negative' );
    return;
elseif ecc < 1 % Circle or ellipse (closed orbit -- all theta possible)
    thetaMax = 180;
elseif ecc == 1
    thetaMax = 160; % Parabola: can't go quite to 180
else
    % Hyperbola
    % 1 + ecc * cosd(thetaMax) = 0
    % thetaMax = acosd( -1/ecc )
    thetaMax = 0.9 * acosd( -1/ecc ); % avoid divide by zero
end
thetaMin = -thetaMax;

% Look for a Hohmann Xfer
if strcmp(type, 'xfer')
    thetaMin = 0;
end

% independent variable theta for drawing shape
nSegments = 100;
theta = thetaMin:(thetaMax-thetaMin)/nSegments:thetaMax;

% polar to cartesian equations
r = conicRadius( rp, ecc, theta );
x = r .* cosd(theta);
y = r .* sind(theta);

% plot
plot( x, y, style)
hold on % Allows multiple plots on same figure
plot( 0, 0, '*' )
axis equal
axis square

