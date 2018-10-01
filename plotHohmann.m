function plotHohmann(hi,hf,ip,ia)
%PLOTHOHMANN Plots a Hohmann transfer
%   Plots Hohmann transfer and plane changes. Outputs figure as well as
%   perigee speed, apogee speed, delta v1, and delta v2.

if nargin == 2
    ip = 0;
    ia = 0;
end;

ri = (hi + 6378)/6378;
rf = (hf + 6378)/6378;

hold on % Now we can do multiple plot commands in the same window
axis equal % The figure will be square
axis( [-10 10 -10 10] ) % Go from -10 to 10 Earth radii (R0 = Earth radius)

% Draw the Earth as a black circle with cross
angles = 0:10:360; % whole circle, every 10 degrees
xvals = cosd(angles);
yvals = sind(angles);
plot( xvals, yvals, 'k-' ) % Draw circle; 'k-' does black solid lines -- do 'help plot'
plot( ri*xvals, ri*yvals*cosd(ip), 'b')%Draw first orbit with plane change 
plot( rf*xvals, rf*yvals*cosd(ip), 'r') %Draw second orbit with plane change
plot( [0 0], [-1 1], 'k-' ) % Draw vertical stroke of cross
plot( [-1 1], [0 0], 'k-' ) % Draw horizontal stroke of cross

ecc = (rf - ri)/(rf + ri);
plotConic( ri, ecc, 'g-', 'xfer');

axis( [-10 10 -10 10] )

a = 6378*(rf + ri)/2;
pSpeed = sqrt(398600*((2/(6378*ri))-(1/a)));%Perigee Speed
aSpeed = sqrt(398600*((2/(6378*rf))-(1/a)));%Apogee Speed
if nargin == 4
    oiSpeed = sqrt(398600/(6378*ri));%Initial Orbit speed
    ofSpeed = sqrt(398600/(6378*rf));%Final Orbit speed
    dv1 = sqrt((pSpeed^2)+(oiSpeed^2)-(2*pSpeed*oiSpeed*cosd(ip)));
    dv2 = sqrt((aSpeed^2)+(ofSpeed^2)-(2*aSpeed*ofSpeed*cosd(ia)));
else
    dv1 = pSpeed - sqrt(398600/(6378*ri));
    dv2 = sqrt(398600/(6378*rf)) - aSpeed;
end;

fprintf('Perigee Speed: %f\nApogee Speed: %f\nDelta V1: %f\nDelta V2: %f\n',pSpeed,aSpeed,dv1,dv2);




end


