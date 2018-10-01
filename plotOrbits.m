function plotOrbits

hold on % Now we can do multiple plot commands in the same window
axis equal % The figure will be square
axis( [-10 10 -10 10] ) % Go from -10 to 10 Earth radii (R0 = Earth radius)

% Draw the Earth as a black circle with cross
angles = 0:10:360; % whole circle, every 10 degrees
xvals = cosd(angles);
yvals = sind(angles);
plot( xvals, yvals, 'k-' ) % Draw circle; 'k-' does black solid lines -- do 'help plot'
plot( [0 0], [-1 1], 'k-' ) % Draw vertical stroke of cross
plot( [-1 1], [0 0], 'k-' ) % Draw horizontal stroke of cross

% Draw a family of Earth orbits, all with the same periapsis radius

rp = 1.5; % Periapsis radius, in units of R0

for ecc=0:0.25:2 % Do a circle, a few ellipses, a parabola, and some hyperbolas
 
    % Plot the orbit, in magenta for circle, blue for ellipse, red for
    % parabola, green for hyperbola
    if ecc==0
        style = 'm-';
    elseif ecc<1
        style = 'b-';
    elseif ecc==1
        style = 'r-';
    else
        style = 'g-';
    end
    plotConic( rp, ecc, style, 'xfer' );
        
end
