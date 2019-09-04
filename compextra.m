%We use Gm instead of m to convenience in km^3/s^2
Gmsun = 132712e6;
Gme = 398600;
Gmsat = 37.931e6;
%Earth to Sun sphere of influence in km
re2sun = 1496e5; 
soie2sun = re2sun/(((Gme/Gmsun)^(-2/5))+1);
%Saturn to Sun sphere of influence in km
rsat2sun = 14294e5;
soisat2sun = rsat2sun/(((Gmsat/Gmsun)^(-2/5))+1);
%Case I
disp('Case I: Elliptical earth orbit')
ra = 46045.953; %From printout 1 in  km
d1 = soie2sun-ra; %km
if d1 <= 0
    disp('The spacecraft will enter a heliocentric orbit');
else
    disp('The spacecraft will stay inside the sphere of influence')
    fprintf('Range = %6.3f km\n',d1);
end
%Case III
disp('Case III: Elliptical earth orbit')
a3 = -1400856.745; %From printout 3 in  km
Energy3 = -Gmsat/(2*a3); %km
v3 = sqrt(Gmsat*((2/soisat2sun)-(1/a3))); %km/s
fprintf('Cassini would travel at %6.3f km/s at the edge of Saturn SOI in Saturn frame\n', v3)
%Case V
disp('Case V: Elliptical earth orbit')
ra = 6890.566; %From printout 5 in km
d3 = soie2sun-ra; %km
if d3 <= 0
    disp('The spacecraft will enter a heliocentric orbit');
else
    disp('The spacecraft will stay inside the sphere of influence')
    fprintf('Range = %6.3f km\n',d3);
end