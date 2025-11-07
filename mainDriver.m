%------------------------------------------
% First digit is max camber as percentage of chord length
% Second digit is location of max camber in tenths of chord
% last two are max thickness as percentage of the chord
%------------------------------------------


% user input and parsing
c = input('enter chord length: ');
Naca = input('1st Digit, 2nd Digit,3rd and 4th: ','s');
panels = input('Numbers of panels: ');
C = strsplit(Naca,',');
airfoil = strcat(C{1},C{2},C{3});
% airfoil2 = str2num(airfoil);
Parts = str2double(C);

% defining variables to pass into function
m = ((Parts(1)/100)*c); 
p = ((Parts(2)/10));
t = ((Parts(3)/100)*c);

% function call to generate airfoil plot 
[XB,YB] = airfoilGen(m,p,t,c,panels,airfoil);


% % user input for alpha 
% alpha = input('Enter AoA in degrees');



