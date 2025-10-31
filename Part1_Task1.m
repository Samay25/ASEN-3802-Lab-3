clc
close all
clear all


m = input('Enter max camber value: ');
p = input('Enter location of max camber: ');
t = input('Enter thickness ratio value: ');
c = input('Enter chord length: ');
panels = input('Enter number of panels: ');

y_t = zeros;
y_c = zeros;
dy_c = zeros;

for x = 0:(p*c)
    
    y_c = m * (x/(p)^2) * ((2*p) - (x/c));
    y_t = c * (t/0.2) * ((0.2969 * sqrt(x/c)) - (0.1260 * (x/c)) - (0.3516 * (x/c)^2) + (0.2843 * (x/c)^3) - (0.1036 * (x/c)^4));
    dy_c = (m/(c*p)^3) * (-(2 * x * p) + (2 * c * p^2));
    y_U = y_c + (y_t * cos(xi));
    y_L = y_c - (y_t * cos(xi));

end

for x = (p*c):c
    
    y_c = m * (((c-x)/(1-p)^2) * (1 + (x/c) - (2*p)));
    y_t = c * (t/0.2) * ((0.2969 * sqrt(x/c)) - (0.1260 * (x/c)) - (0.3516 * (x/c)^2) + (0.2843 * (x/c)^3) - (0.1036 * (x/c)^4));
    dy_c = (m/(c*(1-p)^3)) * ((-2*x) + (2*x*p) + (2*c*p) - (2*c*p^2));
    y_U = y_c + (y_t * cos(xi));
    y_L = y_c - (y_t * cos(xi));

end

xi = zeros(n);
xi = atan(dy_c);

x_first_final = c;
y_first_final = 0;

x_U = x - (y_t * sin(xi));
x_L = x + (y_t * sin(xi));
y_U = y_c + (y_t * cos(xi));
y_L = y_c - (y_t * cos(xi));