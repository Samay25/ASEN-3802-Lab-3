clc
close all
clear all


m = input('Enter max camber value: ');
p = input('Enter location of max camber: ');
t = input('Enter thickness ratio value: ');
c = input('Enter chord length: ');
panels = input('Enter number of panels: ');

y_t = zeros(panels,1);
y_c = zeros(panels,1);
dy_c = zeros(panels,1);
xi = zeros(panels,1);
x_U = zeros(panels,1);
x_L = zeros(panels,1);
y_U = zeros(panels,1);
y_L = zeros(panels,1);
x_Vals = linspace(0,c,panels);

y_t_clean = rmmissing(y_t);
y_c_clean = rmmissing(y_c);
dy_c_clean = rmmissing(dy_c);
xi_clean = rmmissing(xi);
x_U_clean = rmmissing(x_U);
x_L_clean = rmmissing(x_L);
y_U_clean = rmmissing(y_U);
y_L_clean = rmmissing(y_L);


for i = 1:panels

    x = x_Vals(i);
    
    if x < (p*c)
        y_c(i) = m * (x/(p)^2) * ((2*p) - (x/c));
        dy_c(i) = (m/(c*p)^3) * (-(2 * x * p) + (2 * c * p^2));

    else
        y_c(i) = m * (((c-x)/(1-p)^2) * (1 + (x/c) - (2*p)));
        dy_c(i) = (m/(c*(1-p)^3)) * ((-2*x) + (2*x*p) + (2*c*p) - (2*c*p^2));

    end
    sum = ((0.2969 * sqrt(x/c)) - (0.1260 * (x/c)) - (0.3516 * (x/c)^2) + (0.2843 * (x/c)^3) - (0.1036 * (x/c)^4));
    y_t(i) = c * (t/0.2) * sum;
  
    xi(i) = atan(dy_c(i));
    x_U(i) = x - (y_t(i) * sin(xi(i)));
    x_L(i) = x + (y_t(i) * sin(xi(i)));
    y_U(i) = y_c(i) + (y_t(i) * cos(xi(i)));
    y_L(i) = y_c(i) - (y_t(i) * cos(xi(i)));

end

x_first_final = c;
y_first_final = 0;

airfoil_x_coords = [x_first_final; x_L; x_U; x_first_final];
airfoil_y_coords = [y_first_final; y_L; y_U; y_first_final];

figure;
plot(airfoil_x_coords,airfoil_y_coords);
ylim([-8,8]);

