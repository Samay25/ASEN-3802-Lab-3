function [e,c_L,c_Di] = PLLT2(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
%{ 
PLLT takes in where e is the span efficiency factor (to be computed and returned), C_L is the coefficient of lift (to be
computed and returned), C_Di is the induced coefficient of drag (to be computed and reported), 

b is the span (in feet), a0_t is the cross-sectional lift slope at the tips (per radian), a0_r is the cross-sectional lift
slope at the root (per radian), c_t is the chord at the tips (in feet), c_r is the chord at the root (in feet),
aero_t is the zero-lift angle of attack at the tips (in radians), aero_r is the zero-lift angle of attack at the
root (in radians), geo_t is the geometric angle of attack (geometric twist + α) at the tips (in radians),
geo_r is the geometric angle of attack (geometric twist + α) at the root (in radians), and N is the number
7 
%}

% theta linspace
theta_vals = zeros(N,1);

for i = 1: (N)
    theta_vals(i) = ((i)*pi)/(2*N);
end

% y map to theta
y = -b/2.*cos(theta_vals);

% linear variation of lift slopes 
ma0 = (a0_t - a0_r)/-(b/2);
a0_y = ma0*y + a0_r;

% linear variation of chord length
mc = (c_t - c_r)/-(b/2);
c_y = mc*y + c_r;

% linear variation of zero lift AoA
maero = (aero_t - aero_r)/-(b/2);
aero_y = maero*y + aero_r;

% linear variation of geometric twist
mgeo = (geo_t - geo_r)/-(b/2);
geo_y = mgeo*y + geo_r;

% create matricies 

% preallocation 
M = zeros(N,N);
d = zeros(N,1);

for i = 1: N 
    for j = 1: N
      M(i,j) = ( (4*b) / (a0_y(i).*c_y(i)) ) .* (sin((2*j - 1).*theta_vals(i))) +...
          ((2*j - 1) .* ((sin((2*j - 1).*theta_vals(i)))/(sin(theta_vals(i))))); 
    end
    d(i) = geo_y(i) - aero_y(i);
end


% solve system for coefficients
x = M\d;

% compute delta

NVals = 2*(1:N) - 1; % only odd NVals needed
del = zeros(N,1);
for i = 2: N
     del(i) = NVals(i) * ((x(i)./x(1))).^2;
end

delta = sum(del);

% Computing Aspect Ratio

    % Compute Taper Ratio
    taperRatio = c_t/c_r;
    
    % Compute Planform Aera
    S = (b/2)*c_r*(1+taperRatio);
    
    AR = b^2 / S; 


% Compute cl 
c_L = x(1)* pi .* AR;

% Compute cd

c_Di = (c_L.^2 ./(pi.*AR)) .* (1+delta);

% span efficiency 
e = 1 / (1 + delta);

 

end


