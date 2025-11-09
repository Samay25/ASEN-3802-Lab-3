function [xcoords,ycoords] = airfoilGen(m,p,t,c,panels,airfoil,plots)
%AIRFOILGEN creates a figure plotting x,y coordinates of the airfoil by taking in
%naca four digit info. chord length and number of panels
%Put paramter plots equal to one if you want plot to show, if not 0 

if m == 0 
    y_t = zeros(panels,1);
    y_c = zeros(panels,1);
    dy_c = zeros(panels,1);
    xi = zeros(panels,1);
    x_U = zeros(panels,1);
    x_L = zeros(panels,1);
    y_U = zeros(panels,1);
    y_L = zeros(panels,1);
    x_Vals = linspace(0,c,panels);
    
    
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
    
    flip_x_L = flip(x_L);
    flip_y_L = flip(y_L);
    
    flip_x_L(panels) = [];
    flip_y_L(panels) = [];
    
    airfoil_x_coords = [flip_x_L;x_U];
    airfoil_y_coords = [flip_y_L;y_U];
    
    xcoords = airfoil_x_coords;
    ycoords = airfoil_y_coords;
    
    figtitle = strcat('NACA ' , airfoil);
    if plots == 1
        figure();
        plot(airfoil_x_coords,airfoil_y_coords,LineWidth=.5,DisplayName='Airfoil');
        xlabel('Ditance along chord [m]');
        ylabel('Height from chord line [m]');
        title(figtitle);
        legend('Location','best');
        ylim([-8,8]);
    end
else



        y_t = zeros(panels,1);
        y_c = zeros(panels,1);
        dy_c = zeros(panels,1);
        xi = zeros(panels,1);
        x_U = zeros(panels,1);
        x_L = zeros(panels,1);
        y_U = zeros(panels,1);
        y_L = zeros(panels,1);
        x_Vals = linspace(0,c,panels);
        
        
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
        
        flip_x_L = flip(x_L);
        flip_y_L = flip(y_L);
        
        flip_x_L(panels) = [];
        flip_y_L(panels) = [];
        
        airfoil_x_coords = [flip_x_L;x_U];
        airfoil_y_coords = [flip_y_L;y_U];
        
        xcoords = airfoil_x_coords;
        ycoords = airfoil_y_coords;

        camberline = y_c;
        
        figtitle = strcat('NACA ' , airfoil);
        if plots == 1
            figure();
            plot(airfoil_x_coords,airfoil_y_coords,LineWidth=.5,DisplayName='Airfoil');
            hold on; 
            plot(x_Vals,camberline,LineWidth=.5,DisplayName='Camberline');
            xlabel('Ditance along chord [m]');
            ylabel('Height from chord line [m]');
            title(figtitle);
            legend('Location','best');
            ylim([-8,8]);
        end
    end

end

