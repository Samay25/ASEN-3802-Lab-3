clear;
close all;
clc;

%% Pretty plots set up
% set(0, 'defaultFigureUnits', 'inches', 'defaultFigurePosition', [1 1 8 5]);
% figures are 8" wide and 5" tall, with the bottom left corner of the
%figure beginning 1" up, and 1" to the right from the bottom left corner
%of your screen -- adjust the size of the figure to your liking
set(0,'defaultLineLineWidth',2.5) % sets all line widths of plotted lines
set(0,'DefaultaxesLineWidth', 1.5) % sets line widths of axes
set(0,'DefaultaxesFontSize', 14)
set(0,'DefaultTextFontSize', 14)
set(0,'DefaultaxesFontName', 'Times new Roman')
set(0,'DefaultlegendFontName', 'Times new Roman')
set(0,'defaultAxesXGrid','on')
set(0,'defaultAxesYGrid','on')
%... there are semi-infinite options here!

%------------------------------------------
% First digit is max camber as percentage of chord length
% Second digit is location of max camber in tenths of chord
% last two are max thickness as percentage of the chord
%------------------------------------------

%% Plot suppressions 
PartOneTask1 = 0;
PartOneDel1 = 0;
PartOneDel2 = 0;
PartOneTask3 =0;

PartTwoTask2 = 0; 

Part3_all = 1; 



% -----------------------------------------
%% PART 1 TASK 1 
% -----------------------------------------
if PartOneTask1 == 1
    % user input and parsing
c = input('enter chord length: ');
Naca = input('1st Digit, 2nd Digit,3rd and 4th: ','s');
panels = input('Numbers of panels: ');
C = strsplit(Naca,',');
airfoil = strcat(C{1},C{2},C{3});
Parts = str2double(C);

% defining variables to pass into function
m = ((Parts(1)/100)); 
p = ((Parts(2)/10));
t = ((Parts(3)/100));

% function call to generate airfoil plot 
[XB,YB] = airfoilGen(m,p,t,c,panels,airfoil,1);
end

% -----------------------------------------
%% PART 1 TASK 2
% -----------------------------------------
if PartOneDel1 == 1
    % -------------------Deliverable 1----------------------
% Params for NACA0012 at 5 degrees angle of attack
c2 = 5; 
m2 = 0;
p2 = 0; 
t2 =.12; 
alpha = 5; 

% length of numPanels to test
numPanels = 10:1:200;
cl = zeros(length(numPanels),1);

% looping through panels and getting corresponding cl 
for i = 1: length(numPanels)

    % function call to airfoilGen with NACA0012 info
[XB2,YB2] = airfoilGen(m2,p2,t2,c2,numPanels(i),'0012',0);

% function call to vortex panel to 
[cl(i)] = Vortex_Panel(XB2,YB2,alpha);
end

% max CL from cl vs. num panels
CL_Max = max(cl);
Clstr = strcat('Cl = ',num2str(CL_Max));

% Loop to get error and determin number of panels to achieve cl within 1
% percent
for i = 1: length(cl)
err(i) = (abs(CL_Max - cl(i)) / CL_Max)*100;
end

[~,numPanelsForExact] = find(err<1,1);


% plotting cl vs. number of panels 
figure();
plot(numPanels,cl);
xlabel('Number of Panels');
ylabel('Predicted Coeff. Lift');
yline(CL_Max,'--',{'Cl Conv.',Clstr},LineWidth=1.5,LabelHorizontalAlignment='left');
ylim([.56,.63]);
title('Predicted Cl vs. Number of Panels');
% -----------------------------------------
end

if PartOneDel2 == 1
      % ----------------Deliverable 2-------------------------
%Values for NACA 0006, 0012, 0018
c3 = 10; 
m2 = 0;
p2 = 0; 
tVals =[.06,.12,.18]; 
alpha = 5; 
alphaVals = -20:5:20;

% Airfoil codes 
NACA = {'0006', '0012', '0018'};
%struct for each 
airfoils = [];



    % function call to airfoilGen with each NACA airfoil
    for i = 1: numel(NACA)
        airfoils(i).code = NACA{i};
        [airfoils(i).XB3(:,1),airfoils(i).YB3(:,1)] = airfoilGen(m2,p2,tVals(i),c3,100,NACA{i},0);
    end
    % ----------testing vortex panel for cl vs. alpha -----------------
    for i =1: numel(NACA)
        for j = 1 : length(alphaVals)
            [airfoils(i).cl(j)] = Vortex_Panel(airfoils(i).XB3(:,1),airfoils(i).YB3(:,1),alphaVals(j));
        end
    end

    %--------------------------------------------------------------------

    % Plotting each cl vs. alpha
    figure();
    plot(alphaVals,airfoils(1).cl,'DisplayName',strcat('NACA',NACA{1}));
    hold on; 
    plot(alphaVals,airfoils(2).cl,'DisplayName',strcat('NACA',NACA{2}));
    hold on; 
    plot(alphaVals,airfoils(3).cl,'DisplayName',strcat('NACA',NACA{3}));
    yline(0,'-.',{'Zero Lift Line'},'HandleVisibility', 'off');
    xline(0,'-.',{'alpha(L=0) = 0'},'LabelHorizontalAlignment','Left','HandleVisibility', 'off');
    xlabel('AoA');
    ylabel('Cl')
    title('Cl vs. AoA');
    legend('Location','best');



%-----------------------------------------

end


% -----------------------------------------
%% PART 1 TASK 3
% -----------------------------------------

if PartOneTask3 ==1 
    % ------------------------------------------------------------------------

%Values for NACA 0012, 2412, 4412
c3 = 10; 
mVals = [0,.02,.04];
pVals = [0,.4,.4]; 
tVals =[.12,.12,.12]; 
alpha = 5; 
alphaVals = -30:1:30;

% Airfoil codes 
NACA = {'0012','2412','4412'}; % need to be same order as t vals 

%struct for each 
airfoils = [];



    % function call to airfoilGen with each NACA airfoil
    for i = 1: numel(NACA)
        airfoils(i).code = NACA{i};
        [airfoils(i).XB3(:,1),airfoils(i).YB3(:,1)] = airfoilGen(mVals(i),pVals(i),tVals(i),c3,100,NACA{i},0);
    end
    % ----------testing vortex panel for cl vs. alpha -----------------
    for i =1: numel(NACA)
        for j = 1 : length(alphaVals)
            [airfoils(i).cl(j)] = Vortex_Panel(airfoils(i).XB3(:,1),airfoils(i).YB3(:,1),alphaVals(j));
            airfoils(i).a0 = mean(diff(airfoils(i).cl));
        end
    end

    %--------------------------------------------------------------------

    % Plotting each cl vs. alpha
    figure();
    for i = 1: numel(NACA)
       plot(alphaVals,airfoils(i).cl,'DisplayName',strcat('NACA',NACA{i}));
        hold on; 
    end
    xlabel('Sectional Angle of Attack [Degrees]');
    ylabel('Sectional Lift Coefficient')
    yline(0,'-.',{'Zero Lift Line'},'HandleVisibility', 'off');
    xline(-2,'-.',{'alpha(L=0) is -2','for NACA 2412'},'LabelHorizontalAlignment','Right','HandleVisibility', 'off');
    xline(-4,'-.',{'alpha(L=0) is -4','for NACA 4412'},'LabelHorizontalAlignment','Left','HandleVisibility', 'off');
    title('Cl vs. AoA');
    legend('Location','best');



%-----------------------------------------
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if PartTwoTask2 ==1
% -----------------------------------------
% PART 2 TASK 2
% -----------------------------------------

% % units of meters and radians

%% First Debugging Scenario
b1 = 100; a0_t1 = 2*pi; a0_r1 = 2*pi; c_t1 = 10; c_r1 = 10; aero_t1 = 0; aero_r1 = 0; geo_t1 = 5*pi/180; geo_r1 = 5*pi/180; N1 = 5;
[e1, c_L1, c_Di1] = PLLT2(b1,a0_t1,a0_r1,c_t1,c_r1,aero_t1,aero_r1,geo_t1,geo_r1,N1);

% should output the following below ----------

%{
    e = 0.9227
    
    
    c_L = 0.4402
    
    
    c_Di = 0.0067

%}

% ----------------------------------------

%% Second Debugging Scenario
 b2 = 100; a0_t2 = 2*pi; a0_r2 = 2*pi; c_t2 = 8; c_r2 = 10; aero_t2 = 0; aero_r2 = 0; geo_t2 = 5*pi/180; geo_r2 = 5*pi/180; N2 = 5;
[e2, c_L2, c_Di2] = PLLT2(b2,a0_t2,a0_r2,c_t2,c_r2,aero_t2,aero_r2,geo_t2,geo_r2,N2);

% should output the following ----------

%{
e = 0.9430


c_L = 0.4534


c_Di = 0.0062
%}

% ----------------------------------------

%% Third Debugging Scenario
b3 = 100; a0_t3 = 6.3; a0_r3 = 6.5; c_t3 = 8; c_r3 = 10; aero_t3 = 0; aero_r3 = -2*pi/180; geo_t3 = 5*pi/180; geo_r3 = 7*pi/180; N3 = 5;
[e3, c_L3, c_Di3] = PLLT2(b3,a0_t3,a0_r3,c_t3,c_r3,aero_t3,aero_r3,geo_t3,geo_r3,N3);

% should output the following ----------

%{
e = 0.9795


c_L = 0.6674


c_Di = 0.0130

%}

% ----------------------------------------

%% Plotting Induced Drag vs. Taper Ratio

% Aspect ratio space
AR = 4:2:10;
% Taper ratio space
TaperRatio = linspace(0,1);
% Induced drag and induced drag factor space
InducedDragFactor =zeros(length(AR),length(TaperRatio));

% wing characteristics besides chord
bP = 100; a0_tP = 2*pi; a0_rP = 2*pi; aero_tP = 0; aero_rP = 0; geo_tP = 5*pi/180; geo_rP = 5*pi/180; NP = 50;


% loop through to get induced drag values with constant aspect ratios 
for i = 1: length(AR)

    tempAR = AR(i);

    for j =1 : length(TaperRatio) 

        tempTapeRat = TaperRatio(j);

        % compute tip and root chord lengths 
         c_r2P = 2*bP / (tempAR*(1+tempTapeRat));
         c_t2P = tempTapeRat * c_r2P;

        % compute induced drag factor for fixed AR wing with N = 50 terms
        [eP,~,~] = PLLT2(bP,a0_tP,a0_rP,c_t2P,c_r2P,aero_tP,aero_rP,geo_tP,geo_rP,NP);
        InducedDragFactor(i,j) = (1/eP)-1;
    end
end



figure();
for i = 1 : length(AR)
    plot(TaperRatio,InducedDragFactor(i,:),'DisplayName',sprintf('AR = %d',AR(i)),LineWidth=2);
    hold on;
end 
xlabel('$\frac{c_t}{c_r}$','Interpreter','latex');
ylabel('$\delta $','Interpreter','latex');
ylim([0,.16]);
title('Induced Drag Factor vs. Taper Ratio');
legend('Position',[0.374791666666667 0.391435880447585 0.0590277777777777 0.0747938751472319]);

end

if Part3_all == 1 
    % -------------------------------------
    % Part 3 Task 1 
    % -------------------------------------
    
    % Using your NACA airfoil generator, the provided vortex panel code, and your PLLT code, generate a plot
    % of coefficient of lift versus angle of attack for the Cessna 180 aircraft
    
    % The Cessna uses a C180 (Can approx with trapezoidal)
    % Span of 36ft 
    % Straight taper from 5ft 4in root chord to 3ft 7in Tip chord
    % ROOT airfoil NACA 2412 
    % TIP airfoil NACA 0012
    % Geometric AoA varies linearly and is 2 degrees larger at the root vs. tip
    
    % define/find params of wing design to pass into functions 
    
    b = 36; % span
    c_r = 5.3333; % root chord
    c_t = 3.5833; % tip chord 
    aero_t = 0; % zero lift AoA at tip
    aero_r = -.1982; % zero lift AoA at root
    a0_t =  .1051; % zero lift slope at tip 
    a0_r = .1119; % zero lift slope at root
    alphaVals = linspace(-30,30,120); % AoA values
    N = 50; % number of terms 
    
    taperRatio = c_t/c_r;
    planformArea = (b/2)*c_r*(1+taperRatio);
    AR = b^2 / planformArea; 
    
    % pre allocate
    geo_t = zeros(length(alphaVals),1);
    geo_r = zeros(length(alphaVals),1);
    
    
    for i = 1: length(alphaVals)
        geo_t(i) = alphaVals(i);
        geo_r(i) = alphaVals(i) + 2; 
    end
    
    % preallocation
    c_L = zeros(length(alphaVals),1);
    c_D = zeros(length(alphaVals),1);
    e = zeros(length(alphaVals),1);
    
    for i = 1: length(alphaVals)
    
        [e(i),c_L(i),c_D(i)] = PLLT2(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t(i),geo_r(i),N);
    
    end
    
    % figure();
    plot(alphaVals,c_L,'DisplayName','Cl');
    xlabel('$\alpha [^\circ]$','Interpreter','latex');
    ylabel('$ C_{l} $','Interpreter','latex');
    yline(0,'--',{'Zero Lift Line'},'HandleVisibility', 'off')
    legend('Position',[0.334166666666667 0.67978354978355 0.0946428571428571 0.0440476190476189]);
    title(' Coefficient of Lift Vs. Angle of Attack ');
    
    % -------------------------------------
    % Part 3 Task 2
    % -------------------------------------
    
    % using exp data ----------------------
    
    % digitize plots from exp data
    data = readmatrix("Re_6x106.csv");
    
    % sectional values from exp data 
    Cl = data(:,1);
    Cd = data(:,2);
    
    % map Cl to AoA (cl = 2piAoA) and convert to degrees
    alpha = Cl / (2*pi) * (180/pi);
    
    
    % polyfit to get 'function' of cd vs. alpha 
    [p,S] = polyfit(alpha,Cd,2);
    
    % space for new AoA values 
    alphaVals = linspace(-30,30);
    
    % creating function for Cd vs. Alpha 
    y1 = polyval(p,alphaVals); 
    
    % here I was trying to plot with fitted data only (its basically the same as the one below)
    % plotting model alongside experimental data 
    % figure();
    % plot(alpha,Cd,'DisplayName','Exp. Data');
    % hold on; 
    % plot(alphaVals,y1,'--','DisplayName','Model');
    % xlabel('Angle of attack $ \alpha [^\circ] $',Interpreter='latex');
    % ylabel('$ C_{d}$',Interpreter='latex');
    % legend('Position',[0.428214285714286 0.606904761904762 0.163392857142857 0.0797619047619047]);
    % title(' $ C_{d}$ vs. $\alpha $',Interpreter='latex');
    
    % -----------------------------------
    
    % using sectional lift coeff --------
    
    % polyfit to get Cd as function of Cl
    [p2,S2] = polyfit(Cl,Cd,2);
    
    % eval function with sectional lift coeff from PLLT to get sectional drag
    % ('profile' drag)
    y2 = polyval(p2,c_L);
    
    
    % get alpha values from c_L
    alphaVals2 = c_L / (2*pi) * (180/pi);
    
    %plotting sectional drag coeff. vs AoA
    figure();
    plot(alphaVals2,y2,'DisplayName','Model');
    hold on; 
    plot(alpha,Cd,'DisplayName','Exp Data');
    xlabel('Angle of attack $ \alpha [^\circ] $',Interpreter='latex');
    ylabel('$ C_{d}$',Interpreter='latex');
    legend('Position',[0.428214285714286 0.606904761904762 0.163392857142857 0.0797619047619047]);
    title('Coefficient of Drag Vs. Angle of Attack ');
    % -----------------------------------
    
    % -------------------------------------
    % Part 3 Task 3
    % -------------------------------------
    
    % total drag coeff (profile [from fitted data] and induced [from PLLT])
    C_Dtot = y2 + c_D;
    
    
    % plotting total drag, induced and profile vs AoA
    figure();
    plot(alphaVals2,C_Dtot,'DisplayName','Total Drag');
    hold on; 
    plot(alphaVals2,c_D,'DisplayName','Induced Drag');
    plot(alphaVals2,y2,'DisplayName','Profile Drag');
    xlabel('Angle of attack $ \alpha [^\circ] $',Interpreter='latex');
    ylabel('$ C_{d}$',Interpreter='latex');
    legend('Position',[0.408472222222222 0.442443070278759 0.0763888888888887 0.0571260306242636]);
    title('Effect of Angle of Attack on Drag Coefficient');
    
    % -------------------------------------
    % Part 3 Task 4
    % -------------------------------------
    
    % given info 
    weight = 2500; % [lbs]
    altitude = 10000; % [ft]
    
    % getting profile drag coeff at zero AoA
    x = min(abs(y2));
    CD0 = x; 
    
    % linspace for airspeed
    v = linspace(10,500); % [ft/s]
    v2 = v/1.688; % [knots]
    
    % standard day (from anderson)
    rho_inf = 1.7556 * 10^(-3); % [slugs/ft^3]
    p_inf = 1.4556 * 10^3; % [lb/ft^2]
    
    % steady level flight so, lift = Weight and ThrustReq = Drag
    q_inf = (.5 * rho_inf) .* v.^2;
    k = 1 / (pi * e(49) * AR);
    
    thrustReq = (q_inf .* CD0 .* planformArea) + ((k * weight^2) ./ (q_inf .* planformArea));
    thrustReqParasite = (q_inf .* CD0 .* planformArea);
    thrustReqInduced = ((k * weight^2) ./ (q_inf .* planformArea));
    
    [minThrustReq,I] = min(thrustReq); % speed for min thrust required
    Vel_LDmax = v2(I); % speed for L/D max 
    
    fprintf('Minimum Thrust Required is %4.2f [lbf] \n',minThrustReq);
    fprintf('Velocity for Minimum Thrust Required is  %4.2f [knots] \n',Vel_LDmax);
    
    
    
    figure();
    plot(v2,thrustReq,'DisplayName','Thrust Required');
    hold on; 
    plot(v2,thrustReqInduced,'DisplayName','Thrust Required for Induced Drag');
    plot(v2,thrustReqParasite,'DisplayName','Thrust Required for Parasite');
    xlabel('Airspeed [knots]');
    ylabel('Thrust Required [lbf]');
    ylim([0,2.5*10^3]);
    xlim([min(v2),max(v2)]);
    title('Thrust required vs. Airspeed');
    legend('Position',[0.382876984126985 0.449781535700261 0.38125 0.11547619047619]);







end 







