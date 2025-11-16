clear;
close all;
clc;

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







