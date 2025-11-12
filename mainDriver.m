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
PartOneTask3 =1;

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
alphaVals = linspace(-20,20);

% Airfoil codes 
NACA = {'0006', '0012', '0018'};
%struct for each 
airfoils = [];



    % function call to airfoilGen with each NACA airfoil
    for i = 1: numel(NACA)
        airfoils(i).code = NACA{i};
        [airfoils(i).XB3(:,1),airfoils(i).YB3(:,1)] = airfoilGen(m2,p2,tVals(i),c3,100,NACA{i},0);
    end
    % vortex panel for cl vs. alpha -----------------
    for i =1: numel(NACA)
        
        for j = 1 : length(alphaVals)
            [airfoils(i).cl(j)] = Vortex_Panel(airfoils(i).XB3(:,1),airfoils(i).YB3(:,1),alphaVals(j));
        end
    end


    % Plotting each cl vs. alpha
    figure();
     for i = 1: numel(NACA)
       plot(alphaVals,airfoils(i).cl,'DisplayName',strcat('NACA',NACA{i}));
        hold on; 
    end
    xlabel('Sectional Angle of Attack [Degrees]');
    ylabel('Sectional Lift Coefficient')
    yline(0, '-', 'Zero-Lift Line', 'LineWidth', 1.5, ...
    'LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');
    title('Coefficient of Lift vs. Angle of Attack');
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









