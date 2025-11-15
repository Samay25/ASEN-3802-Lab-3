clear; close all; clc; 
%---------------------
% PART 2 TASK 2
%---------------------

% % units of meters and radians
% 
% 
% %% First Debugging Scenario
% b = 100; a0_t = 2*pi; a0_r = 2*pi; c_t = 10; c_r = 10; aero_t = 0; aero_r = 0; geo_t = 5*pi/180; geo_r = 5*pi/180; N = 5;
% [e, c_L, c_Di] = PLLT2(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N);

% % should output the following ----------
% 
% %{
% e = 0.9227
% 
% 
% c_L = 0.4402
% 
% 
% c_Di = 0.0067
% 
% %}
% 
% % ----------------------------------------
% 
% 
% 
% %% Second Debugging Scenario
 % b = 100; a0_t = 2*pi; a0_r = 2*pi; c_t = 8; c_r = 10; aero_t = 0; aero_r = 0; geo_t = 5*pi/180; geo_r = 5*pi/180; N = 5;
% [e2, c_L2, c_Di2] = PLLT2(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,50);
% 
% % should output the following ----------
% 
% %{
% e = 0.9430
% 
% 
% c_L = 0.4534
% 
% 
% c_Di = 0.0062
% %}
% 
% % ----------------------------------------
% %% Third Debugging Scenario
% b = 100; a0_t = 6.3; a0_r = 6.5; c_t = 8; c_r = 10; aero_t = 0; aero_r = -2*pi/180; geo_t = 5*pi/180; geo_r = 7*pi/180; N = 5;
% [e3, c_L3, c_Di3] = PLLT2(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,50);
% 
% % should output the following ----------
% 
% %{
% e = 0.9795
% 
% 
% c_L = 0.6674
% 
% 
% c_Di = 0.0130
% 
% %}
% 
% % ----------------------------------------

%% Plotting Induced Drag vs. Taper Ratio

% Aspect ratio space
AR = 4:2:10;
% Taper ratio space
TaperRatio = linspace(0,1);
% Induced drag and induced drag factor space
c_Di2 = zeros(length(AR),length(TaperRatio));
InducedDragFactor =zeros(length(AR),length(TaperRatio));

% wing characteristics besides chord
b = 100; a0_t = 2*pi; a0_r = 2*pi; aero_t = 0; aero_r = 0; geo_t = 5*pi/180; geo_r = 5*pi/180; N = 5;


% loop through to get induced drag values with constant aspect ratios 
for i = 1: length(AR)

    tempAR = AR(i);

    for j =1 : length(TaperRatio) 

        tempTapeRat = TaperRatio(j);

        % compute tip and root chord lengths 
         c_r2 = 2*b / (tempAR*(1+tempTapeRat));
         c_t2 = tempTapeRat * c_r2;

        % compute induced drag factor for fixed AR wing for 50 terms
        [e,~,~] = PLLT2(b,a0_t,a0_r,c_t2,c_r2,aero_t,aero_r,geo_t,geo_r,50);
        InducedDragFactor(i,j) = (1/e)-1;
    end
end



figure();
for i = 1 : length(AR)
    plot(TaperRatio,InducedDragFactor(i,:),'DisplayName',sprintf('AR = %d',AR(i)));
    hold on;
end 
xlabel('$\lambda = \frac{c_t}{c_r}$','Interpreter','latex');
ylabel('$\delta = \frac{1}{e} - 1 $','Interpreter','latex');
ylim([0,.16]);
title('Induced Drag Factor vs. Taper Ratio');
legend('Location','best');







