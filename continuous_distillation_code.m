% DAY 2: CONTINUOUS DISTILLATION

clear;
clc;

%%
% Use polyfit to find a' b' and c' for all components
global a1 a2 a3 b1 b2 b3 c1 c2 c3 
% Use global variables so polyfit is used once

T1 = [288.15 306.30 352.35];
V1 = [107.47 109.841 116.63];
p = polyfit(T1,V1,2);
a1 = p(3); b1 = p(2); c1 = p(1);
T2 = [273.15 323.15 373.15];
V2 = [143.045 152.303 163.619];
p = polyfit(T2,V2,2);
a2 = p(3); b2 = p(2); c2 = p(1);
T3 = [303.15 353.15 400];
V3 = [107.415 113.717 120.879];
p = polyfit(T3,V3,2);
a3 = p(3); b3 = p(2); c3 = p(1);
%%

%%
% Write in the given and measured compositions, flowrates and Temperatures

Ydistillate = [0.18933462 0.53525315 0.275412231];
ErrDist = [0.137227109  0.06787682  0.205103929];
Xfeed = [.1 .45 .45];
Xbottoms = [0.027268968 0.435751717 0.536979315];
ErrBottoms = [0.005025448   0.004626404 0.000399044];
Xtray = [0.02406309 0.510471928 0.465464983];
ErrTray = [0.002805839  0.00260143  0.000204409];
Tdistillate = 372.2;
ErrTdist = 4.87903679;
Tfeed = 372.2;
ErrTfeed = 0.070710678;
Tbottoms = 372.35;
ErrTbottoms = 5.939696962;
D = 2.4; %Distillate flow rate (mL/min)
Feed = 5.4; %Feed Flow rate (mL/min)
%%

%%
% Find compositions and temperatures in the rectifying section

X = zeros(22,3); %initiate matrix for tray liquid compositions (1 is distillate)
Y = zeros(22,3); %initiate matrix for tray vapor compositions (1 is distillate)
T = zeros(22,1); %initiate vector for tray temperatures (1 is distillate)
i = 1;

Y(i,:) = Ydistillate; %Add known distillate composition
[X(i,:),T(i)] = findXT(Y(i,:)); %Find liquid composition and temp in distillate
V = 6*D/findvL(Y(1,:),T(1)); %Molar flow rate of vapor (mol/min)
L = 5/6*V; %Molar flow rate of liquid in rectifying section (mol/min)
Y(i+1,:) = (Y(i,:)*V + X(i,:)*L - Y(i,:)*L)./V; %Vapor comp leaving tray below distillate

%Loop through all trays in rectifying column stepping from stage 2 to 10
for i = 2:10 
    [X(i,:),T(i)] = findXT(Y(i,:)); %calculate temperature at tray i
    Y(i+1,:) = (Y(i,:)*V + X(i,:)*L - X(i-1,:)*L)./V; %calculate liquid comp at tray i
end
%%

%%
% Find Feed Stage Composition and Temperatures

i = 11;
[X(i,:),T(i)] = findXT(Y(i,:)); %calculate temperature
F = Feed/findvL(Xfeed,T(11));
Lbar = L + F; %Molar flow rate of liquid in stripping section (mol/min)
Vbar = V; %Vapor flow rate in stripping section (mol/min)
Y(i+1,:) = (Y(i,:)*V + X(i,:)*Lbar - X(i-1,:)*L - Xfeed*F)./Vbar; %calculate composition
%%

%%
% Find Stripping Section Compositions and Temperatures

% Loop through all stripping trays from stage 12 to 21
for i = 12:21
    [X(i,:),T(i)] = findXT(Y(i,:)); 
    Y(i+1,:) = (Y(i,:)*Vbar + X(i,:)*Lbar - X(i-1,:)*Lbar)./Vbar;
end
%%

%%
% Find Bottoms Composition and Temperature
i = 22;
[X(i,:),T(i)] = findXT(Y(i,:));
%%

%% 
% Calculate Liquid Murphree Efficiency

CompEfficiency = ((Xbottoms-X(1,:)))./((X(22,:)-X(1,:))); %Liquid Efficiencies of each component
Efficiency = sum(CompEfficiency)/3; %Average liquid efficiency
fprintf("Column Liquid Efficiency is %9.8f ",Efficiency);
%%

%% 
%Calculate Experimental and Theoretical Composition Curves

X_E = X; Y_E = Y;
for i = 2:22
    X_E(i,:) = X_E(i-1,:) + CompEfficiency.*(X(i,:)-X(i-1,:));
    Y_E(i,:) = Y_E(i-1,:) + CompEfficiency.*(Y(i,:)-Y(i-1,:));
end
Stages = 21:-1:0;
%%

%%
% Plot vapor mole fraction versus stage number

figure('DefaultAxesFontSize',12)
plot(Stages,Y(:,1),'--b','LineWidth',2)
hold on
plot(Stages,Y_E(:,1),'-b','LineWidth',2)
hold on
plot(Stages,Y(:,2),'--r','LineWidth',2)
hold on
plot(Stages,Y_E(:,2),'-r','LineWidth',2)
hold on
plot(Stages,Y(:,3),'--k','LineWidth',2)
hold on
plot(Stages,Y_E(:,3),'-k','LineWidth',2)
hold on
plot(20.9,Ydistillate(1),'b.','MarkerSize',25)
hold on
plot(21,Ydistillate(2),'r.','MarkerSize',25)
hold on
plot(21.1,Ydistillate(3),'k.','MarkerSize',25)
hold on
errorbar(20.9,Ydistillate(1),ErrDist(1),'b','LineWidth',2)
hold on
errorbar(21,Ydistillate(2),ErrDist(2),'r','LineWidth',2)
hold on
errorbar(21.1,Ydistillate(3),ErrDist(3),'k','LineWidth',2)
xlim([0 21.5])
xlabel('Stage Number')
ylabel('Vapor Mole Fraction (y_i)')
%%
 
%%
% Plot liquid mole fraction versus stage number

figure('DefaultAxesFontSize',12)
plot(Stages,X(:,1),'--b','LineWidth',2)
hold on
plot(Stages,X_E(:,1),'-b','LineWidth',2)
hold on
plot(Stages,X(:,2),'--r','LineWidth',2)
hold on
plot(Stages,X_E(:,2),'-r','LineWidth',2)
hold on
plot(Stages,X(:,3),'--k','LineWidth',2)
hold on
plot(Stages,X_E(:,3),'-k','LineWidth',2)
hold on
plot(0,Xbottoms(1),'b.','MarkerSize',25)
hold on
plot(0,Xbottoms(2),'r.','MarkerSize',25)
hold on
plot(0,Xbottoms(3),'k.','MarkerSize',25)
hold on
plot(1,Xtray(1),'b.','MarkerSize',25)
hold on
plot(1,Xtray(2),'r.','MarkerSize',25)
hold on
plot(1,Xtray(3),'k.','MarkerSize',25)
hold on
errorbar(0,Xbottoms(1),ErrBottoms(1),'b','LineWidth',2)
hold on
errorbar(0,Xbottoms(2),ErrBottoms(2),'r','LineWidth',2)
hold on
errorbar(0,Xbottoms(3),ErrBottoms(3),'k','LineWidth',2)
hold on
errorbar(1,Xtray(1),ErrTray(1),'b','LineWidth',2)
hold on
errorbar(1,Xtray(2),ErrTray(2),'r','LineWidth',2)
hold on
errorbar(1,Xtray(3),ErrTray(3),'k','LineWidth',2)
xlim([-.5 21])
xlabel('Stage Number')
ylabel('Liquid Mole Fraction (x_i)')
%%

%% 
% Plot Theoretical and Experimental Temperature Profiles
T_E = T;
for i = 1:length(T_E)
    T_E(i) = FindT(X_E(i,1),X_E(i,2),X_E(i,3));
end
figure('DefaultAxesFontSize',12)
plot(Stages,T-273.15,'--b','LineWidth',2)
hold on
plot(Stages,T_E-273.15,'-r','LineWidth',2)
hold on
plot(0,Tbottoms-273.15,'r.','MarkerSize',25)
hold on
plot(21,Tdistillate-273.15,'r.','MarkerSize',25)
hold on
plot(11,Tfeed-273.15,'r.','MarkerSize',25)
hold on
errorbar(0,Tbottoms-273.15,ErrTbottoms,'r','LineWidth',2)
hold on
errorbar(21,Tdistillate-273.15,ErrTdist,'r','LineWidth',2)
hold on
errorbar(11,Tfeed-273.15,ErrTfeed,'r','LineWidth',2)
hold on
xlim([-0.25 21.25])
xlabel('Stage Number')
ylabel('Temperature (Celsius)')
%%

%%
% The functions defined below are called on to find the values needed above

%This Calculates composition using Raoult's law and antoine's eqn
function [X,T] = findXT(Y) 
    P = 760;
    Tguess = 365;
    solver = @(X) [Y(2)*P - X(2)*findgamma(2,X(1),X(2),X(3),X(4))*antoine(2,X(4));...
        Y(3)*P - X(3)*findgamma(3,X(1),X(2),X(3),X(4))*antoine(3,X(4));...
        P - X(1)*findgamma(1,X(1),X(2),X(3),X(4))*antoine(1,X(4))...
        - X(2)*findgamma(2,X(1),X(2),X(3),X(4))*antoine(2,X(4))...
        - X(3)*findgamma(3,X(1),X(2),X(3),X(4))*antoine(3,X(4));...
        1 - X(1) - X(2) - X(3)];
    options = optimset('Diagnostics','off','Display','off');
    X1guess = Y(1)*P/antoine(1,Tguess);
    X2guess = Y(2)*P/antoine(2,Tguess);
    X3guess = Y(3)*P/antoine(3,Tguess);
    Start = [X1guess X2guess X3guess Tguess];
    Answers = fsolve(solver,Start,options);
    X(1) = Answers(1);
    X(2) = Answers(2);
    X(3) = Answers(3);
    T = Answers(4);
end
function T = FindT(x1,x2,x3) %Calculates temperature so P = 1atm
    err = 100;
    Tl = 200;
    Tu = 500;
    while err > 0.01
        Ta = (Tl+Tu)/2;
        P = x1*findgamma(1,x1,x2,x3,Ta)*antoine(1,Ta)...
        + x2*findgamma(2,x1,x2,x3,Ta)*antoine(2,Ta)...
        + x3*findgamma(3,x1,x2,x3,Ta)*antoine(3,Ta) - 760;
        err = abs(P);
        if P > 0
            Tu = Ta;
        else
            Tl = Ta;
        end
    end
    T = Ta;
end
function MolVol = findvL(X,T)
    global a1 a2 a3 b1 b2 b3 c1 c2 c3
    vL = zeros(1,3);
    vL(1) = a1 + b1*T + c1*T^2;
    vL(2) = a2 + b2*T + c2*T^2;
    vL(3) = a3 + b3*T + c3*T^2;
    MolVol = X(1)*vL(1) + X(2)*vL(2) + X(3)*vL(3);
end
function gamma = findgamma(i,x1,x2,x3,T) %Use summations to calculate gamma of each species
    X = [x1 x2 x3];
    sum1 = 0; sum3 = 0;
    for j = 1:3
        sum1 = sum1 + X(j)*findlambda(i,j,T);
    end
    for k = 1:3
        sum2 = 0;
        for  j = 1:3
            sum2 = sum2 + X(j)*findlambda(k,j,T);
        end
        sum3 = sum3 + (X(k)*findlambda(k,i,T)/sum2);
    end
    loggamma = 1 - log(sum1) - sum3;
    gamma = exp(loggamma);
end
function Psat = antoine(i,T) %Calculate Saturated pressure of each component
    T = T - 273.15; %Input is in Kelvin, but antoine uses celsius
    A1 = 6.84498; B1 = 1203.526; C1 = 222.863; 
    A2 = 6.90240; B2 = 1268.115; C2 = 216.900;
    A3 = 6.95334; B3 = 1343.943; C3 = 219.377;
    Pisat = [10^(A1-(B1/(C1+T))) 10^(A2-(B2/(C2+T))) 10^(A3-(B3/(C3+T)))];
    Psat = Pisat(i);
end
function lambda = findlambda(i,j,T) %Calculates desired lambda for 2 species at one T
    global a1 a2 a3 b1 b2 b3 c1 c2 c3
    vL = zeros(1,3);
    vL(1) = a1 + b1*T + c1*T^2;
    vL(2) = a2 + b2*T + c2*T^2;
    vL(3) = a3 + b3*T + c3*T^2;
    small = zeros(3,3);
    small(1,2) = 866.47;
    small(1,3) = -416.33;
    small(2,1) = -629.14;
    small(2,3) = -45.06;
    small(3,1) = 911.48;
    small(3,2) = 264.44;
    lambdas = zeros(3,3);
    for x = 1:3
        for y = 1:3
            lambdas(x,y) = (vL(y)/vL(x))*exp(-small(x,y)/(1.9872*T));
        end
    end
    lambda = lambdas(i,j);
end
