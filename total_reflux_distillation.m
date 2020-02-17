%%DAY 1: TOTAL REFLUX DISTILLATION

clear;
clc;
%% Use polyfit to find a' b' and c' for all components
global a1 a2 a3 b1 b2 b3 c1 c2 c3 %Use global variables so polyfit is used once
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
%% Given and measured compositions
Xbottoms = [.05 0.46 0.49];
ErrBottoms = [0.172719196 0.07875159 0.093967603];
Ydistillate = [0.948263174  0.030486186 0.02125064];
ErrDist = [0.003723788  0.002591978 0.00113181];
Xtray = [0.128760483    0.500736977 0.370502539];
ErrTray = [0.042341171  0.042996174 0.000655002];
%% Finding all tray compositions and temperatures
X = zeros(22,3); %Create matrix for tray liquid compositions (1 is reboiler)
Y = zeros(22,3); %Create matrix for tray vapor compositions (1 is reboiler)
T = zeros(22,1); %Create vector for tray temperatures (1 is reboiler)
X(1,:) = Xbottoms; %Add known bottoms composition
T(1) = FindT(X(1,1),X(1,2),X(1,3)); %Find Temperature in reboiler
Y(1,:) = findy(X(1,1),X(1,2),X(1,3),T(1)); %Find vapor composition leaving reboiler
X(2,:) = Y(1,:); %Liquid composition of tray 1 from mole balance
for i = 2:21 %Loop through all trays after reboiler
    T(i) = FindT(X(i,1),X(i,2),X(i,3)); %Find temperature of tray i
    Y(i,:) = findy(X(i,1),X(i,2),X(i,3),T(i)); %Find vapor comp of tray i
    X(i+1,1) = X(i,1) + Y(i,1) - Y(i-1,1); %Mole balance to find liquid comp of tray i+1
    X(i+1,2) = X(i,2) + Y(i,2) - Y(i-1,2);
    X(i+1,3) = X(i,3) + Y(i,3) - Y(i-1,3);
end
T(22) = FindT(X(22,1),X(22,2),X(22,3));
Y(22,:) = findy(X(22,1),X(22,2),X(22,3),T(22));
CompEfficiency = (Ydistillate-Y(1,:))./(Y(21,:)-Y(1,:)); %Vapor Efficiencies of each component
Efficiency = sum(CompEfficiency)/3; %Average vapor efficiency
fprintf("Column Vapor Efficiency is %9.8f ",Efficiency);
%% Calculate Experimental and Theoretical Composition Curves
X_E = X; Y_E = Y;
for i = 2:22
    X_E(i,:) = X_E(i-1,:) + CompEfficiency.*(X(i,:)-X(i-1,:));
    Y_E(i,:) = Y_E(i-1,:) + CompEfficiency.*(Y(i,:)-Y(i-1,:));
end
Stages = 0:1:21;
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
plot(21,Ydistillate(1),'b.','MarkerSize',25)
hold on
plot(21,Ydistillate(2),'r.','MarkerSize',25)
hold on
plot(21,Ydistillate(3),'k.','MarkerSize',25)
hold on
errorbar(21,Ydistillate(1),ErrDist(1),'b','LineWidth',2)
hold on
errorbar(21,Ydistillate(2),ErrDist(2),'r','LineWidth',2)
hold on
errorbar(21,Ydistillate(3),ErrDist(3),'k','LineWidth',2)
xlim([0 21.5])
ylim([0 1])
xlabel('Stage Number')
ylabel('Vapor Mole Fraction (y_i)')
 
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
ylim([0 1])
xlabel('Stage Number')
ylabel('Liquid Mole Fraction (x_i)')
%% Calculate Theoretical and Experimental Temperature Profiles
T_E = T;
for i = 1:length(T_E)
    T_E(i) = FindT(X_E(i,1),X_E(i,2),X_E(i,3));
end
figure('DefaultAxesFontSize',12)
plot(Stages,T-273,'--b','LineWidth',2)
hold on
plot(Stages,T_E-273,'-r','LineWidth',2)
hold on
plot(21,81.1,'r.','MarkerSize',25)
hold on
plot(11,85.2,'r.','MarkerSize',25)
hold on
plot(0,102.2,'r.','MarkerSize',25)
xlim([0 21])
xlabel('Stage Number')
ylabel('Temperature (Celsius)')
% legend({'Theoretical T Profile','Experimental T Profile'},'FontSize',12)
%% Functions called for finding values needed above
function y = findy(x1,x2,x3,T) %Calculates vapor mole fraction
    X = [x1 x2 x3];
    y = zeros(1,3);
    for i = 1:3
        y(i) = (X(i)*findgamma(i,x1,x2,x3,T)*antoine(i,T))/760;
    end
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
