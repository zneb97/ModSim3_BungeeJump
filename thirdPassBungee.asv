function [ T, R ] = thirdPassBungee()
%Spring with air resistance, tension in one direction

%universal variables
A = .5; %Cross-sectional area of human lying down (m^2)
rho = 1.2; %Density of air (kg/m^3)
y = 70; %Current height of the jumper m
Cd = 3.3; %Drag coeffiecent of human lying down (unitless)(data from Google Drive)
g = 9.81; %Acceleration due to gravity m/s^2
startH = y; %Starting height of the jumper m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Jumper 1
k1 = 100; %Spring constant of the bungee cord N/m
lo1 =17.5;%Resting length of the cord with no mass on it m
vy1 = 0; %Current velocity of the jumper m/s
mass1 = 60; %Mass of jumper kg
restingL1 = lo1 +((mass1*g)/k1); %Resting length of cord with mass on it m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%jumper 2

k2 = 100; %Spring constant of the bungee cord N/m
lo2 =17.5;%Resting length of the cord with no mass on it m
vy2 = 0; %Current velocity of the jumper m/s
mass2 = 90; %Mass of jumper kg
restingL2 = lo2 +((mass2*g)/k2); %Resting length of cord with mass on it m

Data1 = [y, vy1 ];
Data2 = [y, vy2];

   % function kopt= jumperOptimization(dataNum)
     %   changingValues(~,dataNum);
   % end



function res = changingValues1(~, Data)
    y = Data(1);
    vy = Data(2);
    dydt = vy;
    
    %Acceleration
    if(y <= (startH - lo1)) %Below tension point
        if(vy >= 0) %Going up STEP 3
            dvydt = (-(mass1*g) + (k1 *(startH-y-lo1)) - (.5 *rho* Cd* A *vy^2 ))/mass1;
        else %Going down STEP 2
            dvydt = (-(mass1*g) + (k1 *(startH-y-lo1)) + (.5 *rho* Cd* A *vy^2 ))/mass1;
        end
        
    else %Above tension point
        if(vy >= 0) %Going up STEP 4
            dvydt = (-(mass1*g) -(.5 *rho* Cd* A *vy^2 ))/mass1;
        else %Going down STEP 1
            dvydt = (-(mass1*g) + (.5 *rho* Cd* A *vy^2 ))/mass1;
        end
    end
     if(dvydt >(3.5*9.8) || dvydt <-(3.5*9.8))
        disp(dvydt);
     end
    res = [dydt; dvydt];
    
end

function res = changingValues2(~, Data)
    y = Data(1);
    vy = Data(2);
    dydt = vy;
    
    %Acceleration
    if(y <= (startH - lo2)) %Below tension point
        if(vy >= 0) %Going up STEP 3
            dvydt = (-(mass2*g) + (k2 *(startH-y-lo2)) - (.5 *rho* Cd* A *vy^2 ))/mass2;
        else %Going down STEP 2
            dvydt = (-(mass2*g) + (k2 *(startH-y-lo2)) + (.5 *rho* Cd* A *vy^2 ))/mass2;
        end
        
    else %Above tension point
        if(vy >= 0) %Going up STEP 4
            dvydt = (-(mass2*g) -(.5 *rho* Cd* A *vy^2 ))/mass2;
        else %Going down STEP 1
            dvydt = (-(mass2*g) + (.5 *rho* Cd* A *vy^2 ))/mass2;
        end
    end
   
    if(dvydt >(3.5*9.8) || dvydt <-(3.5*9.8))
        disp(dvydt);
    end
    res = [dydt; dvydt];
    
end
%for k=5:5:50
 %   for lo = 10:5:50
[T1, R1] = ode45(@changingValues1, [0 60], Data1);
[T2, R2] = ode45(@changingValues2, [0 60], Data2);
hold on
figure(1);
Y1 = R1(:,1);
plot(T1,Y1);
Y2 = R2(:,1);
hold on
plot(T2,Y2);
title('Height of Bungee Jumper vs Time');
xlabel('Time (Seconds)');
ylabel('Height (Meters)');

figure(2);
V1 = R1(:,2);
plot(T1,V1);
V2 = R2(:,2);
%f=fit(V1,V2,'poly2')
hold on;
%plot(f,V1,V2) 
hold on
plot(T2,V2);
title('Velocity of Bungeer vs Time');
xlabel('Time (Seconds)');
ylabel('Velocity (m/s)');

figure(3)
t = (1:50)';
 X = ones(50,3);
 X(:,2) = cos((2*pi)/50*t);
 X(:,3) = sin((2*pi)/50*t);
 y = 2*cos((2*pi)/50*t-pi/4)+randn(size(t));
 y = y(:);
 beta = X\y;
 yhat = beta(1)+beta(2)*cos((2*pi)/50*t)+beta(3)*sin((2*pi)/50*t);
 plot(t,y,'b');
 hold on
 plot(t,yhat,'r','linewidth',2);
end
