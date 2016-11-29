function [ T, R ] = secondPassBungee()
%Spring with air resistance, tension in one direction

A = .5; %Cross-sectional area of human lying down (m^2)
rho = 1.2; %Density of air (kg/m^3)
Cd = 44.4; %Drag coeffiecent of human lying down (unitless)(data from Google Drive)
k = 100; %Spring constant of the bungee cord N/m
lo =17.5;%Resting length of the cord with no mass on it m
y = 70; %Current height of the jumper m
startH = y; %Starting height of the jumper m
vy = 0; %Current velocity of the jumper m/s
g = 9.81; %Acceleration due to gravity m/s^2
mass1 = 60; %Mass of jumper kg
restingL = lo +((mass1*g)/k); %Resting length of cord with mass on it m

Data = [y, vy];

function res = changingValues(~, Data)
    y = Data(1);
    vy = Data(2);
    
    dydt = vy;
    
    %Acceleration
    if(y <= (startH - lo)) %Below tension point
        if(vy >= 0) %Going up STEP 3
            dvydt = (-(mass1*g) + (k *(startH-y-lo)) - (.5 *rho* Cd* A *vy^2 ))/mass1;
        else %Going down STEP 2
            dvydt = (-(mass1*g) + (k *(startH-y-lo)) + (.5 *rho* Cd* A *vy^2 ))/mass1;
        end
        
    else %Above tension point
        if(vy >= 0) %Going up STEP 4
            dvydt = (-(mass1*g) -(.5 *rho* Cd* A *vy^2 ))/mass1;
        else %Going down STEP 1
            dvydt = (-(mass1*g) + (.5 *rho* Cd* A *vy^2 ))/mass1;
        end
    end
    
    res = [dydt; dvydt];
    
end
for k=5:5:50
    for lo = 10:5:50
[T, R] = ode45(@changingValues, [0 60], Data);
hold on
Y = R(:,1);
plot(T,Y);
title('Height of Bungee Jumper vs Time');
xlabel('Time (Seconds)');
ylabel('Height (Meters)');
    end
end
end
