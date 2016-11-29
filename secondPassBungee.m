function [ T, R ] = secondPassBungee()
%Spring with air resistance, tension in one direction

A = .5; %Cross-sectional area of human lying down (m^2)
rho = 1.2; %Density of air (kg/m^3)
Cd = 44.4; %Drag coeffiecent of human lying down (unitless)(data from Google Drive)
k = 52; %Spring constant of the bungee cord N/m
lo =30;%Resting length of the cord with no mass on it m
y = 70; %Current height of the jumper m
startH = y; %Starting height of the jumper m
vy = 0; %Current velocity of the jumper m/s
g = 9.8; %Acceleration due to gravity m/s^2
mass1 = 60; %Mass of jumper kg
restingL = lo +((mass1*g)/k); %Resting length of cord with mass on it m

Data = [y, vy];

function res = changingValues(~, Data)
    y = Data(1);
    vy = Data(2);
    
    dydt = vy;
    if(y < (startH - lo))
       dvydt = (-(mass1*g) - (k *(startH-y-lo))-(.5 *rho* Cd* A *vy^2 ))/mass1;
    else
        dvydt = (-g*mass1 - (.5 *rho* Cd* A *vy^2 ))/mass1;
    end
    vy
    res = [dydt; dvydt];
end

[T, R] = ode45(@changingValues, [0 60], Data);

Y = R(:,1);
figure(2)
comet(T,Y);
end
