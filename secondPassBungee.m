function [ T, R ] = secondPassBungee()
%Spring with air resistance, tension in one direction

A = .7; %Cross-sectional area of human lying down (m^2)
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
    if
       dvydt = (-(mass1*g) - (k *(y - restingL))+(.5 *rho* Cd* A *vy^2 ))/mass1;
    else
        dvydt = (-g*mass1 +(.5 *rho* Cd* A *vy^2 ))/mass1;;
    end
    end
   
    res = [dydt; dvydt];
end
[T, R] = ode45(@changingValues, [0 60], Data);

Y = R(:,1);
figure(2)
comet(T,Y);
end
%Constants
% 
%
% 
% 
% A = 0; %Cross-sectional area of human lying down (m^2)
% rho = 1.2; %Density of air (kg/m^3)
% g = 9.81; %Gravity acceleration (m/sec^2)
% Cd = 44.4; %Drag coeffiecent of human lying down (unitless)(data from Google Drive)
% h = 91.44; %Height of bridge (m) (Sacramento California)
% 
% 
% 
% 
% L0 = 45
% k = 52.2; %Starting k value for bungee (N/kg)