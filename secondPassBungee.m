function [ T, R ] = secondPassBungee()
%spring with air resistance

A = .07; %Cross-sectional area of human lying down (m^2)
rho = 1.2; %Density of air (kg/m^3)
Cd = 44.4; %Drag coeffiecent of human lying down (unitless)(data from Google Drive)
k = 52; % N/m
lo =30;%m
y = 70; %m
vy = 0; %m/s
g = 9.8; %m/s^2
mass1 = 60; %Mass of jumper 60 kg
restingL = lo +((-mass1*g)/k)
Data = [y, vy];
function res = changingValues(~, Data)
    y = Data(1);
    vy = Data(2);
    dydt = vy;
    if( y <= (70- restingL))
       dvydt = (-(mass1*g) - (k *(y - restingL))+(.5 *rho* Cd* A *vy^2 ))/mass1;
    else dvydt = -g;
    end
   
    res = [dydt; dvydt];
end
[T, R] = ode45(@changingValues, [0 60], Data);

Y = R(:,1);
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