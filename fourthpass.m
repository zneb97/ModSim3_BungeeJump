function [ T, R ] = fourthpass()
%Spring with air resistance, tension in one direction

%universal variables
A = .5; %Cross-sectional area of human lying down (m^2)
rho = 1.2; %Density of air (kg/m^3)
y = 70; %Current height of the jumper m
Cd = 3.3; %Drag coeffiecent of human lying down (unitless)(data from Google Drive)
g = 9.81; %Acceleration due to gravity m/s^2
startH = y; %Starting height of the jumper m
%%%%%%%%%%%%%%%

% Jumper 1
k1 = 100; %Spring constant of the bungee cord N/m
lo1 =17.5;%Resting length of the cord with no mass on it m
vy1 = 0; %Current velocity of the jumper m/s
mass1 = 100; %Mass of jumper kg
restingL1 = lo1 +((mass1*g)/k1); %Resting length of cord with mass on it m
%%%%%%%%%%%%%%%
%jumper 2

k2 = 100; %Spring constant of the bungee cord N/m
lo2 =17.5;%Resting length of the cord with no mass on it m
vy2 = 0; %Current velocity of the jumper m/s
mass2 = 50; %Mass of jumper kg
restingL2 = lo2 +((mass2*g)/k2); %Resting length of cord with mass on it m

Data1 = [y, vy1 ,y, vy2];

% function kopt= jumperOptimization(dataNum)
%   changingValues(~,dataNum);
% end

    function res = changingValues1(~, Data)
        y = Data(1);
        vy = Data(2);
        dydt = vy;
        y2 = Data(3);
        vy2 = Data(4);
        dydt2 = vy2;
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
        %Acceleration
        if(y2 <= (startH - lo2)) %Below tension point
            if(vy2 >= 0) %Going up STEP 3
                dvydt2 = (-(mass2*g) + (k2 *(startH-y2-lo2)) - (.5 *rho* Cd* A *vy2^2 ))/mass2;
            else %Going down STEP 2
                dvydt2 = (-(mass2*g) + (k2 *(startH-y2-lo2)) + (.5 *rho* Cd* A *vy2^2 ))/mass2;
            end
            
        else %Above tension point
            if(vy2 >= 0) %Going up STEP 4
                dvydt2 = (-(mass2*g) -(.5 *rho* Cd* A *vy^2 ))/mass2;
            else %Going down STEP 1
                dvydt2 = (-(mass2*g) + (.5 *rho* Cd* A *vy^2 ))/mass2;
            end
        end
        
        if(dvydt2 >(3.5*9.8) || dvydt2 <-(3.5*9.8))
            disp(dvydt2);
        end
        res = [dydt; dvydt; dydt2; dvydt2];
        
    end

[T1, R1] = ode45(@changingValues1, [0 60], Data1);
hold on
figure(1);
Y1 = R1(:,1);
plot(T1,Y1);
Y2 = R1(:,3);
hold on
plot(T1,Y2);
title('Height of Bungee Jumper vs Time');
xlabel('Time (Seconds)');
ylabel('Height (Meters)');
hold on;
figure(2);
V1 = R1(:,2);
plot(T1,V1);
V2 = R1(:,4);
hold on
plot(T1,V2);
title('Velocity of Bungeer vs Time');
xlabel('Time (Seconds)');
ylabel('Velocity (m/s)');

Diff= Y1- Y2;
intersect(1) = 0;
indexOfInterest(1) = 1;
j = 1;
    function inte = findIntersect()
        for i = 2: length(T1)
            if(sign(Diff(i))~= sign(Diff(i-1)))
                intersect(i) = T1(i);
                indexOfInterest(j) = i;
                j = j+1;
            end
        end
        inte = intersect;
    end


ValuesOfIntersection = findIntersect();
Velocity1(1)=0;
Velocity2(1) = 0;
for i = 1: length(indexOfInterest)
    Velocity1(i)= Y1(i);
    Velocity2(i)= Y2(i);
end

velocityAtIntersections = [Velocity1;Velocity2]
results = length(indexOfInterest);
end