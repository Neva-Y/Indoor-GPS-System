clc
clear all
close all

file = "Data/oldeng2.mat";
measuringAlgo = "integrate";

[cumDistance, timeStep] = distance_estimate(file, measuringAlgo);
[vectorLen, x, lat, long, gyroAzimuth, azimuth, Orientationx, t, K, P] = heading_direction(file);


longditude = zeros(vectorLen,1);
latitude = zeros(vectorLen,1);
longditude(1) = 0;
latitude(1) = 0;

if measuringAlgo == "weinberg"
    j = 1;
    for i = 2:vectorLen
        if i == timeStep(j) && j < length(timeStep)
            j=j+1;
            longditude(i) = longditude(i-1) + (cumDistance(i)-cumDistance(i-1))*sind(x(i));
            latitude(i) = latitude(i-1) + (cumDistance(i)-cumDistance(i-1))*cosd(x(i));
        else
            longditude(i) = longditude(i-1);
            latitude(i) = latitude(i-1);
        end
    end

elseif measuringAlgo == "integrate"
    for i = 2:vectorLen
        longditude(i) = longditude(i-1) + (cumDistance(i)-cumDistance(i-1))*sind(x(i));
        latitude(i) = latitude(i-1) + (cumDistance(i)-cumDistance(i-1))*cosd(x(i));
    end
end
        
        
        
figure();
plot(t, gyroAzimuth, 'LineWidth', 1);
hold on 
plot(t, azimuth, 'LineWidth', 1);
plot(t, Orientationx, 'LineWidth', 1);
plot(t, x, 'k','linewidth',2);
title("Heading Direction");
xlabel('Time (sec)');
ylabel('Heading Angle (deg)');
legend('Gyro','Magnetometer','True Heading', 'Kalman Prediction')



figure();
img = imread('oldeng.png');
% %walking
%imagesc([-77 70], [-39 98], flipud(img));

%running
%imagesc([-67 42], [-34 83], flipud(img));

hold on
[trueX, trueY] = Spherical2Azimuth(lat, long, lat(1), long(1), 0, 0, 6371000*pi);
trueCoords = [trueX, trueY];
A_new = zeros(vectorLen,2) ;
xi = (1:vectorLen)' ; 
x = linspace(1,vectorLen,length(trueX))' ; 
for i = 1:2
    Ai = interp1(x,trueCoords(:,i),xi) ; 
    A_new(:,i) = Ai ; 
end


plot(A_new(:,1), A_new(:,2), 'b-*','linewidth',1.5);
hold on
latitude = -latitude;
longditude = -longditude;
plot(latitude, longditude,  'r-*','linewidth',1.5)

set(gca,'ydir','normal');

xlabel("Latitude (m)")
ylabel("Longitude (m)")
legend("GNSS location", "Pedestrian Dead Reckoning")
title("Walking Path")

MSE = immse(A_new, [latitude longditude])

