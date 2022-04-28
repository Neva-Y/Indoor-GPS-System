clc
clear all
close all

file = "Data/oldeng2.mat";

[stepLength, timeStep] = distance_estimate(file);
[vectorLen, x, lat, long, gyroAzimuth, azimuth, Orientationx, t] = heading_direction(file);


longditude = zeros(vectorLen,1);
latitude = zeros(vectorLen,1);
longditude(1) = 0;
latitude(1) = 0;

j = 1;
for i = 2:vectorLen
    if i == timeStep(j) && j < length(timeStep)
        j=j+1;
        longditude(i) = longditude(i-1) + stepLength*sind(x(i));
        latitude(i) = latitude(i-1) + stepLength*cosd(x(i));
    else
        longditude(i) = longditude(i-1);
        latitude(i) = latitude(i-1);
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
% img = imread('oldeng.png');
% img = imrotate(img,0,'bilinear','crop');
% imagesc([-77 70], [-39 98], flipud(img));
% 
% 
% hold on
[trueX, trueY] = Spherical2Azimuth(lat, long, lat(1), long(1), 0, 0, 6371000*pi);

plot(trueX, trueY, 'b-*','linewidth',1.5);
hold on
plot(-latitude, -longditude,  'r-*','linewidth',1.5)

set(gca,'ydir','normal');

xlabel("Latitude (m)")
xlabel("Longitude (m)")
legend("GNSS location", "Pedestrian Dead Reckoning")
title("Walking Path")
