 clc
 clear all
 close all
 distance_estimate
 heading_direction

longditude = zeros(vectorLen);
latitude = zeros(vectorLen);
longditude(1) = 0;
latitude(1) = 0;

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

figure();
plot(latitude, longditude)
hold on
[trueX, trueY] = Spherical2Azimuth(Position.latitude, Position.longitude, Position.latitude(1), Position.longitude(1), 0, 0, 6371000*pi);
plot(-trueX, -trueY);
xlabel("Latitude (m)")
xlabel("Longitude (m)")
legend("Pedestrian Dead Reckoning", "GNSS location")
title("Walking Path")
