close all
clear all
clc
load('Data/heading_testing.mat')
% img = imread('oldeng2.png');
% img = imrotate(img,0,'bilinear','crop');
% %imagesc([-37.80005 -37.7988], [144.9609 144.9626],img);
% imagesc([-77 70], [-39 98], flipud(img));
% 
% 
% hold on
% [trueX, trueY] = Spherical2Azimuth(Position.latitude, Position.longitude, Position.latitude(1), Position.longitude(1), 0, 0, 6371000*pi);
% 
% plot(trueX, trueY, 'b-*','linewidth',1.5);
% set(gca,'ydir','normal');
time = MagneticField.Timestamp;

% Conversion of date-time format to seconds starting from 0
initial_min = double(minute(time(1)));
initial_sec = double(second(time(1)));

t = zeros(length(time));
for i=1:length(time)
    currentmin = double(minute(time(i)));
    if currentmin == initial_min + 1
        t(i) = double(second(time(i))) + 60 - initial_sec;
    elseif currentmin == initial_min + 2
        t(i) = double(second(time(i))) + 120 - initial_sec;
    elseif currentmin == initial_min + 3
        t(i) = double(second(time(i))) + 180 - initial_sec;
    elseif currentmin == initial_min + 4
        t(i) = double(second(time(i))) + 240 - initial_sec;    
    else
        t(i) = double(second(time(i)))-initial_sec;
    end
end

Bx = MagneticField.X;
By = MagneticField.Y;
Bz = MagneticField.Z;
B = [Bx By Bz];
B = B*10^-6;

[A,b] = magcal(B);
B_cal = (B-b)*A;
Bx = B_cal(:,1);
Bx = Bx(1:100) 
By = B_cal(:,2);
By = By(1:100);
azimuth = (2*pi - atan2(By, Bx));
azimuth = unwrap(azimuth) * 180/pi;
var(azimuth)
plot(t(1:100),azimuth)